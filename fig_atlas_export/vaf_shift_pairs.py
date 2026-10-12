"""샘플 쌍 사이 VAF 변화의 모양: 클론 구성 변화(양방향으로 퍼짐)인지 측정 편향(한 방향 이동)인지.

사이트: 품질 맞춘 재호출(chr1–3)에서 12샘플 중 하나라도 PASS 이형접합 biallelic SNV로 호출한 자리
read 수준 AD: 12개 품질 맞춘 BAM을 한 번에 bcftools mpileup (-q 1 -Q 13)
쌍마다 두 샘플 모두 depth ≥ 15이고 합친 VAF가 0.05–0.95인 사이트에서
  z = (VAF_y - VAF_x) / sqrt(p(1-p)(1/n_x + 1/n_y)),  p = 합친 VAF
  mean_z  : 한 방향 이동(측정 편향이면 0에서 벗어남)
  var_z   : 이항 표집으로 기대되는 1보다 큰 만큼이 실제 VAF 변화(구성 변화면 커짐)
  up3/down3 : z > 3, z < -3 사이트 비율
쌍 유형: within-run 다른 조건(예: A–D), across-run 같은 조건(예: A1–A2, 웰 사이 기본 변동), across-run 다른 조건
부분집합: all, nonrepeat(RepeatMasker·simpleRepeat 밖)
출력: OUT.sites.tsv.gz(사이트 × 샘플 alt/depth), OUT.pairs.tsv, stdout 요약
사용: python vaf_shift_pairs.py QM_DIR OUT_PREFIX
"""
import gzip
import itertools
import subprocess
import sys

import numpy as np
import pandas as pd

BCF = "/home/ygkim/program/anaconda3/envs/wgs-tertiary-pipeline/bin/bcftools"
RES = "/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39"
REPO = "/data_4tb/hifi-human-wgs-wdl-custom"
QM, OUT = sys.argv[1], sys.argv[2]
CHR = ["chr1", "chr2", "chr3"]
REF = f"{RES}/mouse_GRCm39.fasta"
meta = pd.read_csv(f"{REPO}/base_mod_per_sample/samples.tsv", sep="\t")
SIDS = meta.sample_id.tolist()
GRP = dict(zip(meta.sample_id, meta.group))
RUN = dict(zip(meta.sample_id, meta.run))

# 1) 사이트 합집합
sites = set()
for s in SIDS:
    res = subprocess.run([BCF, "query", "-t", ",".join(CHR), "-i", 'FILTER="PASS" && TYPE="snp"', "-f",
                          "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", f"{QM}/dv/{s}/{s}.vcf.gz"],
                         capture_output=True, text=True, check=True).stdout
    for line in res.splitlines():
        c, p, ref, alt, gt = line.split("\t")
        al = gt.replace("|", "/").split("/")
        if "," not in alt and len(al) == 2 and al[0] != al[1]:
            sites.add((c, int(p), ref, alt))
sites = sorted(sites)
print(f"union het SNV sites: {len(sites)}", flush=True)
reg = f"{OUT}.regions.tsv"
pd.DataFrame([(c, p) for c, p, _, _ in sites]).drop_duplicates().to_csv(reg, sep="\t", index=False, header=False)

# 2) BAM마다 mpileup을 병렬로 (12개 동시 mpileup은 느림)
want = {(c, p): a for c, p, _, a in sites}
procs = {}
for s in SIDS:
    cmd = (f"{BCF} mpileup -f {REF} -R {reg} -a AD -d 10000 -q 1 -Q 13 -Ou {QM}/bams/{s}.qm.chr1-3.bam "
           f"| {BCF} query -f '%CHROM\\t%POS\\t%ALT[\\t%AD]\\n' > {OUT}.{s}.ad.tsv")
    procs[s] = subprocess.Popen(["bash", "-c", cmd])
for s, pr in procs.items():
    if pr.wait() != 0:
        sys.exit(f"mpileup failed: {s}")
keys = [(c, p) for c, p, _, _ in sites]
keys = list(dict.fromkeys(keys))
alt_n, dep = {}, {}
for s in SIDS:
    av, dv = {}, {}
    for line in open(f"{OUT}.{s}.ad.tsv"):
        f = line.rstrip("\n").split("\t")
        k = (f[0], int(f[1]))
        if k not in want:
            continue
        alts = f[2].split(",")
        v = [int(x) if x != "." else 0 for x in f[3].split(",")]
        av[k] = v[1 + alts.index(want[k])] if want[k] in alts else 0
        dv[k] = sum(v)
    alt_n[s] = [av.get(k, 0) for k in keys]
    dep[s] = [dv.get(k, 0) for k in keys]
    subprocess.run(["rm", "-f", f"{OUT}.{s}.ad.tsv"])
subprocess.run(["rm", "-f", reg])
A = pd.DataFrame(alt_n, index=pd.MultiIndex.from_tuples(keys, names=["chrom", "pos"]))
D = pd.DataFrame(dep, index=A.index)

# 반복서열 마스크
fai = pd.read_csv(f"{REF}.fai", sep="\t", header=None, usecols=[0, 1], names=["c", "n"]).set_index("c")
mask = {c: np.zeros(int(fai.loc[c, "n"]) + 1, bool) for c in CHR}
for path, cols in ((f"{RES}/rmsk.mm39.txt.gz", (5, 6, 7)), (f"{RES}/simpleRepeat.txt.gz", (1, 2, 3))):
    with gzip.open(path, "rt") as fh:
        for line in fh:
            x = line.split("\t")
            if x[cols[0]] in mask:
                mask[x[cols[0]]][int(x[cols[1]]):int(x[cols[2]])] = True
rep = np.array([mask[c][p - 1] for c, p in A.index])
out = pd.concat([A.add_suffix(".alt"), D.add_suffix(".dp")], axis=1)
out["in_repeat"] = rep
out.to_csv(f"{OUT}.sites.tsv.gz", sep="\t")

# 3) 쌍별 z
rows = []
for x, y in itertools.combinations(SIDS, 2):
    if RUN[x] == RUN[y]:
        kind = "within_run_diff_cond"
    elif GRP[x] == GRP[y]:
        kind = "across_run_same_cond"
    else:
        kind = "across_run_diff_cond"
    # 방향을 조건 순서 A<B<C<D로 고정(y가 뒤 조건): z>0이면 뒤 조건에서 VAF가 높음
    if GRP[x] > GRP[y] or (GRP[x] == GRP[y] and RUN[x] > RUN[y]):
        x, y = y, x
    ax, ay, nx, ny = A[x].to_numpy(float), A[y].to_numpy(float), D[x].to_numpy(float), D[y].to_numpy(float)
    p = (ax + ay) / np.maximum(nx + ny, 1)
    ok = (nx >= 15) & (ny >= 15) & (p >= 0.05) & (p <= 0.95)
    for sub, m in (("all", ok), ("nonrepeat", ok & ~rep)):
        z = (ay[m] / ny[m] - ax[m] / nx[m]) / np.sqrt(p[m] * (1 - p[m]) * (1 / nx[m] + 1 / ny[m]))
        rows.append(dict(x=x, y=y, gx=GRP[x], gy=GRP[y], rx=RUN[x], ry=RUN[y], kind=kind, subset=sub, n=int(m.sum()),
                         mean_z=z.mean(), var_z=z.var(ddof=1), up3=(z > 3).mean(), down3=(z < -3).mean(),
                         mean_dvaf=float(np.mean(ay[m] / ny[m] - ax[m] / nx[m]))))
P = pd.DataFrame(rows)
P.to_csv(f"{OUT}.pairs.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
for sub in ("all", "nonrepeat"):
    s = P[P.subset == sub].copy()
    s["pair"] = s.gx + s.rx.astype(str) + "-" + s.gy + s.ry.astype(str)
    print(f"\n## {sub}: pair summary")
    print(s.groupby("kind")[["n", "mean_z", "var_z", "up3", "down3"]].mean().round(4).to_string())
    print("\n# within-run pairs")
    print(s[s.kind == "within_run_diff_cond"][["pair", "n", "mean_z", "var_z", "up3", "down3", "mean_dvaf"]]
          .round(4).to_string(index=False))
    print("\n# across-run same-condition pairs")
    print(s[s.kind == "across_run_same_cond"][["pair", "n", "mean_z", "var_z", "up3", "down3", "mean_dvaf"]]
          .round(4).to_string(index=False))
    t = s.groupby(["gx", "gy"]).var_z.mean().unstack()
    print("\n# mean var_z by condition pair (rows x, cols y)")
    print(t.round(3).to_string())
