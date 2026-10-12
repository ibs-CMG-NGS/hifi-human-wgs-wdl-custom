"""품질 맞춘 재호출(chr1–3)에서 이형접합 결실의 성격 비교: 시퀀싱/DNA artifact인지 유전적 차이인지.

샘플마다 PASS 이형접합 결실(REF가 ALT보다 긴 biallelic indel)에 대해
  length = len(REF) - len(ALT)
  hp_run = 결실 위치에서 시작하는 참조 homopolymer 길이 (삭제된 첫 염기와 같은 염기가 연속된 수)
  in_hp4 = hp_run >= 4,  in_repeat = RepeatMasker 또는 simpleRepeat 안
  vaf    = FORMAT/VAF
  shared = 같은 run의 다른 3개 샘플 중 하나라도 같은 자리·같은 allele을 호출(PASS, 유전형 무관)
출력: 샘플별 요약(stdout + OUT.tsv), 길이·homopolymer·VAF 분포
사용: python deletion_character.py QM_DIR RUN OUT_PREFIX
"""
import gzip
import subprocess
import sys

import numpy as np
import pandas as pd

BCF = "/home/ygkim/program/anaconda3/envs/wgs-tertiary-pipeline/bin/bcftools"
RES = "/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39"
REPO = "/data_4tb/hifi-human-wgs-wdl-custom"
QM, RUNSEL, OUT = sys.argv[1], int(sys.argv[2]), sys.argv[3]
CHR = ["chr1", "chr2", "chr3"]
meta = pd.read_csv(f"{REPO}/base_mod_per_sample/samples.tsv", sep="\t")
meta = meta[meta.run == RUNSEL]

fai = pd.read_csv(f"{RES}/mouse_GRCm39.fasta.fai", sep="\t", header=None, usecols=[0, 1, 2, 3, 4],
                  names=["c", "n", "off", "lb", "lw"]).set_index("c")
seq, mask = {}, {}
for c in CHR:
    r = fai.loc[c]
    with open(f"{RES}/mouse_GRCm39.fasta", "rb") as f:
        f.seek(r.off)
        raw = f.read(int(r.n + (r.n + r.lb - 1) // r.lb * (r.lw - r.lb)))
    seq[c] = raw.replace(b"\n", b"").upper()[: r.n]
    mask[c] = np.zeros(int(r.n) + 1, bool)
for path, cols in ((f"{RES}/rmsk.mm39.txt.gz", (5, 6, 7)), (f"{RES}/simpleRepeat.txt.gz", (1, 2, 3))):
    with gzip.open(path, "rt") as f:
        for line in f:
            x = line.split("\t")
            if x[cols[0]] in mask:
                mask[x[cols[0]]][int(x[cols[1]]):int(x[cols[2]])] = True


def calls(sid):
    res = subprocess.run([BCF, "query", "-t", ",".join(CHR), "-i", 'FILTER="PASS"', "-f",
                          "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%VAF]\n", f"{QM}/dv/{sid}/{sid}.vcf.gz"],
                         capture_output=True, text=True, check=True).stdout
    rows = []
    for line in res.splitlines():
        c, p, ref, alt, gt, vaf = line.split("\t")
        rows.append((c, int(p), ref, alt, gt, vaf))
    return rows


allcalls = {s: calls(s) for s in meta.sample_id}
keys = {s: {(c, p, ref, alt) for c, p, ref, alt, _, _ in rows} for s, rows in allcalls.items()}
summ, dist = [], []
for s, rows in allcalls.items():
    others = set().union(*[keys[o] for o in keys if o != s])
    recs = []
    for c, p, ref, alt, gt, vaf in rows:
        if "," in alt or len(ref) <= len(alt):
            continue
        al = gt.replace("|", "/").split("/")
        if len(al) != 2 or al[0] == al[1] or "." in al:
            continue
        L = len(ref) - len(alt)
        p0 = p  # VCF POS(1-based)의 다음 염기(0-based p)부터 삭제
        b = seq[c][p0:p0 + 1]
        hp = 0
        while p0 + hp < len(seq[c]) and seq[c][p0 + hp:p0 + hp + 1] == b:
            hp += 1
        try:
            v = float(vaf)
        except ValueError:
            v = np.nan
        recs.append(dict(length=L, hp_run=hp, in_hp4=hp >= 4, in_repeat=bool(mask[c][p0]), vaf=v,
                         shared=(c, p, ref, alt) in others))
    d = pd.DataFrame(recs)
    grp = meta.set_index("sample_id").loc[s, "group"]
    summ.append(dict(sample_id=s, condition=grp, n_het_del=len(d), frac_1bp=(d.length == 1).mean(),
                     frac_hp4=d.in_hp4.mean(), frac_repeat=d.in_repeat.mean(), frac_shared=d.shared.mean(),
                     vaf_median=d.vaf.median(), frac_vaf_lt0_3=(d.vaf < 0.3).mean(),
                     n_private=int((~d.shared).sum()),
                     private_frac_1bp=(d[~d.shared].length == 1).mean(), private_frac_hp4=d[~d.shared].in_hp4.mean(),
                     private_vaf_median=d[~d.shared].vaf.median(),
                     shared_vaf_median=d[d.shared].vaf.median(), shared_frac_hp4=d[d.shared].in_hp4.mean()))
    for col, bins in (("length", [1, 2, 3, 4, 6, 11, 51, 10**6]), ("hp_run", [0, 1, 2, 4, 6, 8, 11, 16, 10**6]),
                      ("vaf", [0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.01])):
        h = pd.cut(d[col], bins, right=False).value_counts(sort=False)
        for k, n in h.items():
            dist.append(dict(sample_id=s, condition=grp, variable=col, bin=str(k), n=int(n)))
S = pd.DataFrame(summ)
S.to_csv(f"{OUT}.summary.tsv", sep="\t", index=False)
pd.DataFrame(dist).to_csv(f"{OUT}.dist.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(S.sort_values("condition").to_string(index=False, float_format=lambda x: f"{x:.3f}"))
dd = pd.DataFrame(dist)
for var in ("length", "hp_run", "vaf"):
    t = dd[dd.variable == var].pivot_table(index="bin", columns="condition", values="n", sort=False)
    print(f"\n## {var} distribution (counts)\n" + t.to_string())
