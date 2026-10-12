"""무처치(A)의 이형접합 SNV 과잉이 실제 유전 구성(개체 유래 변이)인지 라이브러리 artifact인지 가른다.

품질 맞춘 재호출(chr1–3)에서 run마다
  A_only = A에서 PASS 이형접합 SNV로 호출되고 같은 run의 D에서는 어떤 유전형으로도 호출되지 않은 SNV
  D_only = 반대 방향(대조 집합; 두 집합의 크기 차이가 과잉)
각 집합에 대해
  in_mgp      : Mouse Genomes Project REL2021 SNP(같은 위치·같은 allele)에 있는가 → 알려진 계통 다형성
  rep_other_A : 다른 두 run의 같은 조건 샘플(A_only면 A, D_only면 D) 중 몇 개에서도 호출되었는가 (0–2)
  in_BC       : 같은 run의 B 또는 C에서 호출되었는가
  other_alt / other_dp / other_vaf : 반대쪽 샘플(A_only면 D) 품질 맞춘 BAM의 read 수준 alt 근거
  own_vaf     : 호출된 샘플의 VAF
  in_repeat   : RepeatMasker 또는 simpleRepeat 안,  near_indel : 같은 샘플의 PASS indel에서 10 bp 이내
  sbs         : 피리미딘 기준 치환 유형
출력: OUT.sites.tsv.gz(사이트별), OUT.summary.tsv(run × 집합 요약), stdout
사용: python snv_excess_origin.py QM_DIR MGP_BCF OUT_PREFIX
"""
import gzip
import subprocess
import sys

import numpy as np
import pandas as pd

BCF = "/home/ygkim/program/anaconda3/envs/wgs-tertiary-pipeline/bin/bcftools"
RES = "/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39"
REPO = "/data_4tb/hifi-human-wgs-wdl-custom"
QM, MGP, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
CHR = ["chr1", "chr2", "chr3"]
REF = f"{RES}/mouse_GRCm39.fasta"
meta = pd.read_csv(f"{REPO}/base_mod_per_sample/samples.tsv", sep="\t")
SID = {(r.run, r.group): r.sample_id for r in meta.itertuples()}
COMP = str.maketrans("ACGT", "TGCA")

mask = {}
fai = pd.read_csv(f"{REF}.fai", sep="\t", header=None, usecols=[0, 1], names=["c", "n"]).set_index("c")
for c in CHR:
    mask[c] = np.zeros(int(fai.loc[c, "n"]) + 1, bool)
for path, cols in ((f"{RES}/rmsk.mm39.txt.gz", (5, 6, 7)), (f"{RES}/simpleRepeat.txt.gz", (1, 2, 3))):
    with gzip.open(path, "rt") as f:
        for line in f:
            x = line.split("\t")
            if x[cols[0]] in mask:
                mask[x[cols[0]]][int(x[cols[1]]):int(x[cols[2]])] = True


def load(sid):
    """PASS 호출: biallelic SNV(유전형, VAF)와 indel 위치."""
    res = subprocess.run([BCF, "query", "-t", ",".join(CHR), "-i", 'FILTER="PASS"', "-f",
                          "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%VAF]\n", f"{QM}/dv/{sid}/{sid}.vcf.gz"],
                         capture_output=True, text=True, check=True).stdout
    snv, indel = {}, {c: [] for c in CHR}
    for line in res.splitlines():
        c, p, ref, alt, gt, vaf = line.split("\t")
        al = gt.replace("|", "/").split("/")
        if "1" not in al and "2" not in al:
            continue
        if len(ref) == 1 and len(alt) == 1:
            snv[(c, int(p), ref, alt)] = ("het" if al[0] != al[1] else "hom", float(vaf) if vaf != "." else np.nan)
        else:
            indel[c].append(int(p))
    return snv, {c: np.array(sorted(v)) for c, v in indel.items()}


def near(indel, c, p, w=10):
    a = indel[c]
    i = np.searchsorted(a, p)
    return any(0 <= j < len(a) and abs(a[j] - p) <= w for j in (i - 1, i))


def pileup_ad(sid, sites):
    """반대쪽 샘플 BAM에서 사이트별 (alt read, depth)."""
    tsv = f"{OUT}.{sid}.regions.tsv"
    pd.DataFrame(sorted(sites), columns=["c", "p", "r", "a"])[["c", "p"]].to_csv(tsv, sep="\t", index=False, header=False)
    cmd = (f"{BCF} mpileup -f {REF} -R {tsv} -a AD -d 10000 -q 1 -Q 13 -Ou {QM}/bams/{sid}.qm.chr1-3.bam "
           f"| {BCF} query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n'")
    out = {}
    for line in subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, check=True).stdout.splitlines():
        c, p, ref, alts, ad = line.split("\t")
        alts, ad = alts.split(","), [int(x) for x in ad.split(",")]
        out[(c, int(p))] = (dict(zip(alts, ad[1:])), sum(ad))
    subprocess.run(["rm", "-f", tsv])
    res = {}
    for c, p, r, a in sites:
        d, dp = out.get((c, p), ({}, 0))
        res[(c, p, r, a)] = (d.get(a, 0), dp)
    return res


def in_mgp(sites):
    if MGP == "none":
        return None
    tsv = f"{OUT}.mgp_regions.tsv"
    pd.DataFrame(sorted({(c, p) for c, p, _, _ in sites}), columns=["c", "p"]).to_csv(tsv, sep="\t", index=False, header=False)
    res = subprocess.run([BCF, "query", "-R", tsv, "-f", "%CHROM\t%POS\t%REF\t%ALT\n", MGP],
                         capture_output=True, text=True, check=True).stdout
    subprocess.run(["rm", "-f", tsv])
    known = set()
    for line in res.splitlines():
        c, p, ref, alts = line.split("\t")
        for a in alts.split(","):
            known.add((c, int(p), ref, a))
    return known


calls = {s: load(s) for s in meta.sample_id}
rows = []
for run in (1, 2, 3):
    for name, own, other in (("A_only", "A", "D"), ("D_only", "D", "A")):
        so, sx = SID[(run, own)], SID[(run, other)]
        snv_o, indel_o = calls[so]
        snv_x = calls[sx][0]
        sites = [k for k, (gt, _) in snv_o.items() if gt == "het" and k not in snv_x]
        rep_ids = [SID[(r, own)] for r in (1, 2, 3) if r != run]
        bc = [SID[(run, g)] for g in "BC"]
        ad = pileup_ad(sx, sites)
        for k in sites:
            c, p, r, a = k
            alt_n, dp = ad[k]
            r2, a2 = (r, a) if r in "CT" else (r.translate(COMP), a.translate(COMP))
            rows.append(dict(run=run, set=name, sample=so, other=sx, chrom=c, pos=p, ref=r, alt=a,
                             sbs=f"{r2}>{a2}", own_vaf=snv_o[k][1],
                             rep_other_runs=sum(k in calls[s][0] for s in rep_ids),
                             in_BC=any(k in calls[s][0] for s in bc),
                             other_alt=alt_n, other_dp=dp, other_vaf=alt_n / dp if dp else np.nan,
                             in_repeat=bool(mask[c][p - 1]), near_indel=near(indel_o, c, p)))
d = pd.DataFrame(rows)
known = in_mgp(set(zip(d.chrom, d.pos, d.ref, d.alt)))
d["in_mgp"] = [k in known for k in zip(d.chrom, d.pos, d.ref, d.alt)] if known is not None else False
d.to_csv(f"{OUT}.sites.tsv.gz", sep="\t", index=False)

# 같은 run에서 A와 D 모두 이형접합인 SNV(공유 배경)의 MGP 비율: 기준선
base = []
for run in (1, 2, 3):
    a, dd = calls[SID[(run, "A")]][0], calls[SID[(run, "D")]][0]
    base += [k for k, (gt, _) in a.items() if gt == "het" and dd.get(k, ("", 0))[0] == "het"]
kb = in_mgp(set(base))
if kb is not None:
    print(f"baseline (het in both A and D of the same run): n={len(base)}, in_mgp={np.mean([k in kb for k in base]):.3f}")

summ = d.groupby(["run", "set"]).agg(
    n=("pos", "size"), in_mgp=("in_mgp", "mean"), rep_ge1=("rep_other_runs", lambda x: (x >= 1).mean()),
    rep_2=("rep_other_runs", lambda x: (x == 2).mean()), in_BC=("in_BC", "mean"),
    own_vaf_median=("own_vaf", "median"), other_alt_ge2=("other_alt", lambda x: (x >= 2).mean()),
    other_vaf_median=("other_vaf", "median"), other_dp_median=("other_dp", "median"),
    in_repeat=("in_repeat", "mean"), near_indel=("near_indel", "mean")).reset_index()
summ.to_csv(f"{OUT}.summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(summ.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print("\n## SBS fraction by set (all runs)")
print(pd.crosstab(d.sbs, d.set, normalize="columns").round(3).to_string())
print("\n## A_only: MGP vs replication across runs")
print(pd.crosstab(d[d.set == "A_only"].rep_other_runs, d[d.set == "A_only"].in_mgp, margins=True).to_string())
print("\n## other-sample alt read evidence (other_alt bins)")
d["oa_bin"] = pd.cut(d.other_alt, [-1, 0, 1, 2, 4, 1000], labels=["0", "1", "2", "3-4", "5+"])
print(pd.crosstab([d.set, d.in_mgp], d.oa_bin, normalize="index").round(3).to_string())
