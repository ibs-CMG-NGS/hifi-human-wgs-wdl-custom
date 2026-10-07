"""Fig 3 클론성 지표의 견고성 점검 (indel·시퀀싱 품질 교란 확인).

샘플마다 phased small-variant VCF(상염색체, PASS)를 한 번 훑어 다음을 계산한다.
  region = all | nonrepeat (UCSC mm39 RepeatMasker + simpleRepeat 밖)
  - het_snv, het_del, het_ins, het_hom_peak_frac_snv (SNV VAF > 0.9 비율)
  - het SNV VAF 구간 비율: <0.2, 0.2–0.35, 0.35–0.65, 0.65–0.9, >0.9
  - haplotype 방향 불균형 (wgs-tertiary-analysis scripts/haplotype_imbalance.py와 같은 규칙:
    GT 0|1 / 1|0, PS 블록 안 사이트 ≥ 10개, phase-0 allele VAF > 0.5 부호의 양측 이항검정 p < 0.05)
    를 variant 집합 snv / all(SNV+indel) 두 가지로
출력: OUT_DIR/<sample>.robust.tsv (sample_id, region, metric, value)
사용: python clonality_robustness.py samples.tsv OUT_DIR
"""
import gzip
import math
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

BCF = "/home/ygkim/program/anaconda3/envs/wgs-tertiary-pipeline/bin/bcftools"
RES = "/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39"
AUTOS = [f"chr{i}" for i in range(1, 20)]
meta = pd.read_csv(sys.argv[1], sep="\t")
outdir = sys.argv[2]

fai = pd.read_csv(f"{RES}/mouse_GRCm39.fasta.fai", sep="\t", header=None, usecols=[0, 1], names=["c", "n"]).set_index("c").n
mask = {c: np.zeros(int(fai[c]) + 1, bool) for c in AUTOS}
for path, cols in ((f"{RES}/rmsk.mm39.txt.gz", (5, 6, 7)), (f"{RES}/simpleRepeat.txt.gz", (1, 2, 3))):
    with gzip.open(path, "rt") as f:
        for line in f:
            x = line.split("\t")
            c = x[cols[0]]
            if c in mask:
                mask[c][int(x[cols[1]]):int(x[cols[2]])] = True
print("repeat mask built:", {c: round(mask[c].mean(), 3) for c in AUTOS[:3]}, flush=True)


def binom_two_sided(k, n):
    if n == 0:
        return 1.0
    logp = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + n * math.log(0.5) for i in range(n + 1)]
    pk = logp[k]
    return min(1.0, sum(math.exp(lp) for lp in logp if lp <= pk + 1e-9))


for _, r in meta.iterrows():
    o = r.wgs_haplotagged_bam.split("/out/")[0] + "/out"
    vcf = f"{o}/phased_small_variant_vcf/{r.sample_id}.GRCm39.small_variants.phased.vcf.gz"
    cmd = [BCF, "query", "-t", ",".join(AUTOS), "-i", 'FILTER="PASS"', "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%PS\t%VAF]\n", vcf]
    stats = {reg: defaultdict(float) for reg in ("all", "nonrepeat")}
    vafs = {reg: [] for reg in ("all", "nonrepeat")}
    blocks = {(reg, vs): defaultdict(list) for reg in ("all", "nonrepeat") for vs in ("snv", "all")}
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as p:
        for line in p.stdout:
            c, pos, ref, alt, gt, ps, vaf = line.rstrip("\n").split("\t")
            al = gt.replace("|", "/").split("/")
            if len(al) != 2 or al[0] == al[1] or "." in al:
                continue
            try:
                v = max(float(x) for x in vaf.split(",") if x != ".")
            except ValueError:
                continue
            alts = alt.split(",")
            snv = len(ref) == 1 and all(len(a) == 1 for a in alts)
            regs = ["all"] if mask[c][int(pos) - 1] else ["all", "nonrepeat"]
            for reg in regs:
                st = stats[reg]
                if snv:
                    st["het_snv"] += 1
                    vafs[reg].append(v)
                elif len(ref) > max(len(a) for a in alts):
                    st["het_del"] += 1
                else:
                    st["het_ins"] += 1
                if gt in ("0|1", "1|0") and ps not in (".", ""):
                    ph0 = v if gt == "1|0" else 1 - v
                    blocks[(reg, "all")][ps].append(ph0)
                    if snv:
                        blocks[(reg, "snv")][ps].append(ph0)
    rows = []
    for reg in ("all", "nonrepeat"):
        st = stats[reg]
        va = np.array(vafs[reg])
        for k in ("het_snv", "het_del", "het_ins"):
            rows.append((reg, k, int(st[k])))
        rows.append((reg, "het_hom_peak_frac_snv", float((va > 0.9).mean())))
        for lo, hi, nm in ((0, 0.2, "lt0.2"), (0.2, 0.35, "0.2-0.35"), (0.35, 0.65, "0.35-0.65"), (0.65, 0.9, "0.65-0.9")):
            rows.append((reg, f"het_snv_vaf_frac_{nm}", float(((va >= lo) & (va < hi)).mean())))
        rows.append((reg, "het_snv_vaf_frac_gt0.9", float((va > 0.9).mean())))
        for vs in ("snv", "all"):
            tested = biased = 0
            for fr in blocks[(reg, vs)].values():
                if len(fr) < 10:
                    continue
                tested += 1
                k = sum(1 for f in fr if f > 0.5)
                n = sum(1 for f in fr if f != 0.5)
                if n and binom_two_sided(k, n) < 0.05:
                    biased += 1
            rows += [(reg, f"hap_blocks_tested_{vs}", tested), (reg, f"hap_blocks_biased_{vs}", biased),
                     (reg, f"hap_imbalance_block_frac_{vs}", biased / tested if tested else float("nan"))]
    pd.DataFrame(rows, columns=["region", "metric", "value"]).assign(sample_id=r.sample_id).to_csv(
        f"{outdir}/{r.sample_id}.robust.tsv", sep="\t", index=False)
    print(r.sample_id, {m: v for reg, m, v in rows if reg == "all" and m in ("het_snv", "het_del", "hap_imbalance_block_frac_snv", "hap_imbalance_block_frac_all")}, flush=True)
