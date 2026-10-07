"""샘플 하나의 phased small-variant VCF / SV VCF를 한 번 훑어 fig-atlas용 샘플 지표와 VAF 분포를 만든다.

상염색체(chr1–19)만, FILTER=PASS.
- snv_pass_autosome / indel_pass_autosome: PASS 레코드를 SNV(REF·ALT 모두 1염기)와 그 외로 분류
- het_snv_autosome: GT가 서로 다른 두 allele인 PASS SNV
- het_hom_peak_frac: het SNV 중 FORMAT/VAF(멀티알렐릭은 최댓값) > 0.9 비율
  (wgs-tertiary-analysis scripts/sample_level_summary.py와 같은 정의, 단 chrX/Y 외 contig도 제외)
- het_hom_peak_frac_ds: 같은 het SNV 집합에서 AD를 표본 평균 depth → 25×로 이항 thinning
  (ref, alt 각각 Binomial(n, 25/mean_depth); seed 고정) 후 VAF > 0.9 비율. 유전형 재호출은 하지 않음
- VAF 분포: het SNV의 FORMAT/VAF를 0.01 bin, depth bin(<20, 20-30, 30-40, 40+)별 개수
- sv_pass_autosome: SV VCF PASS 레코드 수
출력: OUT_PREFIX.metrics.tsv (metric, value), OUT_PREFIX.vaf_bins.tsv
사용: python vcf_sample_pass.py SAMPLE SMALL_VCF SV_VCF MEAN_DEPTH OUT_PREFIX
"""
import subprocess
import sys

import numpy as np

BCF = "/home/ygkim/program/anaconda3/envs/wgs-tertiary-pipeline/bin/bcftools"
AUTOS = ",".join(f"chr{i}" for i in range(1, 20))
TARGET_DEPTH = 25.0
sid, vcf, sv_vcf, mean_depth, out = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
rng = np.random.default_rng(20261008)
f_thin = min(1.0, TARGET_DEPTH / mean_depth)

n_snv = n_indel = 0
het_vaf, het_ref, het_alt, het_dp = [], [], [], []
cmd = [BCF, "query", "-t", AUTOS, "-i", 'FILTER="PASS"', "-f", "%REF\t%ALT[\t%GT\t%VAF\t%AD\t%DP]\n", vcf]
with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as p:
    for line in p.stdout:
        ref, alt, gt, vaf, ad, dp = line.rstrip("\n").split("\t")
        alts = alt.split(",")
        is_snv = len(ref) == 1 and all(len(a) == 1 for a in alts)
        if is_snv:
            n_snv += 1
        else:
            n_indel += 1
            continue
        al = gt.replace("|", "/").split("/")
        if len(al) != 2 or al[0] == al[1] or "." in al:
            continue
        try:
            v = max(float(x) for x in vaf.split(",") if x != ".")
        except ValueError:
            continue
        ads = [int(x) if x != "." else 0 for x in ad.split(",")]
        het_vaf.append(v)
        het_ref.append(ads[0])
        het_alt.append(max(ads[1:]) if len(ads) > 1 else 0)
        het_dp.append(int(dp) if dp != "." else sum(ads))
if p.returncode:
    sys.exit(f"bcftools failed on {vcf}")

n_sv = int(subprocess.run(f"{BCF} view -H -t {AUTOS} -f PASS {sv_vcf} | wc -l", shell=True,
                          capture_output=True, text=True, check=True).stdout.strip())

v = np.array(het_vaf)
r = np.array(het_ref)
a = np.array(het_alt)
d = np.array(het_dp)
r2 = rng.binomial(r, f_thin)
a2 = rng.binomial(a, f_thin)
tot2 = r2 + a2
ok = tot2 > 0
hom_ds = ((a2[ok] / tot2[ok]) > 0.9).sum() / ok.sum()

with open(f"{out}.metrics.tsv", "w") as fo:
    fo.write("sample_id\tmetric\tvalue\tunit\n")
    for m, val, u in (("snv_pass_autosome", n_snv, "count"), ("indel_pass_autosome", n_indel, "count"),
                      ("sv_pass_autosome", n_sv, "count"), ("het_snv_autosome", len(v), "count"),
                      ("het_hom_peak_frac", (v > 0.9).sum() / len(v), "fraction"),
                      ("het_hom_peak_frac_ds", hom_ds, "fraction"),
                      ("het_hom_peak_frac_ds_n_sites", int(ok.sum()), "count"),
                      ("het_hom_peak_frac_ds_thinning_factor", f_thin, "fraction")):
        val = int(val) if u == "count" else float(val)  # numpy 스칼라 repr(np.float64(...)) 방지
        fo.write(f"{sid}\t{m}\t{val!r}\t{u}\n")

edges = np.round(np.arange(0, 1.0001, 0.01), 2)
bins = np.clip(np.floor(v * 100).astype(int), 0, 99)
dbin = np.select([d < 20, d < 30, d < 40], ["lt20", "20-30", "30-40"], "40+")
with open(f"{out}.vaf_bins.tsv", "w") as fo:
    fo.write("sample_id\tdepth_bin\tbin_start\tbin_end\tn_sites\n")
    for db in ("all", "lt20", "20-30", "30-40", "40+"):
        sel = np.ones(len(v), bool) if db == "all" else dbin == db
        cnt = np.bincount(bins[sel], minlength=100)
        for i in range(100):
            if cnt[i]:
                fo.write(f"{sid}\t{db}\t{edges[i]}\t{edges[i + 1]}\t{cnt[i]}\n")
print(f"{sid}: snv={n_snv} indel={n_indel} sv={n_sv} het={len(v)} hom_peak={(v > 0.9).sum() / len(v):.6f} ds={hom_ds:.6f}")
