"""Fig 3h 내보내기: 조건별 '3중 겹침' hom-peak 사이트와 조건 간 공유 패턴.

정의 (H2O2_LOH_hotspot_cross_validation_2026-07-22.md 방법 1–5와 동일):
  1. 샘플마다 상염색체 PASS 이형접합 SNV 중 FORMAT/VAF > 0.9인 사이트
  2. 조건마다 3개 반복 모두에 있는 사이트(3중 겹침)
  3. 사이트를 조건 간에 비교해 pattern = 그 사이트가 3중 겹침인 조건 조합(A_only, B_D, A_B_C_D …)
  4. 위치 ±5 kb 안의 GENCODE vM36 유전자(gene_name; 없으면 gene_id)
n_supporting_samples = 12샘플 중 그 사이트가 het VAF > 0.9인 샘플 수.
같은 패턴의 사이트를 10 kb 이내로 묶은 locus_id도 함께 준다.
출력: OUT.tsv
사용: python loh_hotspots.py samples.tsv OUT
"""
import subprocess
import sys

import pandas as pd

BCF = "/home/ygkim/program/anaconda3/envs/wgs-tertiary-pipeline/bin/bcftools"
GTF = "/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/gencode.vM36.annotation.gtf.gz"
AUTOS = ",".join(f"chr{i}" for i in range(1, 20))
meta = pd.read_csv(sys.argv[1], sep="\t")
out = sys.argv[2]

sites = {}
for _, r in meta.iterrows():
    o = r.wgs_haplotagged_bam.split("/out/")[0] + "/out"
    vcf = f"{o}/phased_small_variant_vcf/{r.sample_id}.GRCm39.small_variants.phased.vcf.gz"
    res = subprocess.run([BCF, "query", "-t", AUTOS, "-i", 'FILTER="PASS"', "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%VAF]\n", vcf],
                         capture_output=True, text=True, check=True).stdout
    s = set()
    for line in res.splitlines():
        c, p, ref, alt, gt, vaf = line.split("\t")
        if len(ref) != 1 or any(len(a) != 1 for a in alt.split(",")):
            continue
        al = gt.replace("|", "/").split("/")
        if len(al) != 2 or al[0] == al[1] or "." in al:
            continue
        try:
            v = max(float(x) for x in vaf.split(",") if x != ".")
        except ValueError:
            continue
        if v > 0.9:
            s.add((c, int(p)))
    sites[r.sample_id] = s
    print(f"{r.sample_id}: {len(s)} het VAF>0.9 SNVs", flush=True)

tri = {}
for g, sub in meta.groupby("group"):
    ss = [sites[x] for x in sub.sample_id]
    tri[g] = ss[0] & ss[1] & ss[2]
    print(f"condition {g}: {len(tri[g])} triple-overlap sites", flush=True)
allsites = sorted(set().union(*tri.values()))

gtf = pd.read_csv(GTF, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 8],
                  names=["chrom", "type", "start", "end", "attr"])
gtf = gtf[gtf.type == "gene"]
gtf["name"] = gtf.attr.str.extract(r'gene_name "([^"]+)"', expand=False).fillna(
    gtf.attr.str.extract(r'gene_id "([^"]+)"', expand=False))
rows = []
for c, p in allsites:
    pat = "_".join(g for g in "ABCD" if (c, p) in tri[g])
    gg = gtf[(gtf.chrom == c) & (gtf.start <= p + 5000) & (gtf.end >= p - 5000)]
    rows.append(dict(chrom=c, pos=p, pattern=pat, n_conditions=len(pat.split("_")),
                     n_supporting_samples=sum((c, p) in s for s in sites.values()),
                     genes=";".join(dict.fromkeys(gg.name))))
d = pd.DataFrame(rows).sort_values(["pattern", "chrom", "pos"])
lid, prev = 0, None
ids = []
for _, r in d.iterrows():
    if prev is None or r.pattern != prev[0] or r.chrom != prev[1] or r.pos - prev[2] > 10_000:
        lid += 1
    ids.append(f"L{lid:03d}")
    prev = (r.pattern, r.chrom, r.pos)
d["locus_id"] = ids
d.to_csv(out, sep="\t", index=False)
print(d.pattern.value_counts().to_string())
