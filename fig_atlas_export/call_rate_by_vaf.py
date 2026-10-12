"""같은 read 수준 VAF에서 DeepVariant가 이형접합으로 호출하는 비율(호출 민감도)을 샘플별로 비교.
vaf_shift_pairs.py의 사이트 × 샘플 alt/depth 표를 쓰고, 품질 맞춘 VCF의 PASS het/hom 호출을 붙인다.
사용: python call_rate_by_vaf.py QM_DIR SITES_TSV_GZ OUT_TSV
"""
import subprocess, sys
import numpy as np, pandas as pd
BCF = "/home/ygkim/program/anaconda3/envs/wgs-tertiary-pipeline/bin/bcftools"
QM, SITES, OUT = sys.argv[1:4]
meta = pd.read_csv("/data_4tb/hifi-human-wgs-wdl-custom/base_mod_per_sample/samples.tsv", sep="\t")
t = pd.read_csv(SITES, sep="\t", index_col=[0, 1])
rows = []
for r in meta.itertuples():
    s = r.sample_id
    res = subprocess.run([BCF, "query", "-t", "chr1,chr2,chr3", "-i", 'FILTER="PASS" && TYPE="snp"', "-f", "%CHROM\t%POS[\t%GT]\n",
                          f"{QM}/dv/{s}/{s}.vcf.gz"], capture_output=True, text=True, check=True).stdout
    called = {}
    for line in res.splitlines():
        c, p, gt = line.split("\t")
        al = gt.replace("|", "/").split("/")
        called[(c, int(p))] = "het" if al[0] != al[1] else "hom"
    a, d = t[f"{s}.alt"].to_numpy(float), t[f"{s}.dp"].to_numpy(float)
    v = np.where(d > 0, a / np.maximum(d, 1), np.nan)
    st = np.array([called.get(k, "none") for k in t.index])
    for lo, hi in ((0.05, 0.15), (0.15, 0.25), (0.25, 0.35), (0.35, 0.45), (0.45, 0.6)):
        m = (d >= 20) & (d <= 60) & (v >= lo) & (v < hi)
        rows.append(dict(sample_id=s, condition=r.group, run=r.run, vaf_bin=f"{lo}-{hi}", n=int(m.sum()),
                         het_call_rate=float((st[m] == "het").mean()) if m.sum() else np.nan))
o = pd.DataFrame(rows)
o.to_csv(OUT, sep="\t", index=False)
print(o.pivot_table(index=["vaf_bin"], columns=["run", "condition"], values="het_call_rate").round(3).to_string())
print(o.pivot_table(index=["vaf_bin"], columns=["condition"], values="n", aggfunc="mean").round(0).to_string())
