"""Fig 6b: 5hmC 수준을 같은 샘플의 5mC 수준 구간별로 비교 — 전역 5hmC 차이가 생물학적인지 보정(calibration) 차이인지 판별.

근거: 5hmC는 5mC의 산화로만 생기므로, 5mC가 거의 없는 CpG(비메틸화 CGI 등)의 5hmC 호출은 대부분 위양성이다.
  조건 차이가 5mC 0–5 % 구간에서도 같은 비율로 나타나면 호출 보정 차이(기술적), 고메틸 구간에만 있으면 생물학적 신호에 가깝다.
입력: base_modification_per_sample.wdl 출력 5mC·5hmC CpG TSV (chr, pos, N_valid, N_mod), 같은 modkit pileup에서 나온 값
구간: 샘플 자신의 5mC 비율로 CpG를 나눔 (coverage ≥ 20)
출력: OUT.tsv (sample, run, group, mc_bin, n_cpg, hmc_frac)
사용: python hmc_by_mc_context.py samples.tsv OUT.tsv [CHROM ...]
"""
import sys

import numpy as np
import pandas as pd

HMC_DIR = "/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/5hmC"
BINS = [0, 0.05, 0.2, 0.5, 0.8, 0.9, 1.0001]
LABELS = ["0-5", "5-20", "20-50", "50-80", "80-90", "90-100"]
MIN_COV = 20

meta = pd.read_csv(sys.argv[1], sep="\t")
out = sys.argv[2]
chroms = set(sys.argv[3:]) if len(sys.argv) > 3 else {f"chr{i}" for i in range(1, 20)}


def load(sid, mod):
    df = pd.read_csv(f"{HMC_DIR}/{sid}/_LAST/call-modbam_cpg_pileup/work/{sid}.{mod}.tsv.gz", sep="\t", header=None,
                     names=["chrom", "pos", "cov", mod], dtype={"chrom": str})
    return df[df.chrom.isin(chroms)].set_index(["chrom", "pos"])


rows = []
for _, r in meta.iterrows():
    m, h = load(r.sample_id, "5mC"), load(r.sample_id, "5hmC")
    d = m.join(h["5hmC"], how="inner")
    d = d[d["cov"] >= MIN_COV]
    d["bin"] = pd.cut(d["5mC"] / d["cov"], BINS, labels=LABELS, right=False)
    g = d.groupby("bin", observed=True).agg(n_cpg=("cov", "size"), cov=("cov", "sum"), hmc=("5hmC", "sum"))
    for b, x in g.iterrows():
        rows.append(dict(sample_id=r.sample_id, run=r.run, group=r.group, mc_bin=b, n_cpg=x.n_cpg,
                         hmc_frac=x.hmc / x["cov"]))
    print(f"loaded {r.sample_id}", flush=True)
res = pd.DataFrame(rows)
res.to_csv(out, sep="\t", index=False)
piv = res.pivot_table(index=["run", "group"], columns="mc_bin", values="hmc_frac")[LABELS]
print(piv.round(5).to_string())
# 각 run에서 C/B 비 (같은 세포 풀 안 대비)
for b in LABELS:
    ratios = [piv.loc[(k, "C"), b] / piv.loc[(k, "B"), b] for k in sorted(res.run.unique())]
    print(f"C/B ratio {b:>6}: " + " ".join(f"{x:.3f}" for x in ratios))
print(res.pivot_table(index="sample_id", columns="mc_bin", values="n_cpg")[LABELS].to_string())
