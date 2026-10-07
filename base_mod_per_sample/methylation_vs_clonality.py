"""Fig 4f: 샘플 단위로 메틸화 이질성과 유전적 클론성 지표의 관계.

메틸화 지표(cpg_dispersion.pairs.tsv, min_cov=10, depth 보정 제곱거리):
  dist_all     : 다른 11개 샘플까지 평균 거리
  dist_within  : 같은 조건 반복 2개까지 평균 거리
  dist_between : 다른 run의 다른 조건 샘플까지 평균 거리 (run 효과 배제)
유전 지표(sample_metrics_combined.tsv): het_hom_peak_frac, n_het_snv
Spearman 상관 + run 효과를 뺀(run 평균 차감) 뒤의 상관
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

pairs = pd.read_csv(sys.argv[1], sep="\t")
met = pd.read_csv(sys.argv[2], sep="\t").set_index("sample")
out = sys.argv[3]
p = pairs[pairs.min_cov == 10]
both = pd.concat([p.rename(columns={"s1": "s", "s2": "o", "g1": "g", "g2": "go", "run1": "r", "run2": "ro"}),
                  p.rename(columns={"s2": "s", "s1": "o", "g2": "g", "g1": "go", "run2": "r", "run1": "ro"})])
d = pd.DataFrame({
    "dist_all": both.groupby("s").corrected_sqdist.mean(),
    "dist_within": both[both.g == both.go].groupby("s").corrected_sqdist.mean(),
    "dist_between": both[(both.g != both.go) & (both.r != both.ro)].groupby("s").corrected_sqdist.mean(),
}) * 1e3
d = d.join(met[["group", "het_hom_peak_frac", "n_het_snv"]])
d["run"] = [1 if not s.startswith("H2O2-2") and not s.startswith("H2O2-3") else int(s.split("-")[1]) for s in d.index]
d.to_csv(out, sep="\t")
print(d.round(4).to_string())
lines = []
for x in ("het_hom_peak_frac", "n_het_snv"):
    for y in ("dist_all", "dist_within", "dist_between"):
        r, pv = stats.spearmanr(d[x], d[y])
        rx = d[x] - d.groupby("run")[x].transform("mean")
        ry = d[y] - d.groupby("run")[y].transform("mean")
        r2, pv2 = stats.pearsonr(rx, ry)
        lines.append(f"{x} vs {y}: Spearman rho={r:+.2f} p={pv:.4f} | run-centered Pearson r={r2:+.2f} p={pv2:.4f}")
print("\n".join(lines))
open(out.replace(".tsv", ".stats.txt"), "w").write("\n".join(lines) + "\n")
