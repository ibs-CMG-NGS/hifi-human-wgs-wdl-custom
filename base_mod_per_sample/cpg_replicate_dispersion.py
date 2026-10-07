"""Fig 4d: 조건별 반복 샘플 간 5mC 거리 (클론 균질화와 메틸화 수렴의 연결).

- 입력: pb-CpG-tools cpg_combined_bed, 상염색체, 12샘플 모두 cov>=MIN_COV인 CpG
- 샘플 쌍마다 depth 편향을 뺀 제곱거리:
    D(i,j) = mean[(m_i - m_j)^2 - m_i(1-m_i)/(n_i-1) - m_j(1-m_j)/(n_j-1)]
  (m(1-m)/(n-1)은 이항 표본오차 분산의 불편추정량 → 저 depth 샘플의 거리 부풀림 제거)
- 같은 조건 쌍(서로 다른 run, 조건당 3쌍)과 다른 조건 쌍 중 서로 다른 run인 쌍을 비교
- 민감도: cov 하한 10 / 20 두 가지
출력: 쌍별 거리 tsv, 조건별 요약
사용: python cpg_replicate_dispersion.py samples.tsv OUT_PREFIX
"""
import itertools
import sys

import numpy as np
import pandas as pd

SAMPLES, OUT = sys.argv[1], sys.argv[2]
AUTOSOMES = [f"chr{i}" for i in range(1, 20)]
CHR_IDX = {c: i for i, c in enumerate(AUTOSOMES)}
COV_LEVELS = (10, 20)

meta = pd.read_csv(SAMPLES, sep="\t")


def cpg_path(wgs_bam):
    run_dir = wgs_bam.split("/out/")[0]
    sid = wgs_bam.split("/")[-1].replace(".GRCm39.haplotagged.bam", "")
    return f"{run_dir}/out/cpg_combined_bed/{sid}.GRCm39.cpg_pileup.combined.bed.gz"


series = []
for _, r in meta.iterrows():
    df = pd.read_csv(cpg_path(r.wgs_haplotagged_bam), sep="\t", comment="#",
                     usecols=[0, 1, 3, 5], names=["chrom", "begin", "mod", "cov"], header=None,
                     dtype={"chrom": str}, compression="gzip", low_memory=False)
    df = df[df["chrom"].isin(CHR_IDX)]
    cov = pd.to_numeric(df["cov"], errors="coerce").fillna(0).to_numpy(np.int32)
    keep = cov >= min(COV_LEVELS)
    k = df["chrom"].map(CHR_IDX).to_numpy(np.int64)[keep] * 1_000_000_000 + df["begin"].to_numpy(np.int64)[keep]
    m = df["mod"].astype(np.float32).to_numpy()[keep] / 100.0
    o = np.argsort(k, kind="stable")  # bed는 chr1, chr10, … 사전순이라 키가 단조 증가하지 않음
    series.append((k[o], m[o], cov[keep][o]))
    print(f"loaded {r.sample_id}: {keep.sum():,}", flush=True)

common = series[0][0]
for k, _, _ in series[1:]:
    common = np.intersect1d(common, k, assume_unique=True)
M = np.empty((len(common), len(series)), np.float32)
N = np.empty((len(common), len(series)), np.int32)
for j, (k, m, c) in enumerate(series):
    idx = np.searchsorted(k, common)
    M[:, j], N[:, j] = m[idx], c[idx]
print(f"common CpGs: {len(common):,}", flush=True)

sid = meta["sample_id"].tolist()
grp = meta["group"].tolist()
run = meta["run"].tolist()
rows = []
for cmin in COV_LEVELS:
    sel = (N >= cmin).all(axis=1)
    Ms, Ns = M[sel].astype(np.float64), N[sel].astype(np.float64)
    noise = Ms * (1 - Ms) / (Ns - 1)
    for i, j in itertools.combinations(range(len(sid)), 2):
        raw = (Ms[:, i] - Ms[:, j]) ** 2
        rows.append({"min_cov": cmin, "n_cpg": int(sel.sum()), "s1": sid[i], "s2": sid[j],
                     "g1": grp[i], "g2": grp[j], "run1": run[i], "run2": run[j],
                     "raw_sqdist": raw.mean(), "corrected_sqdist": (raw - noise[:, i] - noise[:, j]).mean(),
                     "mean_abs_diff": np.abs(Ms[:, i] - Ms[:, j]).mean()})
d = pd.DataFrame(rows)
d.to_csv(f"{OUT}.pairs.tsv", sep="\t", index=False)

lines = []
for cmin, g in d.groupby("min_cov"):
    lines.append(f"## min_cov={cmin}  (n_CpG={g['n_cpg'].iat[0]:,})  — 값은 corrected_sqdist ×1e3 (RMS 거리 환산은 sqrt)")
    within = g[(g.g1 == g.g2)]
    between = g[(g.g1 != g.g2) & (g.run1 != g.run2)]
    for c in "ABCD":
        w = within[within.g1 == c]
        vals = ", ".join(f"{v * 1e3:.3f}" for v in w["corrected_sqdist"])
        lines.append(f"  within {c}: mean={w['corrected_sqdist'].mean() * 1e3:.3f}  RMS={np.sqrt(max(w['corrected_sqdist'].mean(), 0)):.4f}  pairs=[{vals}]")
    lines.append(f"  between conditions (different runs): mean={between['corrected_sqdist'].mean() * 1e3:.3f}  "
                 f"range={between['corrected_sqdist'].min() * 1e3:.3f}-{between['corrected_sqdist'].max() * 1e3:.3f}")
    for c in "ABCD":
        b = g[((g.g1 == c) | (g.g2 == c)) & (g.g1 != g.g2) & (g.run1 != g.run2)]
        lines.append(f"  between {c}-vs-others: mean={b['corrected_sqdist'].mean() * 1e3:.3f}")
    # 샘플별 평균 거리 (이상치 샘플 확인)
    per = {s: np.mean(g[(g.s1 == s) | (g.s2 == s)]["corrected_sqdist"]) for s in sid}
    lines.append("  per-sample mean distance to all others: " +
                 ", ".join(f"{s}={v * 1e3:.3f}" for s, v in per.items()))
open(f"{OUT}.summary.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
