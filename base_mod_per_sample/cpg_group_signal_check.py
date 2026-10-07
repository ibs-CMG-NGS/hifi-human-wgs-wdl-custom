"""12샘플 CpG 5mC에서 그룹 신호가 우연 수준보다 큰지 빠르게 확인.

- 입력: pb-CpG-tools cpg_combined_bed (WGS), 상염색체만, 12샘플 모두 cov>=MIN_COV인 CpG
- 1) 샘플별 전역 평균 메틸화
- 2) PCA (가장 변이가 큰 CpG 상위 N개)
- 3) 사이트별 2-way 가법 ANOVA(group + run, 셀당 1개 관측 → 잔차 df=6)에서
     group F-test p<α 사이트 수를, run 안에서 그룹 라벨을 섞은 permutation과 비교
"""
import gzip
import sys

import numpy as np
import pandas as pd
from scipy import stats

SAMPLES = sys.argv[1]
OUT = sys.argv[2]
MIN_COV = 10
TOP_VAR = 200_000
N_TEST = 2_000_000
N_PERM = 100
ALPHAS = (1e-3, 1e-4)
AUTOSOMES = [f"chr{i}" for i in range(1, 20)]
CHR_IDX = {c: i for i, c in enumerate(AUTOSOMES)}

meta = pd.read_csv(SAMPLES, sep="\t")


def cpg_path(wgs_bam):
    run_dir = wgs_bam.split("/out/")[0]
    sid = wgs_bam.split("/")[-1].replace(".GRCm39.haplotagged.bam", "")
    return f"{run_dir}/out/cpg_combined_bed/{sid}.GRCm39.cpg_pileup.combined.bed.gz"


keys, meth = [], []
for _, r in meta.iterrows():
    df = pd.read_csv(cpg_path(r.wgs_haplotagged_bam), sep="\t", comment="#",
                     usecols=[0, 1, 3, 5], names=["chrom", "begin", "mod", "cov"], header=None,
                     dtype={"chrom": str}, compression="gzip", low_memory=False)
    df = df[df["chrom"].isin(CHR_IDX) & (pd.to_numeric(df["cov"], errors="coerce") >= MIN_COV)]
    k = df["chrom"].map(CHR_IDX).to_numpy(np.int64) * 1_000_000_000 + df["begin"].astype(np.int64).to_numpy()
    keys.append(k)
    meth.append(pd.Series(df["mod"].astype(float).to_numpy() / 100.0, index=k))
    print(f"loaded {r.sample_id}: {len(k):,} autosomal CpGs cov>={MIN_COV}", flush=True)

common = keys[0]
for k in keys[1:]:
    common = np.intersect1d(common, k, assume_unique=True)
M = np.column_stack([m.loc[common].to_numpy() for m in meth])  # sites × 12
print(f"common CpGs (all 12 cov>={MIN_COV}): {len(common):,}", flush=True)

sid = meta["sample_id"].tolist()
group = meta["group"].to_numpy()
run = meta["run"].to_numpy()

lines = []
lines.append("## 1. 샘플별 전역 평균 메틸화 (공통 CpG)")
for s, g, rn, v in zip(sid, group, run, M.mean(axis=0)):
    lines.append(f"{s}\tgroup={g}\trun={rn}\tmean5mC={v:.4f}")

# PCA
var = M.var(axis=1)
top = np.argsort(var)[-TOP_VAR:]
X = M[top].T
X = X - X.mean(axis=0)
U, S, Vt = np.linalg.svd(X, full_matrices=False)
pc = U * S
ev = S**2 / np.sum(S**2)
lines.append(f"\n## 2. PCA (변이 상위 {TOP_VAR:,} CpG) — 설명분산 PC1={ev[0]:.1%} PC2={ev[1]:.1%} PC3={ev[2]:.1%}")
for i, s in enumerate(sid):
    lines.append(f"{s}\tgroup={group[i]}\trun={run[i]}\tPC1={pc[i,0]:+.2f}\tPC2={pc[i,1]:+.2f}\tPC3={pc[i,2]:+.2f}")


def group_pvals(Y, grp, rn):
    """가법 2-way ANOVA group F-test p (balanced 4×3, 셀당 1개)."""
    gm = Y.mean(axis=1, keepdims=True)
    G = np.stack([Y[:, grp == g].mean(axis=1) for g in "ABCD"], axis=1)
    R = np.stack([Y[:, rn == x].mean(axis=1) for x in (1, 2, 3)], axis=1)
    gi = np.array(["ABCD".index(g) for g in grp])
    ri = np.array([(1, 2, 3).index(x) for x in rn])
    fitted = G[:, gi] + R[:, ri] - gm
    ss_res = ((Y - fitted) ** 2).sum(axis=1)
    ss_grp = 3 * ((G - gm) ** 2).sum(axis=1)
    F = (ss_grp / 3) / (ss_res / 6 + 1e-12)
    return stats.f.sf(F, 3, 6)


rng = np.random.default_rng(1)
sub = rng.choice(M.shape[0], size=min(N_TEST, M.shape[0]), replace=False)
Y = M[sub]
p_obs = group_pvals(Y, group, run)
obs = {a: int((p_obs < a).sum()) for a in ALPHAS}

perm_counts = {a: [] for a in ALPHAS}
for _ in range(N_PERM):
    g2 = group.copy()
    for x in (1, 2, 3):
        idx = np.where(run == x)[0]
        g2[idx] = rng.permutation(group[idx])
    p = group_pvals(Y, g2, run)
    for a in ALPHAS:
        perm_counts[a].append(int((p < a).sum()))

lines.append(f"\n## 3. group 효과 사이트 수: 실제 라벨 vs run 내 permutation ({N_PERM}회), 무작위 {len(sub):,} CpG")
for a in ALPHAS:
    pc_ = np.array(perm_counts[a])
    emp = (1 + (pc_ >= obs[a]).sum()) / (1 + N_PERM)
    lines.append(f"p<{a:g}\tobserved={obs[a]:,}\tperm_median={np.median(pc_):,.0f}\tperm_95%={np.percentile(pc_,95):,.0f}\t"
                 f"perm_max={pc_.max():,}\tratio_obs/median={obs[a]/max(np.median(pc_),1):.1f}\tempirical_p={emp:.3f}")

open(OUT, "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
