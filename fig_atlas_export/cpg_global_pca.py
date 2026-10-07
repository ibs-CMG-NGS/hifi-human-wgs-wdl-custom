"""Fig 4a·4b 내보내기: 12샘플 공통 CpG(상염색체, 모두 cov ≥ 10)의 샘플별 평균 5mC와 PCA 점수.

base_mod_per_sample/cpg_group_signal_check.py와 같은 CpG 집합·PCA 정의(변이 상위 200,000 CpG,
샘플 간 평균 중심화, SVD)이며, 반올림 없이 PC1–5와 설명분산을 저장한다.
출력: OUT_DIR/global_5mC.tsv, OUT_DIR/pca_scores.tsv, OUT_DIR/pca_variance.tsv
사용: python cpg_global_pca.py samples.tsv OUT_DIR
"""
import sys

import numpy as np
import pandas as pd

meta = pd.read_csv(sys.argv[1], sep="\t")
out = sys.argv[2]
MIN_COV, TOP_VAR, NPC = 10, 200_000, 5
AUTOSOMES = [f"chr{i}" for i in range(1, 20)]
CHR_IDX = {c: i for i, c in enumerate(AUTOSOMES)}


def cpg_path(wgs_bam):
    run_dir = wgs_bam.split("/out/")[0]
    sid = wgs_bam.split("/")[-1].replace(".GRCm39.haplotagged.bam", "")
    return f"{run_dir}/out/cpg_combined_bed/{sid}.GRCm39.cpg_pileup.combined.bed.gz"


series = []
for _, r in meta.iterrows():
    df = pd.read_csv(cpg_path(r.wgs_haplotagged_bam), sep="\t", comment="#", usecols=[0, 1, 3, 5],
                     names=["chrom", "begin", "mod", "cov"], header=None, dtype={"chrom": str},
                     compression="gzip", low_memory=False)
    df = df[df.chrom.isin(CHR_IDX) & (pd.to_numeric(df["cov"], errors="coerce") >= MIN_COV)]
    k = df.chrom.map(CHR_IDX).to_numpy(np.int64) * 1_000_000_000 + df.begin.to_numpy(np.int64)
    o = np.argsort(k, kind="stable")
    series.append((k[o], df["mod"].to_numpy(np.float64)[o] / 100.0))
    print(f"loaded {r.sample_id}", flush=True)
common = series[0][0]
for k, _ in series[1:]:
    common = np.intersect1d(common, k, assume_unique=True)
M = np.column_stack([m[np.searchsorted(k, common)] for k, m in series])
print(f"common CpGs: {len(common)}", flush=True)

sid = meta.sample_id.tolist()
pd.DataFrame({"sample_id": sid, "global_5mC_mean": M.mean(axis=0), "global_5mC_n_cpg": len(common)}).to_csv(
    f"{out}/global_5mC.tsv", sep="\t", index=False, float_format="%.17g")
top = np.argsort(M.var(axis=1))[-TOP_VAR:]
X = M[top].T
X = X - X.mean(axis=0)
U, S, _ = np.linalg.svd(X, full_matrices=False)
pc = U * S
ev = S ** 2 / np.sum(S ** 2)
scores = pd.DataFrame(pc[:, :NPC], columns=[f"PC{i + 1}" for i in range(NPC)])
scores.insert(0, "sample_id", sid)
scores.to_csv(f"{out}/pca_scores.tsv", sep="\t", index=False, float_format="%.17g")
pd.DataFrame({"pc": [f"PC{i + 1}" for i in range(len(ev))], "variance_explained": ev}).to_csv(
    f"{out}/pca_variance.tsv", sep="\t", index=False, float_format="%.17g")
print("done")
