"""Fig 6a: 샘플별 전역 5hmC 수준과 조건 대비 (exact within-run permutation), 라이브러리 교란 지표와의 상관.

입력
  samples.tsv                     : sample_id, run, group
  COUNT_DIR                       : dmr_region_counts.py ... 5hmC 출력 ({set}.meth/cov.tsv.gz)
  TECH                            : str_contraction_chr1.tsv (sample, contraction_idx_2bp)
지표 (샘플별, 상염색체)
  global      : 전체 CpG의 5hmC read 합 / coverage 합 (HMC_DIR의 CpG TSV에서 직접)
  {set}       : 영역 집합 안 5hmC read 합 / coverage 합 (promoter, genebody, cgi, shore, tile1k)
검정: 값을 run 안에서 중심화한 뒤 대비(H2O2=B−A, drug=C−B, chronic=D−C, DvsA=D−A)의 군 평균 차이,
      run 안 조건 라벨 13,824 배치 전수로 양측 p. 교란: run 중심화 값과 logit(STR 수축 지수)의 Pearson r.
출력: OUT_PREFIX.samples.tsv, OUT_PREFIX.tests.tsv
사용: python hmc_global.py samples.tsv COUNT_DIR TECH OUT_PREFIX
"""
import itertools
import sys

import numpy as np
import pandas as pd

HMC_DIR = "/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/5hmC"
AUTO = {f"chr{i}" for i in range(1, 20)}
SETS = ["promoter", "genebody", "cgi", "shore", "tile1k"]
CON = {"H2O2": ("B", "A"), "drug": ("C", "B"), "chronic": ("D", "C"), "DvsA": ("D", "A")}

meta = pd.read_csv(sys.argv[1], sep="\t")[["sample_id", "run", "group"]]
cdir, tech_path, out = sys.argv[2], sys.argv[3], sys.argv[4]

rows = []
for sid in meta.sample_id:
    m = c = 0.0
    for ch in pd.read_csv(f"{HMC_DIR}/{sid}/_LAST/call-modbam_cpg_pileup/work/{sid}.5hmC.tsv.gz", sep="\t",
                          header=None, names=["chrom", "pos", "cov", "mod"], usecols=[0, 2, 3],
                          dtype={"chrom": str}, chunksize=5_000_000):
        ch = ch[ch.chrom.isin(AUTO)]
        m += ch["mod"].sum(); c += ch["cov"].sum()
    rows.append(dict(sample_id=sid, global_=m / c))
    print(f"loaded {sid}", flush=True)
df = meta.merge(pd.DataFrame(rows), on="sample_id").rename(columns={"global_": "global"})
for s in SETS:
    M = pd.read_csv(f"{cdir}/{s}.meth.tsv.gz", sep="\t", index_col=0)
    C = pd.read_csv(f"{cdir}/{s}.cov.tsv.gz", sep="\t", index_col=0)
    df[s] = (M.sum() / C.sum()).reindex(df.sample_id).to_numpy()
tech = pd.read_csv(tech_path, sep="\t").set_index("sample")["contraction_idx_2bp"]
df["tech_logit"] = np.log(tech / (1 - tech)).reindex(df.sample_id).to_numpy()
df.to_csv(f"{out}.samples.tsv", sep="\t", index=False)

# run 안 라벨 배치 전수
runs = [df.index[df.run == r].to_numpy() for r in sorted(df.run.unique())]
lab = df.group.to_numpy()
perms = []
for combo in itertools.product(*[itertools.permutations(ix) for ix in runs]):
    g = lab.copy()
    for ix, p in zip(runs, combo):
        g[ix] = lab[list(p)]
    perms.append(g)
perms = np.array(perms)
print(f"permutations: {len(perms)}")

res = []
tc = df.tech_logit - df.groupby("run").tech_logit.transform("mean")
for col in ["global"] + SETS:
    v = (df[col] - df.groupby("run")[col].transform("mean")).to_numpy()
    for k, (x, y) in CON.items():
        stat = lambda g: v[g == x].mean() - v[g == y].mean()
        obs = stat(lab)
        null = np.array([stat(g) for g in perms])
        res.append(dict(metric=col, contrast=k, diff=obs, rel_diff=obs / df[col].mean(),
                        perm_p=(np.abs(null) >= abs(obs) - 1e-15).mean()))
    res.append(dict(metric=col, contrast="r_tech", diff=np.corrcoef(v, tc)[0, 1], rel_diff=np.nan, perm_p=np.nan))
res = pd.DataFrame(res)
res.to_csv(f"{out}.tests.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(df.round(5).to_string(index=False))
print(df.groupby("group")[["global"] + SETS].mean().round(5).to_string())
print(res.round(5).to_string(index=False))
