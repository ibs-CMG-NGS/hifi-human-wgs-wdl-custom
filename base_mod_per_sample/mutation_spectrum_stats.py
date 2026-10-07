"""Fig 2e 통계: 샘플 고유 SNV 중 C>A 비율의 조건 차이 (아티팩트 정제 집합 포함).

집합:
  all            : 모든 private SNV
  clean          : depth ≤ 샘플 평균 1.5배, 같은 샘플 후보와 10 kb 이상 떨어진 자리
  clean_VAF>=0.2 / clean_VAF<0.3 : clean의 VAF 구간
모델: quasi-binomial GLM  C>A 비율 ~ condition + run (가중치 = 샘플별 private SNV 수)
  대비: B−A(H2O2), C−B(급성에서 약물), D−C(약물 하 만성 vs 급성), D−A
  정확 permutation: run 안에서 조건 라벨을 바꾸는 (4!)^3 = 13,824가지 전부의 F 통계량
  (조건 이름만 바꾼 24가지는 같은 F를 주므로 가능한 최소 p = 24/13,824 = 0.0017)
사용: python mutation_spectrum_stats.py samples.tsv FIGDATA_DIR
"""
import itertools
import sys
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

warnings.filterwarnings("ignore")
meta = pd.read_csv(sys.argv[1], sep="\t")
fd = sys.argv[2]
grp = dict(zip(meta.sample_id, meta.group))
run = dict(zip(meta.sample_id, meta.run.astype(str)))
s = pd.read_csv(f"{fd}/mutspec.private_sites.flagged.tsv.gz", sep="\t")
clean = s[~s.high_dp & ~s.clustered]
SETS = {"all": s, "clean": clean, "clean_VAF>=0.2": clean[clean.vaf >= 0.2], "clean_VAF<0.3": clean[clean.vaf < 0.3]}
KEY = "C(group, Treatment('A'))[T.{}]"


def table(sub):
    t = sub.groupby("sample").agg(CA=("class6", lambda x: (x == "C>A").sum()), tot=("class6", "size"))
    t = t.reindex(meta.sample_id)
    t["frac"] = t.CA / t.tot
    t["group"] = [grp[i] for i in t.index]
    t["run"] = [run[i] for i in t.index]
    return t


def fit(d):
    return smf.glm("frac ~ C(group, Treatment('A')) + run", d, family=sm.families.Binomial(),
                   var_weights=d.tot).fit(scale="X2")


def fstat(d):
    m = fit(d)
    m0 = smf.glm("frac ~ run", d, family=sm.families.Binomial(), var_weights=d.tot).fit()
    return (m0.deviance - m.deviance) / 3 / m.scale, m


runs = [[i for i in meta.sample_id if run[i] == r] for r in ("1", "2", "3")]
perms = [list(itertools.permutations(r)) for r in runs]
rows, lines = [], []
for label, sub in SETS.items():
    d = table(sub)
    F, m = fstat(d)
    p, cov = m.params, m.cov_params()
    con = {}
    for nm, v in (("B-A", {"B": 1}), ("C-B", {"C": 1, "B": -1}), ("D-C", {"D": 1, "C": -1}), ("D-A", {"D": 1})):
        vec = pd.Series(0.0, index=p.index)
        for k, w in v.items():
            vec[KEY.format(k)] = w
        e = vec @ p
        se = np.sqrt(vec @ cov @ vec)
        con[nm] = (np.exp(e), 2 * stats.t.sf(abs(e / se), m.df_resid))
    ge = 0
    for a in perms[0]:
        for b in perms[1]:
            for c in perms[2]:
                mp = {**dict(zip(runs[0], [grp[x] for x in a])), **dict(zip(runs[1], [grp[x] for x in b])),
                      **dict(zip(runs[2], [grp[x] for x in c]))}
                dd = d.copy()
                dd["group"] = [mp[i] for i in dd.index]
                ge += fstat(dd)[0] >= F - 1e-9
    n_perm = len(perms[0]) * len(perms[1]) * len(perms[2])
    frac = {g: d[d.group == g].CA.sum() / d[d.group == g].tot.sum() for g in "ABCD"}
    row = dict(set=label, n_sites=int(d.tot.sum()), **{f"CA_frac_{g}": frac[g] for g in "ABCD"},
               overdispersion=m.scale, F=F, p_F=stats.f.sf(F, 3, m.df_resid), p_perm_exact=ge / n_perm,
               **{f"OR_{k}": v[0] for k, v in con.items()}, **{f"p_{k}": v[1] for k, v in con.items()})
    rows.append(row)
    lines.append(f"## {label} (n={row['n_sites']:,}): C>A " + " ".join(f"{g}={frac[g]:.3f}" for g in "ABCD")
                 + f" | F={F:.2f} p_F={row['p_F']:.4f} exact perm p={row['p_perm_exact']:.4f} | "
                 + " ".join(f"{k} OR={v[0]:.2f} p={v[1]:.3f}" for k, v in con.items()))
    print(lines[-1], flush=True)
    d.to_csv(f"{fd}/mutspec.CA_per_sample.{label.replace('>=', 'ge').replace('<', 'lt')}.tsv", sep="\t")
pd.DataFrame(rows).to_csv(f"{fd}/mutspec.CA_stats.tsv", sep="\t", index=False)
open(f"{fd}/mutspec.CA_stats.txt", "w").write("\n".join(lines) + "\n")
