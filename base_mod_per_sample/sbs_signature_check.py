"""Fig 2e 보강: private SNV 96-context가 산화 손상 signature(SBS18)와 닮았는지 확인.

- 그룹별 private SNV 96-context를 이형접합 SNV 10만 개당으로 정규화
- 처치 조건의 초과 프로파일 = 조건 − 무처치(A), 음수는 0으로
- C>A 16채널과 96채널 전체에서 COSMIC v3.4(mm39 정규화) 서명과 코사인 유사도
- 기준선: 생식세포(12샘플 공통) 프로파일, 무처치 private 프로파일 자체
- 무작위 기준: 초과 프로파일의 C>A 16채널을 섞었을 때(1만 회) SBS18 코사인 분포
사용: python sbs_signature_check.py FIGDATA_DIR COSMIC_TXT samples_metrics.tsv
"""
import sys

import numpy as np
import pandas as pd

figdir, cosmic_txt, metrics = sys.argv[1:4]
cos = pd.read_csv(cosmic_txt, sep="\t").set_index("Type")
ctx = pd.read_csv(f"{figdir}/mutspec.ctx96.tsv", sep="\t")
met = pd.read_csv(metrics, sep="\t").set_index("sample")
ctx["group"] = ctx["sample"].map(met["group"])
CA = [c for c in cos.index if "[C>A]" in c]
FOCUS = ["SBS18", "SBS36", "SBS45", "SBS4", "SBS24", "SBS29", "SBS1", "SBS5", "SBS40a"]
rng = np.random.default_rng(1)


def cosine(a, b):
    return float(a @ b / np.sqrt((a @ a) * (b @ b)))


def report(name, v):
    v = v.reindex(cos.index).fillna(0).clip(lower=0)
    sims96 = {s: cosine(v.to_numpy(), cos[s].to_numpy()) for s in cos.columns}
    va = v[CA].to_numpy()
    sims16 = {s: cosine(va, cos.loc[CA, s].to_numpy()) for s in cos.columns}
    top96 = sorted(sims96.items(), key=lambda x: -x[1])[:5]
    top16 = sorted(sims16.items(), key=lambda x: -x[1])[:5]
    print(f"\n## {name}  (C>A share of profile = {va.sum() / max(v.sum(), 1e-9):.2f})")
    print("  96ch top: " + ", ".join(f"{k}={x:.2f}" for k, x in top96))
    print("  C>A16 top: " + ", ".join(f"{k}={x:.2f}" for k, x in top16))
    print("  focus (96ch / C>A16): " + ", ".join(f"{s}={sims96[s]:.2f}/{sims16[s]:.2f}" for s in FOCUS if s in cos.columns))
    return va


het = met.groupby("group")["n_het_snv"].sum()
priv = ctx[ctx["set"] == "private"].groupby(["group", "ctx96"])["n"].sum().unstack(fill_value=0)
norm = priv.div(het, axis=0) * 1e5
germ = ctx[ctx["set"] == "shared_all"].groupby("ctx96")["n"].sum() / 12

report("germline (shared by all 12)", germ)
report("private, untreated A", norm.loc["A"])
for g, lab in (("B", "acute H2O2"), ("C", "acute H2O2 + drug"), ("D", "chronic H2O2 + drug")):
    report(f"private, {lab} {g}", norm.loc[g])
for g, lab in (("B", "acute H2O2"), ("C", "acute H2O2 + drug"), ("D", "chronic H2O2 + drug"), ("BD", "B+D pooled")):
    if g == "BD":
        exc = (priv.loc[["B", "D"]].sum() / het[["B", "D"]].sum() * 1e5) - norm.loc["A"]
    else:
        exc = norm.loc[g] - norm.loc["A"]
    exc = exc.clip(lower=0)
    va = report(f"EXCESS over A: {lab}", exc)
    ref = cos.loc[CA, "SBS18"].to_numpy()
    obs = cosine(va, ref)
    null = np.array([cosine(rng.permutation(va), ref) for _ in range(10000)])
    print(f"  SBS18 C>A16 cosine={obs:.2f}; shuffled-channel null median={np.median(null):.2f}, "
          f"95%={np.quantile(null, .95):.2f}, p={(1 + (null >= obs).sum()) / 10001:.4f}; excess C>A per 100k het={va.sum():.1f}")
# C>A 16채널 그룹별 표 (그림용)
tab = norm[CA].T
tab["SBS18"] = cos.loc[CA, "SBS18"] / cos.loc[CA, "SBS18"].sum()
tab.to_csv(f"{figdir}/mutspec.CA16_per100k_het.tsv", sep="\t")
print(f"\nwrote {figdir}/mutspec.CA16_per100k_het.tsv")
print(tab.round(3).to_string())
