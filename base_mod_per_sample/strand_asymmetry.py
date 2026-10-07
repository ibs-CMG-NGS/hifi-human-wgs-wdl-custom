"""② 전사 가닥 비대칭: 샘플 고유 SNV(아티팩트 정제 집합)가 전사 가닥(T)과 비전사 가닥(U) 중
어느 쪽에 치우치는지, 조건별·치환 클래스별로 검정.

규칙 (SigProfiler 관례): 치환은 피리미딘(C/T) 기준으로 표기.
  유전자가 + 가닥이고 참조 염기가 피리미딘(C/T) → 피리미딘이 coding 가닥 = U(untranscribed)
  유전자가 + 가닥이고 참조 염기가 퓨린(G/A)     → 피리미딘이 template 가닥 = T(transcribed)
  (− 가닥 유전자는 반대). 양쪽 가닥 유전자가 겹치는 자리는 제외.
  C>A(U 우세) = 비전사 가닥에 G>T가 많음 = 8-oxoG가 비전사 가닥 G에서 남았다는 뜻
              (전사 연관 수선 TC-NER이 전사 가닥 손상을 우선 제거하면 기대되는 방향)
발현 근사: RNA-seq이 없어, 12샘플 평균 프로모터(TSS ± 1 kb) CpG 메틸화로 유전자를
  active(<0.2) / intermediate(0.2–0.6) / silent(≥0.6)로 구분.
내부 대조: C>A 외 클래스(대부분 생식세포 변이 → 가닥 편향 기대 없음)
사용: python strand_asymmetry.py samples.tsv FIGDATA_DIR
"""
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

meta = pd.read_csv(sys.argv[1], sep="\t")
fd = sys.argv[2]
GTF = "/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/gencode.vM36.annotation.gtf.gz"
AUTOSOMES = [f"chr{i}" for i in range(1, 20)]
PROM = 1000

g = pd.read_csv(GTF, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                names=["chrom", "type", "start", "end", "strand", "attr"])
g = g[(g.type == "gene") & g.chrom.isin(AUTOSOMES)]
g["gene_type"] = g.attr.str.extract(r'gene_type "([^"]+)"', expand=False)
g["gene_name"] = g.attr.str.extract(r'gene_name "([^"]+)"', expand=False)
g = g[g.gene_type.isin(["protein_coding", "lncRNA"])].reset_index(drop=True)
g["start"] -= 1
g["tss"] = np.where(g.strand == "+", g.start, g.end - 1)

# ---- 프로모터 메틸화 (발현 근사) ----
def cpg_path(wgs_bam):
    run_dir = wgs_bam.split("/out/")[0]
    sid = wgs_bam.split("/")[-1].replace(".GRCm39.haplotagged.bam", "")
    return f"{run_dir}/out/cpg_combined_bed/{sid}.GRCm39.cpg_pileup.combined.bed.gz"


prom_sum = np.zeros(len(g))
prom_n = np.zeros(len(g))
for _, r in meta.iterrows():
    c = pd.read_csv(cpg_path(r.wgs_haplotagged_bam), sep="\t", comment="#", usecols=[0, 1, 3, 5],
                    names=["chrom", "pos", "mod", "cov"], header=None, dtype={"chrom": str},
                    compression="gzip", low_memory=False)
    c = c[c.chrom.isin(AUTOSOMES) & (pd.to_numeric(c["cov"], errors="coerce") >= 5)]
    for ch, cc in c.groupby("chrom"):
        pos = np.sort(cc.pos.to_numpy())
        m = cc.set_index("pos")["mod"].astype(float).reindex(pos).to_numpy() / 100
        cs = np.concatenate([[0], np.cumsum(m)])
        idx = np.where(g.chrom == ch)[0]
        lo = np.searchsorted(pos, g.tss.to_numpy()[idx] - PROM)
        hi = np.searchsorted(pos, g.tss.to_numpy()[idx] + PROM + 1)
        prom_sum[idx] += cs[hi] - cs[lo]
        prom_n[idx] += hi - lo
    print(f"promoter CpGs loaded: {r.sample_id}", flush=True)
g["prom_meth"] = np.where(prom_n > 0, prom_sum / np.maximum(prom_n, 1), np.nan)
g["activity"] = pd.cut(g.prom_meth, [-0.01, 0.2, 0.6, 1.01], labels=["active", "intermediate", "silent"])
g[["chrom", "start", "end", "strand", "gene_name", "gene_type", "prom_meth", "activity"]].to_csv(
    f"{fd}/gene_promoter_methylation.tsv.gz", sep="\t", index=False)

# ---- 사이트 → 유전자 가닥 ----
s = pd.read_csv(f"{fd}/mutspec.private_sites.flagged.tsv.gz", sep="\t")
s = s[~s.high_dp & ~s.clustered].copy()
strand_of, act_of = [], []
for _, row in s.iterrows():
    gg = g[(g.chrom == row.chrom) & (g.start < row.pos) & (g.end >= row.pos)]
    st = set(gg.strand)
    if len(st) != 1:
        strand_of.append(None)
        act_of.append(None)
        continue
    strand_of.append(st.pop())
    # 겹치는 유전자가 여러 개면 가장 활성인 것
    a = gg.activity.dropna()
    act_of.append(None if a.empty else sorted(a.astype(str), key=["active", "intermediate", "silent"].index)[0])
s["gene_strand"] = strand_of
s["activity"] = act_of
s = s[s.gene_strand.notna()].copy()
pyr_ref = s.ref.isin(["C", "T"])
s["ts"] = np.where((s.gene_strand == "+") == pyr_ref, "U", "T")
s["is_CA"] = s.class6 == "C>A"
s.to_csv(f"{fd}/strand_asymmetry.sites.tsv.gz", sep="\t", index=False)

lines = []
def binom_line(label, d):
    u, t = int((d.ts == "U").sum()), int((d.ts == "T").sum())
    p = stats.binomtest(u, u + t, 0.5).pvalue if u + t else np.nan
    return f"{label}: U={u} T={t} U/(U+T)={u / max(u + t, 1):.3f} binom p={p:.4f}"

lines.append("## 1. 클래스별 (전 조건 합산)")
for c6, d in s.groupby("class6"):
    lines.append("  " + binom_line(c6, d))
lines.append("\n## 2. C>A, 조건별")
for gr, d in s[s.is_CA].groupby("group"):
    lines.append("  " + binom_line(f"group {gr}", d))
lines.append("\n## 3. C>A vs 다른 클래스 (내부 대조), 조건별 Fisher")
for gr, d in s.groupby("group"):
    a = d[d.is_CA]
    b = d[~d.is_CA]
    tab = [[(a.ts == "U").sum(), (a.ts == "T").sum()], [(b.ts == "U").sum(), (b.ts == "T").sum()]]
    orr, p = stats.fisher_exact(tab)
    lines.append(f"  group {gr}: C>A U/T={tab[0]}, others U/T={tab[1]}, OR={orr:.2f} p={p:.4f}")
lines.append("\n## 4. C>A, 발현 근사(프로모터 메틸화) 구간별 (전 조건 합산 / 처치 B+D / 무처치 A)")
for act in ["active", "intermediate", "silent"]:
    d = s[s.is_CA & (s.activity == act)]
    lines.append("  " + binom_line(f"{act} all", d) + " | " + binom_line("B+D", d[d.group.isin(["B", "D"])])
                 + " | " + binom_line("A", d[d.group == "A"]))
lines.append("\n## 5. 로지스틱 회귀 (C>A만): U ~ condition + activity")
d = s[s.is_CA & s.activity.notna()].copy()
d["U"] = (d.ts == "U").astype(int)
try:
    m = smf.logit("U ~ C(group, Treatment('A')) + C(activity, Treatment('silent'))", d).fit(disp=0)
    lines.append(m.summary2().tables[1].round(3).to_string())
except Exception as e:  # 표본이 작으면 수렴 실패 가능
    lines.append(f"  fit failed: {e}")
out = "\n".join(lines)
open(f"{fd}/strand_asymmetry.summary.txt", "w").write(out + "\n")
print(out)
