"""Fig 2e IGV 확인용 C>A 후보 목록과 IGV batch 파일 생성.

선택: 샘플별 private C>A 중 VAF ≥ 0.2, alt read 많은 순 상위 N개 (조건 비교를 위해 A 포함 12샘플 전부)
표시: SBS18 주요 문맥(SBS18 가중치 상위 4: T[C>A]T, G[C>A]A, C[C>A]A, T[C>A]A) 여부
비교 BAM: 같은 run의 다른 3개 샘플 (해당 자리에 alt read가 없어야 함)
사용: python igv_ca_candidates.py samples.tsv FIGDATA_DIR N [VAF_MIN VAF_MAX SBS18_ONLY(0/1) SUFFIX]
  예) 기본(판독 확인용):          ... 10
      산화 손상 후보(저 VAF):     ... 15 0.1 0.3 1 lowVAF_SBS18
"""
import sys

import pandas as pd

meta = pd.read_csv(sys.argv[1], sep="\t")
fd, n = sys.argv[2], int(sys.argv[3])
vmin, vmax = (float(sys.argv[4]), float(sys.argv[5])) if len(sys.argv) > 5 else (0.2, 1.01)
sbs18_only = len(sys.argv) > 6 and sys.argv[6] == "1"
suffix = f"_{sys.argv[7]}" if len(sys.argv) > 7 else ""
SBS18_TOP = {"T[C>A]T", "G[C>A]A", "C[C>A]A", "T[C>A]A"}
LABEL = {"A": "untreated", "B": "acute_H2O2", "C": "acute_H2O2+drug", "D": "chronic_H2O2+drug"}

# 아티팩트 의심 표시(high_dp: depth > 샘플 평균의 1.5배, clustered: 같은 샘플 후보가 10 kb 안에 있음)가
# 있는 파일을 쓰고, 저 VAF 목록(SUFFIX 지정 시)에서는 표시된 자리를 제외
s = pd.read_csv(f"{fd}/mutspec.private_sites.flagged.tsv.gz", sep="\t")
s = s[(s.class6 == "C>A") & (s.vaf >= vmin) & (s.vaf < vmax)]
if sbs18_only:
    s = s[s.ctx96.isin(SBS18_TOP)]
if suffix:
    s = s[~s.high_dp & ~s.clustered]
s = s.sort_values(["sample", "alt_reads", "vaf"], ascending=[True, False, False]).groupby("sample").head(n)
bam = dict(zip(meta.sample_id, meta.wgs_haplotagged_bam))
run = dict(zip(meta.sample_id, meta.run))
s["condition"] = s.group.map(LABEL)
s["run"] = s["sample"].map(run)
s["sbs18_top_context"] = s.ctx96.isin(SBS18_TOP)
s["locus"] = s.chrom + ":" + (s.pos - 50).astype(str) + "-" + (s.pos + 50).astype(str)
s["sample_bam"] = s["sample"].map(bam)
s["same_run_other_bams"] = [";".join(bam[o] for o in meta.sample_id if run[o] == run[x] and o != x) for x in s["sample"]]
order = {"B": 0, "D": 1, "C": 2, "A": 3}
s = s.sort_values(["group", "run", "alt_reads"], key=lambda c: c.map(order) if c.name == "group" else c,
                  ascending=[True, True, False])
cols = ["condition", "group", "run", "sample", "chrom", "pos", "ref", "alt", "ctx96", "sbs18_top_context",
        "vaf", "dp", "rel_dp", "high_dp", "clustered", "alt_reads", "locus", "sample_bam", "same_run_other_bams"]
s[cols].to_csv(f"{fd}/igv_CA_candidates{suffix}.tsv", sep="\t", index=False)

# IGV batch: run별로 4개 BAM을 한 번에 띄우고 해당 run 후보를 차례로 스냅샷
lines = ["new", "genome mm39", "preference SAM.SHOW_SOFT_CLIPPED false", "maxPanelHeight 400"]
for r in sorted(s.run.unique()):
    lines += ["new", "genome mm39"]
    for x in meta.sample_id[meta.run == r]:
        lines.append(f"load {bam[x]}")
    lines.append(f"snapshotDirectory {fd}/igv_snapshots{suffix}/run{r}")
    for _, row in s[s.run == r].iterrows():
        lines += [f"goto {row.locus}", "sort base", f"snapshot {row.group}_{row['sample']}_{row.chrom}_{row.pos}_{row.ctx96.replace('[', '').replace(']', '').replace('>', 'to')}.png"]
open(f"{fd}/igv_CA_candidates{suffix}.batch.txt", "w").write("\n".join(lines) + "\n")
print(s.groupby(["group", "condition"]).agg(n=("pos", "size"), sbs18_ctx=("sbs18_top_context", "sum"),
                                            median_vaf=("vaf", "median"), median_alt=("alt_reads", "median")))
print(f"\nwrote {fd}/igv_CA_candidates{suffix}.tsv and .batch.txt ({len(s)} sites)")
