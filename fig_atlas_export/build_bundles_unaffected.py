"""fig-atlas 요청서(h2o2_astro_hifi_wgs_figure_bundles)의 번들 중 Fig 3 클론성 지표 재정의와
무관한 것을 만든다: common/sample_sheet, common/reference, fig2b, fig2c, fig2d, fig3f, fig4b, fig4d.

규칙(요청서 §0·§1): UTF-8·LF·CSV, 열 이름 그대로, 반올림 없음, 분율은 0–1, 조건 비교는 상염색체만,
결과 없는 조합은 행을 두지 않음. 번들마다 inputs/data.csv(+추가 파일)와 metadata/metadata.yaml.
사용: python build_bundles_unaffected.py ROOT   (ROOT = .../raw_data/2026-h2o2-wgs)
"""
import gzip
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = sys.argv[1]
REPO = "/data_4tb/hifi-human-wgs-wdl-custom"
TA = "/data_4tb/shared/wgs-tertiary-analysis/h2o2_analysis_results"
TAC = f"{TA}/cohort"
WORK = "/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/fig_atlas_work"
FD = "/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/figure_data"
RES = f"{REPO}/hifi-wdl-resources/GRCm39"
AUTOS = [f"chr{i}" for i in range(1, 20)]
COND_DESC = {"A": "untreated", "B": "acute_H2O2", "C": "acute_H2O2_drug", "D": "chronic_H2O2_drug"}
COMMIT = os.popen(f"git -C {REPO} rev-parse --short HEAD").read().strip()

meta = pd.read_csv(f"{REPO}/base_mod_per_sample/samples.tsv", sep="\t")
GROUP = dict(zip(meta.sample_id, meta.group))


def write_bundle(rel, files, md):
    d = os.path.join(ROOT, rel)
    os.makedirs(f"{d}/inputs", exist_ok=True)
    os.makedirs(f"{d}/metadata", exist_ok=True)
    for name, df in files.items():
        p = f"{d}/inputs/{name}"
        if name.endswith(".gz"):
            df.to_csv(p, index=False, lineterminator="\n", float_format="%.17g", compression="gzip")
        else:
            df.to_csv(p, index=False, lineterminator="\n", float_format="%.17g")
    with open(f"{d}/metadata/metadata.yaml", "w") as f:
        yaml.safe_dump(md, f, sort_keys=False, allow_unicode=True, width=100)
    print(f"wrote {rel}: " + ", ".join(f"{k} ({len(v)} rows)" for k, v in files.items()))


# ---------------- common/sample_sheet ----------------
ss = meta[["sample_id", "group", "run"]].rename(columns={"group": "condition"}).copy()
ss["condition_desc"] = ss.condition.map(COND_DESC)
ss["well_id"] = ""
ss["cell_prep_id"] = ""
ss["movie_id"] = meta.raw_hifi_bam.str.extract(r"(m\d+_\d+_\d+_s\d)", expand=False)
ss["smrt_cell_slot"] = ss.movie_id.str.extract(r"_s(\d)$", expand=False).astype(int)
write_bundle("common/sample_sheet", {"data.csv": ss}, {
    "figure_id": "sample_sheet", "panel": "all", "plot_type": "sample_sheet",
    "data_source": f"base_mod_per_sample/samples.tsv@{COMMIT}",
    "source_files": [f"{REPO}/base_mod_per_sample/samples.tsv"],
    "definition": "12 samples = 4 conditions x 3 independent replicates, one replicate per sequencing run. "
                  "well_id: not recorded. cell_prep_id: left empty until the lab answers whether the four wells "
                  "were seeded from one shared cell preparation. Extra columns: movie_id (PacBio Revio movie) and "
                  "smrt_cell_slot (order of the SMRT cell within its run; condition A was always slot 1 and D slot 4, "
                  "so condition is confounded with cell order/sequencing time within each run).",
})

# ---------------- common/reference ----------------
fai = pd.read_csv(f"{RES}/mouse_GRCm39.fasta.fai", sep="\t", header=None, usecols=[0, 1], names=["chrom", "length"])
fai = fai[fai.chrom.isin(AUTOS + ["chrX", "chrY", "chrM"])]
os.makedirs(f"{ROOT}/common/reference", exist_ok=True)
fai.to_csv(f"{ROOT}/common/reference/chrom_sizes.tsv", sep="\t", index=False, lineterminator="\n")
versions = pd.DataFrame([
    ("sequencing", "PacBio Revio HiFi (on-instrument CCS)", "", "HiFi reads with kinetics (fi/fp/ri/rp)"),
    ("base_modification_onboard", "jasmine (on-instrument)", "2.3.0", "5mC (C+m) and 6mA (A+a/T-a) in HiFi BAM"),
    ("alignment", "pbmm2", "1.17.0", "GRCm39 (mm39); humanwgs_singleton WDL"),
    ("small_variants", "DeepVariant", "1.9.0", "PACBIO model; PASS records"),
    ("joint_genotyping", "GLnexus", "1.4.3", "per condition, 3 replicates"),
    ("structural_variants", "Sawfish", "2.2.1", "single-sample calls (adopted set) and joint calls"),
    ("phasing", "HiPhase", "1.5.0-8af3f1e", "small variants + SVs"),
    ("tandem_repeats", "TRGT", "4.0.0-429afcd", ""),
    ("coverage", "mosdepth", "0.3.9", ""),
    ("cpg_methylation", "pb-CpG-tools (aligned_bam_to_cpg_scores)", "3.0.0", "combined and haplotype-resolved pileups"),
    ("methylation_profile", "MethBat", "0.15.0-dd492c9", ""),
    ("base_modification_5hmC", "jasmine", "26.1.3", "re-called on unaligned raw HiFi BAM"),
    ("base_modification_pileup", "modkit", "0.6.3", "--filter-threshold 0.7; --cpg --combine-strands for 5mC/5hmC"),
    ("variant_tools", "bcftools", "1.23", ""),
    ("variant_annotation", "Ensembl VEP", "110", "mouse cache release 110, GRCm39; canonical transcripts"),
    ("sv_annotation", "svpack consequence", "quay.io/pacbio/svpack:latest", "Ensembl GRCm39 release 113 GFF3"),
    ("dmr", "DSS (Bioconductor)", "Bioconductor 3.18", "DMLtest/callDMR, pairwise (to be replaced by condition + run model)"),
    ("gene_model", "GENCODE", "vM36", "gencode.vM36.annotation.gtf"),
    ("gene_model_sv", "Ensembl", "113", "Mus_musculus.GRCm39.113.chr.gff3"),
    ("repeats", "UCSC RepeatMasker", "mm39 rmsk.txt.gz (downloaded 2026-10-07)", ""),
    ("repeats", "UCSC simpleRepeat (TRF)", "mm39 simpleRepeat.txt.gz (file dated 2020-07-28)", ""),
    ("cpg_islands", "UCSC cpgIslandExt", "mm39", ""),
    ("mutational_signatures", "COSMIC SBS", "v3.4, mm39-normalized (SigProfilerAssignment reference)", "SBS18, SBS36, SBS45 used"),
    ("gene_sets", "MGI GO annotation (GAF 2.2)", "generated 2026-05-21", "20 curated sets, direct annotations only (no GO closure)"),
    ("burden_scripts", "wgs-tertiary-analysis variant_burden.R / sv_burden.R / geneset_burden.R", "", "chi-square with simulated p (presence), Kruskal-Wallis (counts)"),
    ("export_scripts", "hifi-human-wgs-wdl-custom fig_atlas_export/, base_mod_per_sample/", COMMIT, ""),
], columns=["layer", "tool", "version", "parameters"])
versions.to_csv(f"{ROOT}/common/reference/versions.csv", index=False, lineterminator="\n")
print("wrote common/reference")

# ---------------- gene id map ----------------
gid = defaultdict(list)
with gzip.open(f"{RES}/gencode.vM36.annotation.gtf.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        x = line.split("\t", 9)
        if x[2] != "gene":
            continue
        name = re.search(r'gene_name "([^"]+)"', x[8]).group(1)
        g = re.search(r'gene_id "([^"]+)"', x[8]).group(1).split(".")[0]
        if g not in gid[name]:
            gid[name].append(g)

# ---------------- per-sample burden events (autosomes) ----------------
small_ev, sv_ev = [], []
for s in meta.sample_id:
    v = pd.read_csv(f"{TA}/small_variants/{s}.canonical_impacts.tsv", sep="\t", dtype=str)
    v = v[v.IMPACT.fillna("").str.contains("HIGH|MODERATE") & v.SYMBOL.notna() & (v.SYMBOL != ".")]
    small_ev.append(v.assign(sample_id=s, autosome=v.CHROM.isin(AUTOS))[["sample_id", "CHROM", "SYMBOL", "autosome"]])
    sv = pd.read_csv(f"{TA}/structural_variants/{s}.annotated_sv.tsv", sep="\t", comment=None, dtype=str)
    sv.columns = ["CHROM", "POS", "END", "SVTYPE", "SVLEN", "BCSQ", "GT"]
    sv = sv[pd.to_numeric(sv.SVLEN, errors="coerce").abs() <= 1_000_000]  # SVLEN "." (BND 등)은 R의 NA 비교처럼 제외
    for _, r in sv.iterrows():
        if r.BCSQ in (".", "") or pd.isna(r.BCSQ):
            continue
        for e in r.BCSQ.split(","):
            p = e.split("|")
            if len(p) >= 2 and p[0] in ("sv:cds", "sv:utr") and p[1] not in (".", ""):
                sv_ev.append((s, r.CHROM, p[1], r.CHROM in AUTOS))
small_ev = pd.concat(small_ev)
sv_ev = pd.DataFrame(sv_ev, columns=["sample_id", "CHROM", "SYMBOL", "autosome"])
print(f"small events: {len(small_ev)} ({(~small_ev.autosome).sum()} on sex chromosomes); "
      f"SV events: {len(sv_ev)} ({(~sv_ev.autosome).sum()} on sex chromosomes)")


def gene_bundle(stats_csv, ev, cls, panel, rel, test_note):
    st = pd.read_csv(stats_csv)
    ev = ev[ev.autosome]
    if cls == "sv_disrupting":
        ev = ev[~ev.SYMBOL.str.match(r"^(Vmn1r|Mup|Pramel|Ear|Gm)")]
    per = ev.groupby(["SYMBOL", "sample_id"]).size().rename("n_variants").reset_index()
    tested = set(st.gene)
    sexonly = sorted(tested - set(per.SYMBOL))
    d = pd.DataFrame({
        "gene_id": [";".join(gid.get(g, [])) for g in st.gene], "gene_symbol": st.gene, "variant_class": cls,
        "n_variants_total": [int(per[per.SYMBOL == g].n_variants.sum()) for g in st.gene],
        "estimate": st.odds_ratio_extreme, "pvalue": st.pvalue, "fdr": st.fdr,
        "group_max": st.group_max, "group_min": st.group_min, "cramers_v": st.cramers_v,
        **{f"n_samples_{g}": st[f"n_{g}"] for g in "ABCD"}})
    per_out = per[per.SYMBOL.isin(tested)].rename(columns={"SYMBOL": "gene_symbol"})
    write_bundle(rel, {"data.csv": d, "per_sample.csv": per_out}, {
        "figure_id": rel.split("/")[-1], "panel": panel, "plot_type": "gene_burden_rank",
        "data_source": "wgs-tertiary-analysis " + ("variant_burden.R" if cls.startswith("small") else "sv_burden.R"),
        "source_files": [stats_csv.replace(TAC, "TA"), "TA/../small_variants/*.canonical_impacts.tsv" if cls.startswith("small")
                         else "TA/../structural_variants/*.annotated_sv.tsv"],
        "definition": test_note,
        "notes": [f"Genes tested: {len(st)}. Genes whose qualifying events lie only on sex chromosomes "
                  f"(n_variants_total = 0 after the autosome filter): {', '.join(sexonly) if sexonly else 'none'}.",
                  "estimate = Haldane-Anscombe (+0.5) odds ratio of carrier proportions between the conditions with "
                  "the highest and lowest proportion (group_max vs group_min); it is not a fixed contrast.",
                  "pvalue/fdr are copied from the upstream table (chi-square with simulated p because cells < 5; BH)."],
    })
    return sexonly


note2b = ("Per gene, number of samples per condition carrying >= 1 VEP HIGH/MODERATE-impact small variant "
          "(canonical transcript, PASS, SNVs and indels) tested across the four conditions (4x2 chi-square, simulated p). "
          "per_sample.csv counts qualifying variants per gene and sample on autosomes.")
note2c = ("Per gene, number of samples per condition carrying >= 1 single-sample Sawfish SV with svpack consequence "
          "sv:cds or sv:utr and |SVLEN| <= 1 Mb, autosomes only, excluding multi-copy families (Vmn1r, Mup, Pramel, Ear, Gm); "
          "tested as in 2b.")
sex2b = gene_bundle(f"{TAC}/variant_burden_stats.csv", small_ev, "small_HIGH_MODERATE", "2b",
                    "fig2-mutational-burden/fig2b_gene_burden_small", note2b)
gene_bundle(f"{TAC}/sv_burden_stats.csv", sv_ev, "sv_disrupting", "2c", "fig2-mutational-burden/fig2c_gene_burden_sv", note2c)

# ---------------- fig2d pathway burden (recomputed on autosomes) ----------------
import json
GS = json.load(open("/data_4tb/shared/wgs-tertiary-analysis/resources/genesets_curated.json"))
GS_GROUP = {"dna_repair_ber": "DNA repair", "dna_repair_core": "DNA repair", "dna_repair_dsb_nhej": "DNA repair",
            "dna_repair_mmr": "DNA repair", "dna_repair_dsb_hr": "DNA repair", "dna_repair_ner": "DNA repair",
            "dna_damage_checkpoint": "DNA repair", "apoptosis": "cell death", "intrinsic_apoptosis": "cell death",
            "ferroptosis": "cell death", "autophagy": "cell death", "cellular_senescence": "stress response",
            "oxidative_stress_response": "stress response", "glutathione_redox": "stress response",
            "hypoxia_response": "stress response", "er_stress_upr": "stress response", "mitochondrion": "organelle",
            "inflammatory_response": "signalling", "calcium_signaling": "signalling", "gliogenesis_astrocyte": "astrocyte"}
rows, counts_rows = [], []
for cls, ev in (("small", small_ev[small_ev.autosome]), ("sv", sv_ev[sv_ev.autosome])):
    kw = []
    for gs, genes in GS.items():
        e = ev[ev.SYMBOL.isin(genes)]
        cnt = e.groupby("sample_id").size().reindex(meta.sample_id, fill_value=0)
        for s, c in cnt.items():
            counts_rows.append((gs, cls, s, int(c)))
        groups = [cnt[[s for s in meta.sample_id if GROUP[s] == g]].to_numpy() for g in "ABCD"]
        try:
            p = stats.kruskal(*groups).pvalue if len(set(cnt)) > 1 else np.nan
        except ValueError:
            p = np.nan
        kw.append((gs, p))
        tot = {g: groups[i].sum() for i, g in enumerate("ABCD")}
        for comp, (g1, g0) in {"B_vs_A": ("B", "A"), "C_vs_B": ("C", "B"), "D_vs_C": ("D", "C"), "D_vs_A": ("D", "A")}.items():
            k1, k0 = int(tot[g1]), int(tot[g0])
            if k1 + k0 == 0:
                continue  # 결과 없는 조합은 행 없음
            # 같은 샘플 수(3 vs 3)의 조건부 Poisson: k1 ~ Binom(k1+k0, RR/(1+RR))
            lo, hi = stats.binomtest(k1, k1 + k0).proportion_ci(0.95, method="exact")
            rr = k1 / k0 if k0 else np.inf
            rows.append(dict(gene_set=gs, gene_set_group=GS_GROUP.get(gs, "other"), variant_class=cls,
                             n_genes=len(genes), comparison=comp, n_events_1=k1, n_events_0=k0, estimate=rr,
                             ci_low=lo / (1 - lo) if lo < 1 else np.inf, ci_high=hi / (1 - hi) if hi < 1 else np.inf,
                             estimate_pvalue=stats.binomtest(k1, k1 + k0).pvalue, gs_key=(gs, cls)))
    kwd = pd.DataFrame(kw, columns=["gs", "p"])
    ok = kwd.p.notna()
    kwd.loc[ok, "fdr"] = stats.false_discovery_control(kwd.p[ok], method="bh")
    for r in rows:
        if r["gs_key"][1] == cls:
            m = kwd[kwd.gs == r["gs_key"][0]].iloc[0]
            r["pvalue"], r["fdr"] = m.p, m.fdr
d2 = pd.DataFrame(rows).drop(columns="gs_key")
d2 = d2[[c for c in ["gene_set", "gene_set_group", "variant_class", "n_genes", "comparison", "estimate", "ci_low", "ci_high",
                     "estimate_pvalue", "pvalue", "fdr", "n_events_1", "n_events_0"]]]
cnts = pd.DataFrame(counts_rows, columns=["gene_set", "variant_class", "sample_id", "n_events"])
write_bundle("fig2-mutational-burden/fig2d_pathway_burden", {"data.csv": d2, "per_sample.csv": cnts}, {
    "figure_id": "fig2d_pathway_burden", "panel": "2d", "plot_type": "effect_dots", "data_source":
    f"recomputed here (fig_atlas_export/build_bundles_unaffected.py@{COMMIT}) from wgs-tertiary-analysis event files",
    "source_files": ["TA/../small_variants/*.canonical_impacts.tsv", "TA/../structural_variants/*.annotated_sv.tsv",
                     "/data_4tb/shared/wgs-tertiary-analysis/resources/genesets_curated.json"],
    "definition": "Per sample, number of qualifying events (small: VEP HIGH/MODERATE, canonical; sv: svpack sv:cds/sv:utr, "
                  "|SVLEN| <= 1 Mb) in genes of each of 20 curated GO gene sets, autosomes only. "
                  "estimate = event rate ratio between the two conditions named in `comparison` (3 vs 3 samples), "
                  "conditional-Poisson exact (Clopper-Pearson) 95% CI and exact p in estimate_pvalue. "
                  "pvalue/fdr = Kruskal-Wallis on per-sample counts across the four conditions (as upstream), BH over gene sets "
                  "within variant_class.",
    "notes": ["Upstream geneset_burden.R did not exclude chrX/chrY; this export is autosome-only, so p/FDR differ "
              "slightly from geneset_burden_*_curated.csv.",
              "Rows with zero events in both compared conditions are absent; an infinite ci_high means zero events "
              "in the reference condition.",
              "Small-variant events include indels; indel calls are sensitive to read quality, which is lowest in "
              "condition A (see sample_metrics read_qv_median)."],
})

# ---------------- fig3f copy-number instability ----------------
raw = pd.read_csv(f"{TAC}/copynum_locus_raw.csv")
d3 = raw.rename(columns={"locus": "locus_id", "sample": "sample_id"})
d3["gene_symbol"] = d3.locus_id.str.replace(r"-[a-z]$", "", regex=True)
d3 = d3[["locus_id", "chrom", "start", "end", "gene_symbol", "sample_id", "copy_number"]]
det = pd.read_csv(f"{TAC}/copynum_instability_detail.csv")
rank = det.rename(columns={"locus": "locus_id", "group": "condition", "rank": "instability_rank"})[
    ["locus_id", "condition", "range", "instability_rank"]].rename(columns={"range": "copy_number_range"})
write_bundle("fig3-clonal-homogenization/fig3f_copy_number_instability", {"data.csv": d3, "rank.csv": rank}, {
    "figure_id": "fig3f_copy_number_instability", "panel": "3f", "plot_type": "heatmap",
    "data_source": "Sawfish 2.2.1 copy-number bedgraph; wgs-tertiary-analysis copynum_stability_check.py",
    "source_files": ["TA/copynum_locus_raw.csv", "TA/copynum_instability_detail.csv"],
    "definition": "Length-weighted modal integer copy number of the Sawfish depth-based segmentation over each of seven "
                  "multi-copy loci per sample; copy_number_range = max - min across the three replicates of a condition; "
                  "instability_rank = rank of the range among conditions within a locus (ties averaged; higher = less stable).",
    "notes": ["Vmn1r-a/b/c are adjacent sub-regions of one chr7 cluster and are not independent.",
              "Descriptive only; no test is reported."],
})

# ---------------- fig4b PCA ----------------
sc = pd.read_csv(f"{WORK}/cpg/pca_scores.tsv", sep="\t")
var = pd.read_csv(f"{WORK}/cpg/pca_variance.tsv", sep="\t")
write_bundle("fig4-methylation-5mc/fig4b_pca", {"data.csv": sc, "variance.csv": var}, {
    "figure_id": "fig4b_pca", "panel": "4b", "plot_type": "pca_scores",
    "data_source": f"pb-CpG-tools 3.0.0 combined CpG pileup; fig_atlas_export/cpg_global_pca.py@{COMMIT}",
    "source_files": ["{run}/out/cpg_combined_bed/*.cpg_pileup.combined.bed.gz (12 samples)"],
    "definition": "CpGs on chr1-19 with coverage >= 10 in all 12 samples (18,789,241). Methylation = pb-CpG-tools "
                  "modified fraction (mod_score/100). The 200,000 CpGs with the highest across-sample variance were "
                  "centred per CpG and decomposed by SVD; scores = U*S. variance_explained over all 12 components.",
})

# ---------------- fig4d replicate dispersion ----------------
pr = pd.read_csv(f"{FD}/cpg_dispersion.pairs.tsv", sep="\t")
a = pr[pr.min_cov == 10].set_index(["s1", "s2"])
b = pr[pr.min_cov == 20].set_index(["s1", "s2"])
d4 = pd.DataFrame({"sample_a": [i[0] for i in a.index], "sample_b": [i[1] for i in a.index],
                   "same_condition": (a.g1 == a.g2).to_numpy(), "same_run": (a.run1 == a.run2).to_numpy(),
                   "msd_raw": a.raw_sqdist.to_numpy(), "msd_corrected": a.corrected_sqdist.to_numpy(),
                   "msd_corrected_cov20": b.loc[a.index, "corrected_sqdist"].to_numpy(),
                   "n_cpg": a.n_cpg.to_numpy(), "n_cpg_cov20": b.loc[a.index, "n_cpg"].to_numpy(),
                   "mean_abs_diff": a.mean_abs_diff.to_numpy()})
d4["same_condition"] = d4.same_condition.map({True: "true", False: "false"})
d4["same_run"] = d4.same_run.map({True: "true", False: "false"})
write_bundle("fig4-methylation-5mc/fig4d_replicate_dispersion", {"data.csv": d4}, {
    "figure_id": "fig4d_replicate_dispersion", "panel": "4d", "plot_type": "heatmap",
    "data_source": f"pb-CpG-tools 3.0.0; base_mod_per_sample/cpg_replicate_dispersion.py@{COMMIT}",
    "source_files": ["FD/cpg_dispersion.pairs.tsv"],
    "definition": "For each of the 66 sample pairs, over CpGs on chr1-19 covered in all 12 samples: msd_raw = mean of "
                  "(m_a - m_b)^2; msd_corrected subtracts the binomial sampling variance m(1-m)/(n-1) of each sample "
                  "(n = CpG coverage); computed on CpGs >= 10x in all samples (n_cpg) and >= 20x (n_cpg_cov20).",
    "notes": ["Values are in methylation-fraction squared units (the plan quotes them x 1e-3).",
              "Condition A was sequenced at the lowest read quality; the coverage correction removes binomial noise "
              "but not systematic read-quality effects."],
})
print("sex-chromosome-only genes in 2b:", sex2b)
