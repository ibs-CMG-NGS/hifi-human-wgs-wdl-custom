version 1.0

import "wdl-common/wdl/structs.wdl"
import "wdl-common/wdl/tasks/base_modification.wdl" as BM
import "wdl-common/wdl/tasks/cpg_pileup.wdl" as CPG

# PacBio HiFi Kinetics 기반 DNA 염기 변형 발굴 파이프라인
#
# 4개 실험군(Control, Drug A, Drug B, Drug A+B) × 3마리 HiFi BAM에서
# 알려진 변형(5mC, 5hmC, 6mA)과 unknown modification을 genome-wide로 탐색합니다.
#
# 사용 방법:
#   miniwdl run workflows/base_modification_discovery.wdl \
#     --input base_mod_discovery.inputs.json \
#     --dir output/base_mod_discovery --verbose

workflow base_modification_discovery {
  meta {
    description: "Genome-wide base modification discovery from PacBio HiFi kinetics data — 5mC/5hmC/6mA + unknown modifications via ipdSummary two-sample comparison and UMAP/HDBSCAN clustering"
  }

  parameter_meta {
    experiment_id:            { name: "Experiment identifier (used in report headers)" }
    control_bams:             { name: "Array[3] raw HiFi BAMs — Control group (ip/pw kinetics tags preserved)" }
    drug_a_bams:              { name: "Array[3] raw HiFi BAMs — Drug A group (ROS inducer)" }
    drug_b_bams:              { name: "Array[3] raw HiFi BAMs — Drug B group (antioxidant)" }
    drug_ab_bams:             { name: "Array[3] raw HiFi BAMs — Drug A+B combined group" }
    control_ids:              { name: "Array[3] sample IDs for control replicates" }
    drug_a_ids:               { name: "Array[3] sample IDs for Drug A replicates" }
    drug_b_ids:               { name: "Array[3] sample IDs for Drug B replicates" }
    drug_ab_ids:              { name: "Array[3] sample IDs for Drug A+B replicates" }
    ref_fasta:                { name: "GRCm39 reference FASTA" }
    ref_fai:                  { name: "FAI index for ref_fasta" }
    annotation_gtf_bgz:       { name: "Sorted bgzipped Gencode GTF (vM36)" }
    annotation_gtf_tbi:       { name: "Tabix index for annotation_gtf_bgz" }
    cpg_island_tsv:           { name: "CpG island TSV (cpgIslandExt.sorted.mm39.tsv)" }
    min_ipd_ratio:            { name: "IPD ratio cutoff for unknown candidate filtering (default 2.0)" }
    min_mod_score:            { name: "modificationScore cutoff (default 30)" }
    min_coverage:             { name: "Minimum site coverage (default 25)" }
    ipd_pvalue:               { name: "ipdSummary p-value cutoff (default 0.001)" }
    extract_window_bp:        { name: "Window bp for motif sequence extraction (default 250)" }
    annotation_window_bp:     { name: "TSS window for promoter classification (default 2000)" }
    min_cluster_size:         { name: "HDBSCAN min_cluster_size for clustering (default 50)" }
    run_hifimeth_qc:          { name: "Enable hifimeth 5mC QC cross-validation (default false)" }
    default_runtime_attributes: { name: "Runtime attributes" }
  }

  input {
    String        experiment_id

    # Raw HiFi BAMs with kinetics (ip/pw tags)
    Array[File]   control_bams
    Array[File]   drug_a_bams
    Array[File]   drug_b_bams
    Array[File]   drug_ab_bams

    Array[String] control_ids
    Array[String] drug_a_ids
    Array[String] drug_b_ids
    Array[String] drug_ab_ids

    File   ref_fasta
    File   ref_fai

    # Gene annotation + CpG islands
    File   annotation_gtf_bgz
    File   annotation_gtf_tbi
    File   cpg_island_tsv

    # Tuning parameters
    Float  min_ipd_ratio       = 2.0
    Int    min_mod_score       = 30
    Int    min_coverage        = 25
    Float  ipd_pvalue          = 0.001
    Int    extract_window_bp   = 250
    Int    annotation_window_bp = 2000
    Int    min_cluster_size    = 50
    Boolean run_hifimeth_qc   = false

    # ipdSummary reference window 제한 (예: "chr19"). 빈 문자열=전체 게놈.
    # 포유류 WGS 전체 스캔은 매우 느리므로 파일럿/검증 시 특정 염색체로 제한 권장.
    String ipd_reference_window = ""

    RuntimeAttributes default_runtime_attributes
  }

  # ── Step 1: 정렬 (kinetics 보존) — 12샘플 병렬 ─────────────────────────────

  scatter (i in range(length(control_bams))) {
    call BM.pbmm2_align_kinetics as align_ctrl {
      input:
        raw_hifi_bam       = control_bams[i],
        ref_fasta          = ref_fasta,
        ref_fai            = ref_fai,
        sample_id          = control_ids[i],
        runtime_attributes = default_runtime_attributes,
    }
  }

  scatter (i in range(length(drug_a_bams))) {
    call BM.pbmm2_align_kinetics as align_drug_a {
      input:
        raw_hifi_bam       = drug_a_bams[i],
        ref_fasta          = ref_fasta,
        ref_fai            = ref_fai,
        sample_id          = drug_a_ids[i],
        runtime_attributes = default_runtime_attributes,
    }
  }

  scatter (i in range(length(drug_b_bams))) {
    call BM.pbmm2_align_kinetics as align_drug_b {
      input:
        raw_hifi_bam       = drug_b_bams[i],
        ref_fasta          = ref_fasta,
        ref_fai            = ref_fai,
        sample_id          = drug_b_ids[i],
        runtime_attributes = default_runtime_attributes,
    }
  }

  scatter (i in range(length(drug_ab_bams))) {
    call BM.pbmm2_align_kinetics as align_drug_ab {
      input:
        raw_hifi_bam       = drug_ab_bams[i],
        ref_fasta          = ref_fasta,
        ref_fai            = ref_fai,
        sample_id          = drug_ab_ids[i],
        runtime_attributes = default_runtime_attributes,
    }
  }

  # ── Step 1': ipdSummary용 by-strand 정렬 (fi/fp/ri/rp → ip/pw → 정렬) ─────────
  # jasmine가 kinetics를 MM/ML로 소비하므로, ipdSummary에는 ccs-kinetics-bystrandify로
  # ip/pw 태그를 만든 뒤 정렬한 별도 BAM을 공급한다.
  scatter (i in range(length(control_bams))) {
    call BM.ccs_kinetics_bystrandify as bystrand_ctrl {
      input:
        raw_hifi_bam       = control_bams[i],
        sample_id          = control_ids[i],
        runtime_attributes = default_runtime_attributes,
    }
    call BM.pbmm2_align_kinetics as align_ctrl_bystrand {
      input:
        raw_hifi_bam       = bystrand_ctrl.bystrand_bam,
        ref_fasta          = ref_fasta,
        ref_fai            = ref_fai,
        sample_id          = control_ids[i] + ".bystrand",
        runtime_attributes = default_runtime_attributes,
    }
  }

  scatter (i in range(length(drug_a_bams))) {
    call BM.ccs_kinetics_bystrandify as bystrand_drug_a {
      input:
        raw_hifi_bam       = drug_a_bams[i],
        sample_id          = drug_a_ids[i],
        runtime_attributes = default_runtime_attributes,
    }
    call BM.pbmm2_align_kinetics as align_drug_a_bystrand {
      input:
        raw_hifi_bam       = bystrand_drug_a.bystrand_bam,
        ref_fasta          = ref_fasta,
        ref_fai            = ref_fai,
        sample_id          = drug_a_ids[i] + ".bystrand",
        runtime_attributes = default_runtime_attributes,
    }
  }

  scatter (i in range(length(drug_b_bams))) {
    call BM.ccs_kinetics_bystrandify as bystrand_drug_b {
      input:
        raw_hifi_bam       = drug_b_bams[i],
        sample_id          = drug_b_ids[i],
        runtime_attributes = default_runtime_attributes,
    }
    call BM.pbmm2_align_kinetics as align_drug_b_bystrand {
      input:
        raw_hifi_bam       = bystrand_drug_b.bystrand_bam,
        ref_fasta          = ref_fasta,
        ref_fai            = ref_fai,
        sample_id          = drug_b_ids[i] + ".bystrand",
        runtime_attributes = default_runtime_attributes,
    }
  }

  scatter (i in range(length(drug_ab_bams))) {
    call BM.ccs_kinetics_bystrandify as bystrand_drug_ab {
      input:
        raw_hifi_bam       = drug_ab_bams[i],
        sample_id          = drug_ab_ids[i],
        runtime_attributes = default_runtime_attributes,
    }
    call BM.pbmm2_align_kinetics as align_drug_ab_bystrand {
      input:
        raw_hifi_bam       = bystrand_drug_ab.bystrand_bam,
        ref_fasta          = ref_fasta,
        ref_fai            = ref_fai,
        sample_id          = drug_ab_ids[i] + ".bystrand",
        runtime_attributes = default_runtime_attributes,
    }
  }

  # ── Step 2a: jasmine ML 변형 calling — 12샘플 병렬 ──────────────────────────

  scatter (i in range(length(control_bams))) {
    call BM.jasmine_modification_calling as jasmine_ctrl {
      input:
        aligned_bam        = align_ctrl.aligned_bam[i],
        aligned_bai        = align_ctrl.aligned_bai[i],
        sample_id          = control_ids[i],
        runtime_attributes = default_runtime_attributes,
    }
  }

  scatter (i in range(length(drug_a_bams))) {
    call BM.jasmine_modification_calling as jasmine_drug_a {
      input:
        aligned_bam        = align_drug_a.aligned_bam[i],
        aligned_bai        = align_drug_a.aligned_bai[i],
        sample_id          = drug_a_ids[i],
        runtime_attributes = default_runtime_attributes,
    }
  }

  scatter (i in range(length(drug_b_bams))) {
    call BM.jasmine_modification_calling as jasmine_drug_b {
      input:
        aligned_bam        = align_drug_b.aligned_bam[i],
        aligned_bai        = align_drug_b.aligned_bai[i],
        sample_id          = drug_b_ids[i],
        runtime_attributes = default_runtime_attributes,
    }
  }

  scatter (i in range(length(drug_ab_bams))) {
    call BM.jasmine_modification_calling as jasmine_drug_ab {
      input:
        aligned_bam        = align_drug_ab.aligned_bam[i],
        aligned_bai        = align_drug_ab.aligned_bai[i],
        sample_id          = drug_ab_ids[i],
        runtime_attributes = default_runtime_attributes,
    }
  }

  # ── Step 2b: 그룹별 BAM 병합 (~90x) ─────────────────────────────────────────

  call BM.merge_group_bams as merge_ctrl {
    input:
      input_bams         = jasmine_ctrl.modified_bam,
      out_prefix         = "control",
      runtime_attributes = default_runtime_attributes,
  }

  call BM.merge_group_bams as merge_drug_a {
    input:
      input_bams         = jasmine_drug_a.modified_bam,
      out_prefix         = "drug_a",
      runtime_attributes = default_runtime_attributes,
  }

  call BM.merge_group_bams as merge_drug_b {
    input:
      input_bams         = jasmine_drug_b.modified_bam,
      out_prefix         = "drug_b",
      runtime_attributes = default_runtime_attributes,
  }

  call BM.merge_group_bams as merge_drug_ab {
    input:
      input_bams         = jasmine_drug_ab.modified_bam,
      out_prefix         = "drug_ab",
      runtime_attributes = default_runtime_attributes,
  }

  # ── Step 2b': ipdSummary용 by-strand 정렬 BAM 병합 (ip/pw 보존) ───────────────
  call BM.merge_group_bams as merge_ctrl_kin {
    input:
      input_bams         = align_ctrl_bystrand.aligned_bam,
      out_prefix         = "control.kinetics",
      runtime_attributes = default_runtime_attributes,
  }

  call BM.merge_group_bams as merge_drug_a_kin {
    input:
      input_bams         = align_drug_a_bystrand.aligned_bam,
      out_prefix         = "drug_a.kinetics",
      runtime_attributes = default_runtime_attributes,
  }

  call BM.merge_group_bams as merge_drug_b_kin {
    input:
      input_bams         = align_drug_b_bystrand.aligned_bam,
      out_prefix         = "drug_b.kinetics",
      runtime_attributes = default_runtime_attributes,
  }

  call BM.merge_group_bams as merge_drug_ab_kin {
    input:
      input_bams         = align_drug_ab_bystrand.aligned_bam,
      out_prefix         = "drug_ab.kinetics",
      runtime_attributes = default_runtime_attributes,
  }

  # ── Step 3a: ipdSummary — 3개 비교 (A-ctrl, B-ctrl, AB-A) 병렬 ──────────────

  call BM.ipd_summary as ipd_a_vs_ctrl {
    input:
      treated_bam        = merge_drug_a_kin.merged_bam,
      treated_bai        = merge_drug_a_kin.merged_bai,
      control_bam        = merge_ctrl_kin.merged_bam,
      control_bai        = merge_ctrl_kin.merged_bai,
      ref_fasta          = ref_fasta,
      comparison_name    = "drug_a_vs_ctrl",
      ipd_pvalue         = ipd_pvalue,
      min_coverage       = min_coverage,
      reference_window   = ipd_reference_window,
      runtime_attributes = default_runtime_attributes,
  }

  call BM.ipd_summary as ipd_b_vs_ctrl {
    input:
      treated_bam        = merge_drug_b_kin.merged_bam,
      treated_bai        = merge_drug_b_kin.merged_bai,
      control_bam        = merge_ctrl_kin.merged_bam,
      control_bai        = merge_ctrl_kin.merged_bai,
      ref_fasta          = ref_fasta,
      comparison_name    = "drug_b_vs_ctrl",
      ipd_pvalue         = ipd_pvalue,
      min_coverage       = min_coverage,
      reference_window   = ipd_reference_window,
      runtime_attributes = default_runtime_attributes,
  }

  call BM.ipd_summary as ipd_ab_vs_ctrl {
    input:
      treated_bam        = merge_drug_ab_kin.merged_bam,
      treated_bai        = merge_drug_ab_kin.merged_bai,
      control_bam        = merge_ctrl_kin.merged_bam,
      control_bai        = merge_ctrl_kin.merged_bai,
      ref_fasta          = ref_fasta,
      comparison_name    = "drug_ab_vs_ctrl",
      ipd_pvalue         = ipd_pvalue,
      min_coverage       = min_coverage,
      reference_window   = ipd_reference_window,
      runtime_attributes = default_runtime_attributes,
  }

  # Drug A+B vs Drug A (rescue — 가장 중요)
  call BM.ipd_summary as ipd_ab_vs_a {
    input:
      treated_bam        = merge_drug_ab_kin.merged_bam,
      treated_bai        = merge_drug_ab_kin.merged_bai,
      control_bam        = merge_drug_a_kin.merged_bam,
      control_bai        = merge_drug_a_kin.merged_bai,
      ref_fasta          = ref_fasta,
      comparison_name    = "drug_ab_vs_drug_a",
      ipd_pvalue         = ipd_pvalue,
      min_coverage       = min_coverage,
      reference_window   = ipd_reference_window,
      runtime_attributes = default_runtime_attributes,
  }

  # ── Step 3b: modkit pileup — 4개 그룹 병렬 ──────────────────────────────────
  # pileup 전 MM 태그 없는 read 제거 (modkit 0.6.3 deadlock 회피)

  call BM.filter_modbam as filter_ctrl {
    input:
      merged_bam         = merge_ctrl.merged_bam,
      merged_bai         = merge_ctrl.merged_bai,
      group_name         = "control",
      runtime_attributes = default_runtime_attributes,
  }

  call BM.filter_modbam as filter_drug_a {
    input:
      merged_bam         = merge_drug_a.merged_bam,
      merged_bai         = merge_drug_a.merged_bai,
      group_name         = "drug_a",
      runtime_attributes = default_runtime_attributes,
  }

  call BM.filter_modbam as filter_drug_b {
    input:
      merged_bam         = merge_drug_b.merged_bam,
      merged_bai         = merge_drug_b.merged_bai,
      group_name         = "drug_b",
      runtime_attributes = default_runtime_attributes,
  }

  call BM.filter_modbam as filter_drug_ab {
    input:
      merged_bam         = merge_drug_ab.merged_bam,
      merged_bai         = merge_drug_ab.merged_bai,
      group_name         = "drug_ab",
      runtime_attributes = default_runtime_attributes,
  }

  call BM.jasmine_pileup as pileup_ctrl {
    input:
      merged_bam         = filter_ctrl.filtered_bam,
      merged_bai         = filter_ctrl.filtered_bai,
      ref_fasta          = ref_fasta,
      ref_fai            = ref_fai,
      group_name         = "control",
      min_coverage       = min_coverage,
      runtime_attributes = default_runtime_attributes,
  }

  call BM.jasmine_pileup as pileup_drug_a {
    input:
      merged_bam         = filter_drug_a.filtered_bam,
      merged_bai         = filter_drug_a.filtered_bai,
      ref_fasta          = ref_fasta,
      ref_fai            = ref_fai,
      group_name         = "drug_a",
      min_coverage       = min_coverage,
      runtime_attributes = default_runtime_attributes,
  }

  call BM.jasmine_pileup as pileup_drug_b {
    input:
      merged_bam         = filter_drug_b.filtered_bam,
      merged_bai         = filter_drug_b.filtered_bai,
      ref_fasta          = ref_fasta,
      ref_fai            = ref_fai,
      group_name         = "drug_b",
      min_coverage       = min_coverage,
      runtime_attributes = default_runtime_attributes,
  }

  call BM.jasmine_pileup as pileup_drug_ab {
    input:
      merged_bam         = filter_drug_ab.filtered_bam,
      merged_bai         = filter_drug_ab.filtered_bai,
      ref_fasta          = ref_fasta,
      ref_fai            = ref_fai,
      group_name         = "drug_ab",
      min_coverage       = min_coverage,
      runtime_attributes = default_runtime_attributes,
  }

  # ── Step 3c: cpg_pileup (기존 task 재사용 — 5mC QC 비교용) ──────────────────

  call CPG.cpg_pileup as cpg_qc_ctrl {
    input:
      haplotagged_bam       = merge_ctrl.merged_bam,
      haplotagged_bam_index = merge_ctrl.merged_bai,
      ref_fasta             = ref_fasta,
      ref_index             = ref_fai,
      out_prefix            = "control.cpg_qc",
      runtime_attributes    = default_runtime_attributes,
  }

  call CPG.cpg_pileup as cpg_qc_drug_a {
    input:
      haplotagged_bam       = merge_drug_a.merged_bam,
      haplotagged_bam_index = merge_drug_a.merged_bai,
      ref_fasta             = ref_fasta,
      ref_index             = ref_fai,
      out_prefix            = "drug_a.cpg_qc",
      runtime_attributes    = default_runtime_attributes,
  }

  # ── Step 4a: unknown 후보 필터링 — 4개 비교 병렬 ─────────────────────────────

  call BM.filter_and_classify as filter_a_vs_ctrl {
    input:
      modifications_csv  = ipd_a_vs_ctrl.modifications_csv,
      ref_fasta          = ref_fasta,
      group_label        = "drug_a_vs_ctrl",
      min_ipd_ratio      = min_ipd_ratio,
      min_mod_score      = min_mod_score,
      min_coverage       = min_coverage,
      runtime_attributes = default_runtime_attributes,
  }

  call BM.filter_and_classify as filter_b_vs_ctrl {
    input:
      modifications_csv  = ipd_b_vs_ctrl.modifications_csv,
      ref_fasta          = ref_fasta,
      group_label        = "drug_b_vs_ctrl",
      min_ipd_ratio      = min_ipd_ratio,
      min_mod_score      = min_mod_score,
      min_coverage       = min_coverage,
      runtime_attributes = default_runtime_attributes,
  }

  call BM.filter_and_classify as filter_ab_vs_ctrl {
    input:
      modifications_csv  = ipd_ab_vs_ctrl.modifications_csv,
      ref_fasta          = ref_fasta,
      group_label        = "drug_ab_vs_ctrl",
      min_ipd_ratio      = min_ipd_ratio,
      min_mod_score      = min_mod_score,
      min_coverage       = min_coverage,
      runtime_attributes = default_runtime_attributes,
  }

  call BM.filter_and_classify as filter_ab_vs_a {
    input:
      modifications_csv  = ipd_ab_vs_a.modifications_csv,
      ref_fasta          = ref_fasta,
      group_label        = "drug_ab_vs_drug_a",
      min_ipd_ratio      = min_ipd_ratio,
      min_mod_score      = min_mod_score,
      min_coverage       = min_coverage,
      runtime_attributes = default_runtime_attributes,
  }

  # ── Step 4b: known modification 4-group 비교 ─────────────────────────────────

  call BM.compare_known_mods {
    input:
      group_names        = ["control", "drug_a", "drug_b", "drug_ab"],
      fivemC_beds        = [pileup_ctrl.fivemC_bed, pileup_drug_a.fivemC_bed,
                            pileup_drug_b.fivemC_bed, pileup_drug_ab.fivemC_bed],
      fivehmC_beds       = [pileup_ctrl.fivehmC_bed, pileup_drug_a.fivehmC_bed,
                            pileup_drug_b.fivehmC_bed, pileup_drug_ab.fivehmC_bed],
      sixmA_beds         = [pileup_ctrl.sixmA_bed, pileup_drug_a.sixmA_bed,
                            pileup_drug_b.sixmA_bed, pileup_drug_ab.sixmA_bed],
      out_prefix         = experiment_id,
      runtime_attributes = default_runtime_attributes,
  }

  # ── Step 4c: UMAP+HDBSCAN 클러스터링 ─────────────────────────────────────────

  call BM.cluster_unknown_modifications {
    input:
      unknown_tsvs       = [filter_a_vs_ctrl.unknown_candidates_tsv,
                            filter_b_vs_ctrl.unknown_candidates_tsv,
                            filter_ab_vs_ctrl.unknown_candidates_tsv,
                            filter_ab_vs_a.unknown_candidates_tsv],
      group_labels       = ["drug_a_vs_ctrl", "drug_b_vs_ctrl",
                            "drug_ab_vs_ctrl", "drug_ab_vs_drug_a"],
      out_prefix         = experiment_id,
      min_cluster_size   = min_cluster_size,
      runtime_attributes = default_runtime_attributes,
  }

  # ── Step 5: 4-group 통합 비교 + rescue ───────────────────────────────────────

  call BM.multigroup_comparison {
    input:
      unknown_tsvs          = [filter_a_vs_ctrl.unknown_candidates_tsv,
                               filter_b_vs_ctrl.unknown_candidates_tsv,
                               filter_ab_vs_ctrl.unknown_candidates_tsv,
                               filter_ab_vs_a.unknown_candidates_tsv],
      comparison_labels     = ["drug_a_vs_ctrl", "drug_b_vs_ctrl",
                               "drug_ab_vs_ctrl", "drug_ab_vs_drug_a"],
      differential_mods_tsv = compare_known_mods.differential_mods_tsv,
      clustered_tsv         = cluster_unknown_modifications.clustered_tsv,
      out_prefix            = experiment_id,
      min_ipd_ratio         = min_ipd_ratio,
      runtime_attributes    = default_runtime_attributes,
  }

  # ── Step 6: IPD ratio → BigWig (3개 비교 병렬) ────────────────────────────────

  call BM.ipd_to_bigwig as bw_a_vs_ctrl {
    input:
      modifications_csv  = ipd_a_vs_ctrl.modifications_csv,
      chrom_sizes        = ref_fai,
      comparison_name    = "drug_a_vs_ctrl",
      runtime_attributes = default_runtime_attributes,
  }

  call BM.ipd_to_bigwig as bw_b_vs_ctrl {
    input:
      modifications_csv  = ipd_b_vs_ctrl.modifications_csv,
      chrom_sizes        = ref_fai,
      comparison_name    = "drug_b_vs_ctrl",
      runtime_attributes = default_runtime_attributes,
  }

  call BM.ipd_to_bigwig as bw_ab_vs_a {
    input:
      modifications_csv  = ipd_ab_vs_a.modifications_csv,
      chrom_sizes        = ref_fai,
      comparison_name    = "drug_ab_vs_drug_a",
      runtime_attributes = default_runtime_attributes,
  }

  # ── Step 7: Motif analysis ────────────────────────────────────────────────────

  call BM.motif_analysis {
    input:
      rescue_candidates  = multigroup_comparison.rescue_candidates,
      ref_fasta          = ref_fasta,
      drug_a_candidates  = filter_a_vs_ctrl.unknown_candidates_tsv,
      out_prefix         = experiment_id,
      window_bp          = extract_window_bp,
      runtime_attributes = default_runtime_attributes,
  }

  # ── Step 8: 유전체 어노테이션 ────────────────────────────────────────────────

  call BM.annotate_modifications {
    input:
      categorized_tsv    = multigroup_comparison.categorized_tsv,
      annotation_gtf_bgz = annotation_gtf_bgz,
      annotation_gtf_tbi = annotation_gtf_tbi,
      cpg_island_tsv     = cpg_island_tsv,
      out_prefix         = experiment_id,
      window_bp          = annotation_window_bp,
      runtime_attributes = default_runtime_attributes,
  }

  # ── Step 9: Pathway enrichment ────────────────────────────────────────────────

  call BM.pathway_enrichment {
    input:
      annotated_tsv      = annotate_modifications.annotated_tsv,
      out_prefix         = experiment_id,
      runtime_attributes = default_runtime_attributes,
  }

  # ── Step 10: 통합 HTML 보고서 ─────────────────────────────────────────────────

  call BM.modification_report {
    input:
      experiment_id         = experiment_id,
      rescue_candidates     = multigroup_comparison.rescue_candidates,
      fivehmC_rescue_tsv    = multigroup_comparison.fivehmC_rescue_tsv,
      venn_summary          = multigroup_comparison.venn_summary,
      cluster_summary       = cluster_unknown_modifications.cluster_summary,
      annotated_tsv         = annotate_modifications.annotated_tsv,
      go_enrichment         = pathway_enrichment.go_enrichment,
      kegg_enrichment       = pathway_enrichment.kegg_enrichment,
      differential_mods_tsv = compare_known_mods.differential_mods_tsv,
      group_names           = ["control", "drug_a", "drug_b", "drug_ab"],
      run_hifimeth_qc       = run_hifimeth_qc,
      out_prefix            = experiment_id,
      runtime_attributes    = default_runtime_attributes,
  }

  output {
    # ── 정렬 BAM (kinetics 보존) ──────────────────────────────────────────
    Array[File]  ctrl_aligned_bams    = align_ctrl.aligned_bam
    Array[File]  drug_a_aligned_bams  = align_drug_a.aligned_bam
    Array[File]  drug_b_aligned_bams  = align_drug_b.aligned_bam
    Array[File]  drug_ab_aligned_bams = align_drug_ab.aligned_bam

    # ── Pooled BAM (jasmine MM/ML 업데이트) ───────────────────────────────
    File   ctrl_merged_bam     = merge_ctrl.merged_bam
    File   drug_a_merged_bam   = merge_drug_a.merged_bam
    File   drug_b_merged_bam   = merge_drug_b.merged_bam
    File   drug_ab_merged_bam  = merge_drug_ab.merged_bam

    # ── ipdSummary 결과 ───────────────────────────────────────────────────
    File   ipd_a_vs_ctrl_gff   = ipd_a_vs_ctrl.modifications_gff
    File   ipd_a_vs_ctrl_csv   = ipd_a_vs_ctrl.modifications_csv
    File   ipd_b_vs_ctrl_gff   = ipd_b_vs_ctrl.modifications_gff
    File   ipd_ab_vs_a_gff     = ipd_ab_vs_a.modifications_gff
    File   ipd_ab_vs_a_csv     = ipd_ab_vs_a.modifications_csv

    # ── modkit pileup (5mC / 5hmC / 6mA, bgzip 압축) ─────────────────────
    File   ctrl_5mC_bed        = pileup_ctrl.fivemC_bed
    File   drug_a_5hmC_bed     = pileup_drug_a.fivehmC_bed
    File   drug_ab_5hmC_bed    = pileup_drug_ab.fivehmC_bed

    # ── known modification 비교 ────────────────────────────────────────────
    File   differential_known_mods  = compare_known_mods.differential_mods_tsv

    # ── unknown modification 후보 ─────────────────────────────────────────
    File   drug_a_unknown_candidates = filter_a_vs_ctrl.unknown_candidates_tsv
    File   drug_a_ggg_flagged        = filter_a_vs_ctrl.ggg_flagged_tsv

    # ── 클러스터링 ────────────────────────────────────────────────────────
    File   clustered_modifications   = cluster_unknown_modifications.clustered_tsv
    File   cluster_summary           = cluster_unknown_modifications.cluster_summary
    File?  umap_plot                 = cluster_unknown_modifications.umap_plot

    # ── Rescue 분석 ───────────────────────────────────────────────────────
    File   rescue_candidates         = multigroup_comparison.rescue_candidates
    File   fivehmC_rescue            = multigroup_comparison.fivehmC_rescue_tsv
    File   categorized_modifications = multigroup_comparison.categorized_tsv
    File   venn_summary              = multigroup_comparison.venn_summary

    # ── BigWig ────────────────────────────────────────────────────────────
    File   bw_drug_a_vs_ctrl   = bw_a_vs_ctrl.ipdRatio_bw
    File   bw_drug_b_vs_ctrl   = bw_b_vs_ctrl.ipdRatio_bw
    File   bw_drug_ab_vs_a     = bw_ab_vs_a.ipdRatio_bw

    # ── Motif analysis ────────────────────────────────────────────────────
    File   meme_results        = motif_analysis.meme_tar
    File   streme_results      = motif_analysis.streme_tar

    # ── 어노테이션 + enrichment ────────────────────────────────────────────
    File   annotated_tsv       = annotate_modifications.annotated_tsv
    File   go_enrichment       = pathway_enrichment.go_enrichment
    File   kegg_enrichment     = pathway_enrichment.kegg_enrichment
    File?  go_enrichment_plot  = pathway_enrichment.go_plot
    File?  kegg_enrichment_plot = pathway_enrichment.kegg_plot

    # ── 보고서 ───────────────────────────────────────────────────────────
    File   report_html           = modification_report.report_html
    File   rescue_priority_list  = modification_report.rescue_priority_list

    # ── cpg_pileup QC (5mC cross-validation) ─────────────────────────────
    File?  ctrl_cpg_qc_bed     = cpg_qc_ctrl.combined_bed
    File?  drug_a_cpg_qc_bed   = cpg_qc_drug_a.combined_bed
  }
}
