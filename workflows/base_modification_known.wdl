version 1.0

import "wdl-common/wdl/structs.wdl"
import "wdl-common/wdl/tasks/base_modification.wdl" as BM
import "wdl-common/wdl/tasks/cpg_pileup.wdl" as CPG

# PacBio HiFi known base modification 파이프라인 (5mC / 5hmC / 6mA)
#
# base_modification_discovery.wdl에서 ipdSummary 기반 unknown-mod 발굴 경로를
# 제외한 버전. ipdSummary/kineticsTools는 Revio HiFi kinetics(fi/fp/ri/rp)를
# 읽지 못하므로(참조: memory/ipdsummary-hifi-incompatible), HiFi에서 과학적으로
# 유효한 known modification 경로만 산출한다.
#
# 흐름: align(kinetics 보존) → jasmine ML calling(→ MM/ML) → 그룹 병합 →
#       MM 태그 필터 → modkit pileup → 4-group known-mod 비교 + cpg QC
#
# align/jasmine/merge 태스크는 base_modification_discovery.wdl과 해시가 동일하여
# miniwdl call cache가 재사용된다.

workflow base_modification_known {
  meta {
    description: "HiFi known base modification (5mC/5hmC/6mA) via jasmine ML + modkit pileup — ipdSummary unknown-mod 경로 제외"
  }

  input {
    String        experiment_id

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

    Int    min_coverage = 25

    # 다운스트림 필터링/시각화 임계값 (본 분석에서 재조정 시 WDL 수정 없이 inputs.json만 변경)
    Int    known_mods_min_total_coverage   = 10
    Float  known_mods_min_5mC_delta        = 10.0
    Float  known_mods_min_5hmC_ratio_delta = 0.1

    RuntimeAttributes default_runtime_attributes
  }

  # ── Step 1: 정렬 (kinetics 보존) ─────────────────────────────────────────────
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

  # ── Step 2: jasmine ML 변형 calling (kinetics → MM/ML) ──────────────────────
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

  # ── Step 3: 그룹별 BAM 병합 ─────────────────────────────────────────────────
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

  # ── Step 4: MM 태그 필터 (modkit deadlock 회피) ─────────────────────────────
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

  # ── Step 5: modkit pileup (5mC / 5hmC / 6mA) ────────────────────────────────
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

  # ── Step 6: 4-group known modification 비교 ─────────────────────────────────
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

  # ── Step 6.5: differential TSV → Parquet (다운스트림 분석/시각화용) ─────────
  call BM.known_mods_to_parquet {
    input:
      differential_mods_tsv = compare_known_mods.differential_mods_tsv,
      out_prefix             = experiment_id,
      runtime_attributes     = default_runtime_attributes,
  }

  # ── Step 6.6: 필터링 + 시각화 (genome-wide 개요, rescue 산점도, 델타 분포) ──
  call BM.filter_and_visualize_known_mods {
    input:
      differential_mods_parquet = known_mods_to_parquet.differential_mods_parquet,
      ref_fai                    = ref_fai,
      out_prefix                  = experiment_id,
      min_total_coverage          = known_mods_min_total_coverage,
      min_5mC_delta                = known_mods_min_5mC_delta,
      min_5hmC_ratio_delta        = known_mods_min_5hmC_ratio_delta,
      runtime_attributes          = default_runtime_attributes,
  }

  # ── Step 7: cpg_pileup QC (5mC cross-validation) ────────────────────────────
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

  output {
    # modkit pileup — 그룹별 유형별 BED (bgzip 압축)
    File   ctrl_5mC_bed    = pileup_ctrl.fivemC_bed
    File   drug_a_5mC_bed  = pileup_drug_a.fivemC_bed
    File   drug_b_5mC_bed  = pileup_drug_b.fivemC_bed
    File   drug_ab_5mC_bed = pileup_drug_ab.fivemC_bed

    File   ctrl_5hmC_bed    = pileup_ctrl.fivehmC_bed
    File   drug_a_5hmC_bed  = pileup_drug_a.fivehmC_bed
    File   drug_b_5hmC_bed  = pileup_drug_b.fivehmC_bed
    File   drug_ab_5hmC_bed = pileup_drug_ab.fivehmC_bed

    File   ctrl_6mA_bed    = pileup_ctrl.sixmA_bed
    File   drug_a_6mA_bed  = pileup_drug_a.sixmA_bed
    File   drug_b_6mA_bed  = pileup_drug_b.sixmA_bed
    File   drug_ab_6mA_bed = pileup_drug_ab.sixmA_bed

    # known modification 4-group 비교
    File   differential_known_mods         = compare_known_mods.differential_mods_tsv
    File   differential_known_mods_parquet = known_mods_to_parquet.differential_mods_parquet

    # 필터링 + 시각화
    File   significant_known_mods_tsv = filter_and_visualize_known_mods.significant_known_mods_tsv
    File   known_mods_manhattan_5mC   = filter_and_visualize_known_mods.manhattan_5mC_delta_png
    File   known_mods_manhattan_5hmC  = filter_and_visualize_known_mods.manhattan_5hmC_ratio_png
    File   known_mods_rescue_scatter  = filter_and_visualize_known_mods.rescue_scatter_png
    File   known_mods_delta_hist      = filter_and_visualize_known_mods.delta_histograms_png
    File   known_mods_analysis_summary = filter_and_visualize_known_mods.analysis_summary_txt

    # cpg_pileup QC
    File?  ctrl_cpg_qc_bed   = cpg_qc_ctrl.combined_bed
    File?  drug_a_cpg_qc_bed = cpg_qc_drug_a.combined_bed
  }
}
