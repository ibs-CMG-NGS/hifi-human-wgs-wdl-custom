version 1.0

import "wdl-common/wdl/structs.wdl"
import "wdl-common/wdl/tasks/transgene_integration.wdl" as TG

# Transgene Integration Site Analysis Pipeline
#
# De novo assembly 결과(hap1/hap2 FASTA)와 트랜스진 서열로부터
# 게놈 내 트랜스진 통합 위치를 확인합니다.
#
# 사용 방법:
#   miniwdl run workflows/transgene_integration.wdl \
#     --input <sample>.tg_integration.inputs.json \
#     --dir <output_dir> --verbose

workflow transgene_integration {
  meta {
    description: "Identify transgene integration sites from de novo assembly contigs"
  }

  parameter_meta {
    sample_id:           { name: "Sample ID" }
    transgene_name:      { name: "Transgene construct name (e.g. gfa2_CreERT2)" }
    transgene_fasta:     { name: "Transgene sequence FASTA" }
    hap1_fasta:          { name: "Haplotype 1 assembly FASTA (from hifiasm)" }
    hap2_fasta:          { name: "Haplotype 2 assembly FASTA (from hifiasm)" }
    ref_fasta:           { name: "Reference genome FASTA (e.g. GRCm39)" }
    min_match_bp:        { name: "Min matching bases to call chimeric contig (default: 500)" }
    default_runtime_attributes: { name: "Runtime attributes" }
  }

  input {
    String sample_id
    String transgene_name
    File   transgene_fasta
    File   hap1_fasta
    File   hap2_fasta
    File   ref_fasta
    Int    min_match_bp = 500

    RuntimeAttributes default_runtime_attributes
  }

  # Step 3-1: 트랜스진 → assembled contigs 매핑 (hap1, hap2 병렬)
  call TG.map_transgene_to_assembly as map_hap1 {
    input:
      transgene_fasta    = transgene_fasta,
      assembly_fasta     = hap1_fasta,
      out_prefix         = sample_id + ".hap1",
      runtime_attributes = default_runtime_attributes,
  }

  call TG.map_transgene_to_assembly as map_hap2 {
    input:
      transgene_fasta    = transgene_fasta,
      assembly_fasta     = hap2_fasta,
      out_prefix         = sample_id + ".hap2",
      runtime_attributes = default_runtime_attributes,
  }

  # Step 3-2: chimeric contig 추출 (hap1, hap2 병렬)
  call TG.extract_chimeric_contigs as extract_hap1 {
    input:
      tg_paf             = map_hap1.paf,
      assembly_fasta     = hap1_fasta,
      out_prefix         = sample_id + ".hap1",
      min_match_bp       = min_match_bp,
      runtime_attributes = default_runtime_attributes,
  }

  call TG.extract_chimeric_contigs as extract_hap2 {
    input:
      tg_paf             = map_hap2.paf,
      assembly_fasta     = hap2_fasta,
      out_prefix         = sample_id + ".hap2",
      min_match_bp       = min_match_bp,
      runtime_attributes = default_runtime_attributes,
  }

  # Step 4-1: chimeric contig → reference 정렬 (hap1, hap2 병렬)
  call TG.align_chimeric_to_ref as ref_align_hap1 {
    input:
      chimeric_fasta     = extract_hap1.chimeric_fasta,
      ref_fasta          = ref_fasta,
      out_prefix         = sample_id + ".hap1",
      runtime_attributes = default_runtime_attributes,
  }

  call TG.align_chimeric_to_ref as ref_align_hap2 {
    input:
      chimeric_fasta     = extract_hap2.chimeric_fasta,
      ref_fasta          = ref_fasta,
      out_prefix         = sample_id + ".hap2",
      runtime_attributes = default_runtime_attributes,
  }

  # Step 4-2: 통합 위치 리포트 생성
  call TG.integration_report {
    input:
      sample_id          = sample_id,
      transgene_name     = transgene_name,
      hap1_tg_paf        = map_hap1.paf,
      hap2_tg_paf        = map_hap2.paf,
      hap1_ref_paf       = ref_align_hap1.paf,
      hap2_ref_paf       = ref_align_hap2.paf,
      out_prefix         = sample_id + "." + transgene_name,
      runtime_attributes = default_runtime_attributes,
  }

  output {
    # 트랜스진 vs 어셈블리 PAF
    File   hap1_tg_paf          = map_hap1.paf
    File   hap2_tg_paf          = map_hap2.paf

    # chimeric contig 목록 및 서열
    File   hap1_chimeric_fasta  = extract_hap1.chimeric_fasta
    File   hap2_chimeric_fasta  = extract_hap2.chimeric_fasta
    File   hap1_chimeric_list   = extract_hap1.chimeric_contig_list
    File   hap2_chimeric_list   = extract_hap2.chimeric_contig_list

    # chimeric contig vs reference PAF (split alignment)
    File   hap1_ref_paf         = ref_align_hap1.paf
    File   hap2_ref_paf         = ref_align_hap2.paf

    # 리포트
    File   report_txt           = integration_report.report_txt
    File   integration_tsv      = integration_report.integration_tsv

    # 요약 통계
    String stat_hap1_hit_count  = extract_hap1.hit_count
    String stat_hap2_hit_count  = extract_hap2.hit_count
    String stat_integration_chr = integration_report.stat_integration_chr
    String stat_integration_pos = integration_report.stat_integration_pos
  }
}
