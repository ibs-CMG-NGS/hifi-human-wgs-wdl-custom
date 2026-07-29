version 1.0

import "wdl-common/wdl/structs.wdl"
import "wdl-common/wdl/tasks/transgene_integration.wdl" as TG

# Transgene Integration Site Analysis Pipeline
#
# De novo assembly 결과(hap1/hap2 FASTA)와 트랜스진 서열로부터
# 게놈 내 트랜스진 통합 위치를 확인합니다.
# haplotagged_bam을 추가 제공하면 hybrid reference 기반 정밀 분석도 함께 수행합니다.
#
# 사용 방법:
#   miniwdl run workflows/transgene_integration.wdl \
#     --input <sample>.tg_integration.inputs.json \
#     --dir <output_dir> --verbose

workflow transgene_integration {
  meta {
    description: "Identify transgene integration sites from de novo assembly contigs, with optional hybrid-reference fine-mapping and gene annotation"
  }

  parameter_meta {
    sample_id:              { name: "Sample ID" }
    transgene_name:         { name: "Transgene construct name (e.g. gfa2_CreERT2)" }
    transgene_fasta:        { name: "Transgene sequence FASTA" }
    hap1_fasta:             { name: "Haplotype 1 assembly FASTA (from hifiasm)" }
    hap2_fasta:             { name: "Haplotype 2 assembly FASTA (from hifiasm)" }
    ref_fasta:              { name: "Reference genome FASTA (e.g. GRCm39)" }
    min_match_bp:           { name: "Min matching bases to call chimeric contig (default: 500)" }
    haplotagged_bam:        { name: "[Optional] HiFi reads aligned to reference (haplotagged BAM) — enables hybrid reference fine-mapping" }
    haplotagged_bai:        { name: "[Optional] BAI index for haplotagged_bam — required when haplotagged_bam is provided" }
    extract_window_bp:      { name: "Window (bp) around integration site for read extraction (default: 300000)" }
    annotation_gtf_bgz:     { name: "[Optional] Sorted bgzipped Gencode GTF — enables gene annotation" }
    annotation_gtf_tbi:     { name: "[Optional] Tabix index for annotation_gtf_bgz — required when annotation_gtf_bgz is provided" }
    annotation_window_bp:   { name: "Window (bp) around breakpoint for gene annotation query (default: 200000)" }
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

    # Optional: hybrid reference fine-mapping (needs haplotagged BAM)
    File?  haplotagged_bam
    File?  haplotagged_bai
    Int    extract_window_bp = 300000

    # Optional: gene annotation (needs sorted+indexed GTF)
    File?  annotation_gtf_bgz
    File?  annotation_gtf_tbi
    Int    annotation_window_bp = 200000

    RuntimeAttributes default_runtime_attributes
  }

  # ── Step 3-1: 트랜스진 → assembled contigs 매핑 (hap1, hap2 병렬) ────────────

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

  # ── Step 3-2: chimeric contig 추출 (hap1, hap2 병렬) ─────────────────────────

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

  # ── Step 4-1: chimeric contig → reference 정렬 (hap1, hap2 병렬) ─────────────

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

  # ── Step 4-2: 대략적 통합 위치 좌표 추출 (hybrid ref 분석 시작점) ─────────────

  call TG.extract_integration_coords {
    input:
      hap1_ref_paf       = ref_align_hap1.paf,
      hap2_ref_paf       = ref_align_hap2.paf,
      hap1_tg_paf        = map_hap1.paf,
      hap2_tg_paf        = map_hap2.paf,
      runtime_attributes = default_runtime_attributes,
  }

  # ── Step 5: Hybrid Reference 정밀 분석 (haplotagged_bam 제공 시) ──────────────

  if (defined(haplotagged_bam) && defined(haplotagged_bai)) {

    call TG.build_hybrid_ref {
      input:
        ref_fasta          = ref_fasta,
        transgene_fasta    = transgene_fasta,
        integration_chr    = extract_integration_coords.integration_chr,
        out_prefix         = sample_id,
        runtime_attributes = default_runtime_attributes,
    }

    call TG.extract_region_reads {
      input:
        haplotagged_bam    = select_first([haplotagged_bam]),
        haplotagged_bai    = select_first([haplotagged_bai]),
        integration_chr    = extract_integration_coords.integration_chr,
        integration_pos    = extract_integration_coords.integration_pos,
        window_bp          = extract_window_bp,
        out_prefix         = sample_id,
        runtime_attributes = default_runtime_attributes,
    }

    call TG.align_to_hybrid_ref {
      input:
        hybrid_ref_fa      = build_hybrid_ref.hybrid_ref_fa,
        hybrid_ref_fai     = build_hybrid_ref.hybrid_ref_fai,
        region_reads_fq    = extract_region_reads.region_reads_fq,
        out_prefix         = sample_id,
        runtime_attributes = default_runtime_attributes,
    }

    call TG.sort_index_bam {
      input:
        in_sam             = align_to_hybrid_ref.hybrid_sam,
        out_prefix         = sample_id,
        runtime_attributes = default_runtime_attributes,
    }

    call TG.detect_chimeric_reads {
      input:
        hybrid_paf         = align_to_hybrid_ref.hybrid_paf,
        tg_seqname         = build_hybrid_ref.tg_seqname,
        out_prefix         = sample_id,
        runtime_attributes = default_runtime_attributes,
    }

    # ── Step 5-5: 유전자 어노테이션 (annotation GTF 제공 시) ──────────────────

    if (defined(annotation_gtf_bgz) && defined(annotation_gtf_tbi)) {

      call TG.query_gene_annotation {
        input:
          annotation_gtf_bgz = select_first([annotation_gtf_bgz]),
          annotation_gtf_tbi = select_first([annotation_gtf_tbi]),
          breakpoint_chr     = detect_chimeric_reads.breakpoint_chr,
          breakpoint_pos     = detect_chimeric_reads.breakpoint_pos,
          window_bp          = annotation_window_bp,
          out_prefix         = sample_id,
          runtime_attributes = default_runtime_attributes,
      }

    }
  }

  # ── Step 6: 통합 리포트 생성 (항상 실행) ─────────────────────────────────────

  call TG.integration_report {
    input:
      sample_id              = sample_id,
      transgene_name         = transgene_name,
      hap1_tg_paf            = map_hap1.paf,
      hap2_tg_paf            = map_hap2.paf,
      hap1_ref_paf           = ref_align_hap1.paf,
      hap2_ref_paf           = ref_align_hap2.paf,
      out_prefix             = sample_id + "." + transgene_name,
      precise_breakpoint_chr = select_first([detect_chimeric_reads.breakpoint_chr,   "N/A"]),
      precise_breakpoint_pos = select_first([detect_chimeric_reads.breakpoint_pos,   "0"]),
      breakpoint_range       = select_first([detect_chimeric_reads.breakpoint_range, "N/A"]),
      chimeric_count         = select_first([detect_chimeric_reads.chimeric_count,   "0"]),
      bp_support             = select_first([detect_chimeric_reads.bp_support,       "N/A"]),
      tg_strand_vs_ref       = select_first([detect_chimeric_reads.tg_strand,        "N/A"]),
      nearest_gene_name      = select_first([query_gene_annotation.nearest_gene_name,      "N/A"]),
      nearest_gene_type      = select_first([query_gene_annotation.nearest_gene_type,      "N/A"]),
      insertion_feature_type = select_first([query_gene_annotation.insertion_feature_type, "N/A"]),
      nearest_exon_dist_bp   = select_first([query_gene_annotation.nearest_exon_dist_bp,   "N/A"]),
      runtime_attributes     = default_runtime_attributes,
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

    # 대략적 통합 위치 (de novo 기반)
    String stat_integration_chr = extract_integration_coords.integration_chr
    String stat_integration_pos = extract_integration_coords.integration_pos

    # hybrid reference 정밀 분석 BAM (IGV 시각화용)
    File?  hybrid_bam            = sort_index_bam.bam
    File?  hybrid_bai            = sort_index_bam.bai

    # hybrid reference 정밀 분석 결과 (선택적)
    File?  chimeric_reads_tsv   = detect_chimeric_reads.chimeric_tsv
    String stat_breakpoint_chr  = select_first([detect_chimeric_reads.breakpoint_chr,   "N/A"])
    String stat_breakpoint_pos  = select_first([detect_chimeric_reads.breakpoint_pos,   "0"])
    String stat_breakpoint_range= select_first([detect_chimeric_reads.breakpoint_range, "N/A"])
    String stat_chimeric_count  = select_first([detect_chimeric_reads.chimeric_count,   "0"])
    String stat_bp_support      = select_first([detect_chimeric_reads.bp_support,       "N/A"])
    String stat_tg_strand       = select_first([detect_chimeric_reads.tg_strand,        "N/A"])

    # 유전자 어노테이션 결과 (선택적)
    File?  annotation_tsv         = query_gene_annotation.annotation_tsv
    String stat_nearest_gene      = select_first([query_gene_annotation.nearest_gene_name,      "N/A"])
    String stat_insertion_feature = select_first([query_gene_annotation.insertion_feature_type, "N/A"])

    # 리포트
    File   report_txt           = integration_report.report_txt
    File   report_html          = integration_report.report_html
    File   integration_tsv      = integration_report.integration_tsv

    # 요약 통계 (de novo 기반)
    String stat_hap1_hit_count  = extract_hap1.hit_count
    String stat_hap2_hit_count  = extract_hap2.hit_count
  }
}
