version 1.0

import "humanwgs_structs.wdl"
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "process_trgt_catalog/process_trgt_catalog.wdl" as ProcessTrgtCatalog
import "upstream/upstream.wdl" as Upstream
import "downstream/downstream.wdl" as Downstream
import "tertiary/tertiary.wdl" as TertiaryAnalysis
import "wdl-common/wdl/tasks/utilities.wdl" as Utilities
import "wdl-common/wdl/tasks/hifiasm.wdl"   as Hifiasm
import "wdl-common/wdl/tasks/busco.wdl"     as Busco
import "wdl-common/wdl/tasks/merqury.wdl"   as Merqury
import "wdl-common/wdl/tasks/asm_ref.wdl"   as AsmRef

workflow humanwgs_singleton {
  meta {
    description: "PacBio HiFi human whole genome sequencing pipeline for individual samples."
  }

  parameter_meta {
    sample_id: {
      name: "Unique identifier for the sample"
    }
    sex: {
      name: "Sample sex",
      choices: ["MALE", "FEMALE"]
    }
    hifi_reads: {
      name: "Array of paths to HiFi reads in unaligned BAM format."
    }
    fail_reads: {
      name: "Array of paths to failed reads in unaligned BAM format; optional"
    }
    phenotypes: {
      name: "Comma-delimited list of HPO codes for phenotypes"
    }
    ref_map_file: {
      name: "TSV containing reference genome file paths; must match backend"
    }
    tertiary_map_file: {
      name: "TSV containing tertiary analysis file paths and thresholds; must match backend"
    }
    max_reads_per_alignment_chunk: {
      name: "Maximum reads per alignment chunk"
    }
    pharmcat_min_coverage: {
      name: "Minimum coverage for PharmCAT"
    }
    gpu: {
      name: "Use GPU when possible"
    }
    backend: {
      name: "Backend where the workflow will be executed",
      choices: ["GCP", "Azure", "AWS-HealthOmics", "HPC"]
    }
    zones: {
      name: "Zones where compute will take place; required if backend is set to 'GCP'"
    }
    cpuPlatform: {
      help: "Optional minimum CPU platform to use for tasks on GCP"
    }
    gpuType: {
      name: "GPU type to use; required if gpu is set to `true` for cloud backends; must match backend"
    }
    container_registry: {
      name: "Container registry where workflow images are hosted. If left blank, PacBio's public Quay.io registry will be used. Must be set if backend is set to 'AWS-HealthOmics'",
      default: "quay.io/pacbio"
    }
    preemptible: {
      name: "Where possible, run tasks preemptibly"
    }
    debug_version: {
      name: "Debug version for testing purposes"
    }
    run_assembly: {
      name: "Run de novo assembly with hifiasm and QC (BUSCO, Merqury)"
    }
    busco_lineage: {
      name: "BUSCO lineage dataset for assembly completeness evaluation"
    }
  }

  input {
    String sample_id
    String? sex
    Array[File] hifi_reads
    Array[File]? fail_reads

    String phenotypes = "HP:0000001"

    File ref_map_file
    File? tertiary_map_file

    Int max_reads_per_alignment_chunk = 500000
    Int pharmcat_min_coverage = 10

    Boolean gpu = true
    Boolean run_pgx = true      # set to false for non-human references (e.g. mouse)
    Boolean run_assembly = false # run hifiasm de novo assembly + BUSCO + Merqury QC
    String  busco_lineage = "mammalia_odb10"  # BUSCO lineage; change per species

    # Backend configuration
    String backend
    String? zones
    String? cpuPlatform
    String? gpuType
    String? container_registry

    Boolean preemptible = true

    String? debug_version
  }

  call BackendConfiguration.backend_configuration {
    input:
      backend            = backend,
      zones              = zones,
      cpuPlatform        = cpuPlatform,
      gpuType            = gpuType,
      container_registry = container_registry
  }

  RuntimeAttributes default_runtime_attributes = if preemptible then backend_configuration.spot_runtime_attributes else backend_configuration.on_demand_runtime_attributes

  Map[String, String] ref_map = read_map(ref_map_file)

  call ProcessTrgtCatalog.process_trgt_catalog {
    input:
      trgt_catalog               = ref_map["trgt_tandem_repeat_bed"],  # !FileCoercion
      ref_fasta                  = ref_map["fasta"],                   # !FileCoercion
      ref_index                  = ref_map["fasta_index"],             # !FileCoercion
      default_runtime_attributes = default_runtime_attributes
  }

  call Upstream.upstream {
    input:
      sample_id                     = sample_id,
      sex                           = sex,
      hifi_reads                    = hifi_reads,
      fail_reads                    = fail_reads,
      ref_map_file                  = ref_map_file,
      max_reads_per_alignment_chunk = max_reads_per_alignment_chunk,
      trgt_catalog                  = process_trgt_catalog.full_catalog,
      fail_reads_bed                = process_trgt_catalog.include_fail_reads_bed,
      fail_reads_bait_fasta         = process_trgt_catalog.fail_reads_bait_fasta,
      fail_reads_bait_index         = process_trgt_catalog.fail_reads_bait_index,
      single_sample                 = true,
      gpu                           = gpu,
      default_runtime_attributes    = default_runtime_attributes
  }

  call Downstream.downstream {
    input:
      sample_id                  = sample_id,
      small_variant_vcf          = upstream.small_variant_vcf,
      small_variant_vcf_index    = upstream.small_variant_vcf_index,
      sv_vcf                     = select_first([upstream.sv_vcf]),
      sv_vcf_index               = select_first([upstream.sv_vcf_index]),
      trgt_vcf                   = upstream.trgt_vcf,
      trgt_vcf_index             = upstream.trgt_vcf_index,
      trgt_catalog               = process_trgt_catalog.full_catalog,
      aligned_bam                = upstream.out_bam,
      aligned_bam_index          = upstream.out_bam_index,
      pharmcat_min_coverage      = pharmcat_min_coverage,
      run_pgx                    = run_pgx,
      ref_map_file               = ref_map_file,
      default_runtime_attributes = default_runtime_attributes
  }

  if (run_assembly) {
    call Hifiasm.hifiasm {
      input:
        sample_id          = sample_id,
        hifi_reads         = hifi_reads,
        out_prefix         = sample_id + ".hifiasm",
        runtime_attributes = default_runtime_attributes
    }

    call Hifiasm.gfa_to_fasta as gfa_to_fasta_primary {
      input:
        gfa                = hifiasm.primary_gfa,
        out_prefix         = sample_id + ".hifiasm.bp.p_ctg",
        runtime_attributes = default_runtime_attributes
    }

    call Hifiasm.gfa_to_fasta as gfa_to_fasta_hap1 {
      input:
        gfa                = hifiasm.hap1_gfa,
        out_prefix         = sample_id + ".hifiasm.bp.hap1.p_ctg",
        runtime_attributes = default_runtime_attributes
    }

    call Hifiasm.gfa_to_fasta as gfa_to_fasta_hap2 {
      input:
        gfa                = hifiasm.hap2_gfa,
        out_prefix         = sample_id + ".hifiasm.bp.hap2.p_ctg",
        runtime_attributes = default_runtime_attributes
    }

    call Busco.busco as busco_primary {
      input:
        assembly_fasta     = gfa_to_fasta_primary.fasta,
        out_prefix         = sample_id + ".hifiasm.bp.p_ctg",
        lineage            = busco_lineage,
        runtime_attributes = default_runtime_attributes
    }

    call Busco.busco as busco_hap1 {
      input:
        assembly_fasta     = gfa_to_fasta_hap1.fasta,
        out_prefix         = sample_id + ".hifiasm.bp.hap1.p_ctg",
        lineage            = busco_lineage,
        runtime_attributes = default_runtime_attributes
    }

    call Busco.busco as busco_hap2 {
      input:
        assembly_fasta     = gfa_to_fasta_hap2.fasta,
        out_prefix         = sample_id + ".hifiasm.bp.hap2.p_ctg",
        lineage            = busco_lineage,
        runtime_attributes = default_runtime_attributes
    }

    call Merqury.meryl_count {
      input:
        hifi_reads         = hifi_reads,
        out_prefix         = sample_id + ".hifiasm",
        runtime_attributes = default_runtime_attributes
    }

    call Merqury.merqury {
      input:
        meryl_db           = meryl_count.meryl_db,
        hap1_fasta         = gfa_to_fasta_hap1.fasta,
        hap2_fasta         = gfa_to_fasta_hap2.fasta,
        out_prefix         = sample_id + ".hifiasm",
        runtime_attributes = default_runtime_attributes
    }

    # --- assembly-to-reference: minimap2 + SyRI (hap1, hap2) ---
    call AsmRef.minimap2_asm_to_ref as minimap2_hap1 {
      input:
        assembly_fasta     = gfa_to_fasta_hap1.fasta,
        ref_fasta          = ref_map["fasta"],           # !FileCoercion
        out_prefix         = sample_id + ".hifiasm.bp.hap1.p_ctg",
        runtime_attributes = default_runtime_attributes
    }

    call AsmRef.minimap2_asm_to_ref as minimap2_hap2 {
      input:
        assembly_fasta     = gfa_to_fasta_hap2.fasta,
        ref_fasta          = ref_map["fasta"],           # !FileCoercion
        out_prefix         = sample_id + ".hifiasm.bp.hap2.p_ctg",
        runtime_attributes = default_runtime_attributes
    }

    call AsmRef.syri as syri_hap1 {
      input:
        paf                = minimap2_hap1.paf,
        assembly_fasta     = gfa_to_fasta_hap1.fasta,
        ref_fasta          = ref_map["fasta"],           # !FileCoercion
        out_prefix         = sample_id + ".hifiasm.bp.hap1.p_ctg",
        runtime_attributes = default_runtime_attributes
    }

    call AsmRef.syri as syri_hap2 {
      input:
        paf                = minimap2_hap2.paf,
        assembly_fasta     = gfa_to_fasta_hap2.fasta,
        ref_fasta          = ref_map["fasta"],           # !FileCoercion
        out_prefix         = sample_id + ".hifiasm.bp.hap2.p_ctg",
        runtime_attributes = default_runtime_attributes
    }

    # --- assembly-to-reference: QUAST (primary) ---
    call AsmRef.quast {
      input:
        assembly_fasta     = gfa_to_fasta_primary.fasta,
        ref_fasta          = ref_map["fasta"],           # !FileCoercion
        out_prefix         = sample_id + ".hifiasm.bp.p_ctg",
        runtime_attributes = default_runtime_attributes
    }
  }

  Map[String, String] pedigree_sex = {
    "MALE": "1",
    "FEMALE": "2",
    "": "."
  }

  # write sample metadata similar to pedigree format
  # family_id, sample_id, father_id, mother_id, sex, affected
  Array[String] sample_metadata = [
    sample_id, sample_id,
    ".", ".",
    pedigree_sex[upstream.inferred_sex], "2"
  ]

  if (defined(tertiary_map_file)) {
    call TertiaryAnalysis.tertiary_analysis {
      input:
        sample_metadata            = [sample_metadata],
        phenotypes                 = phenotypes,
        is_trio_kid                = [false],
        is_duo_kid                 = [false],
        small_variant_vcf          = downstream.phased_small_variant_vcf,
        small_variant_vcf_index    = downstream.phased_small_variant_vcf_index,
        sv_vcf                     = downstream.phased_sv_vcf,
        sv_vcf_index               = downstream.phased_sv_vcf_index,
        ref_map_file               = ref_map_file,
        tertiary_map_file          = select_first([tertiary_map_file]),
        default_runtime_attributes = default_runtime_attributes
    }
  }

  Array[Array[String]] stats = [
    ['sample_id', sample_id],
    ['read_count', downstream.stat_read_count],
    ['read_length_mean', downstream.stat_read_length_mean],
    ['read_length_median', downstream.stat_read_length_median],
    ['read_length_n50', downstream.stat_read_length_n50],
    ['read_quality_mean', downstream.stat_read_quality_mean],
    ['read_quality_median', downstream.stat_read_quality_median],
    ['mapped_read_count', downstream.stat_mapped_read_count],
    ['mapped_read_percent', downstream.stat_mapped_read_percent],
    ['gap_compressed_identity_mean', downstream.stat_gap_compressed_identity_mean],
    ['gap_compressed_identity_median', downstream.stat_gap_compressed_identity_median],
    ['depth_mean', upstream.stat_depth_mean],
    ['inferred_sex', upstream.inferred_sex],
    ['stat_phased_basepairs', downstream.stat_phased_basepairs],
    ['phase_block_ng50', downstream.stat_phase_block_ng50],
    ['cpg_combined_count', downstream.stat_combined_cpg_count],
    ['cpg_hap1_count', downstream.stat_hap1_cpg_count],
    ['cpg_hap2_count', downstream.stat_hap2_cpg_count],
    ['methbat_methylated_count', downstream.stat_methbat_methylated_count],
    ['methbat_unmethylated_count', downstream.stat_methbat_unmethylated_count],
    ['methbat_asm_count', downstream.stat_methbat_asm_count],
    ['SNV_count', downstream.stat_SNV_count],
    ['TSTV_ratio', downstream.stat_TSTV_ratio],
    ['HETHOM_ratio', downstream.stat_HETHOM_ratio],
    ['INDEL_count', downstream.stat_INDEL_count],
    ['sv_DUP_count', downstream.stat_sv_DUP_count],
    ['sv_DEL_count', downstream.stat_sv_DEL_count],
    ['sv_INS_count', downstream.stat_sv_INS_count],
    ['sv_INV_count', downstream.stat_sv_INV_count],
    ['sv_SWAP_count', downstream.stat_sv_SWAP_count],
    ['sv_BND_count', downstream.stat_sv_BND_count],
    ['trgt_genotyped_count', upstream.stat_trgt_genotyped_count],
    ['trgt_uncalled_count', upstream.stat_trgt_uncalled_count],
    ['assembly_primary_size',  select_first([gfa_to_fasta_primary.stat_assembly_size, "0"])],
    ['assembly_primary_ctgs',  select_first([gfa_to_fasta_primary.stat_contig_count,  "0"])],
    ['assembly_hap1_size',     select_first([gfa_to_fasta_hap1.stat_assembly_size,    "0"])],
    ['assembly_hap1_ctgs',     select_first([gfa_to_fasta_hap1.stat_contig_count,     "0"])],
    ['assembly_hap2_size',     select_first([gfa_to_fasta_hap2.stat_assembly_size,    "0"])],
    ['assembly_hap2_ctgs',     select_first([gfa_to_fasta_hap2.stat_contig_count,     "0"])],
    ['busco_primary_complete', select_first([busco_primary.stat_complete,             "0"])],
    ['busco_primary_single',   select_first([busco_primary.stat_single,               "0"])],
    ['busco_primary_missing',  select_first([busco_primary.stat_missing,              "0"])],
    ['merqury_qv',             select_first([merqury.stat_qv,                         "0"])],
    ['merqury_completeness',   select_first([merqury.stat_completeness,               "0"])],
    ['quast_n50',             select_first([quast.stat_n50,             "0"])],
    ['quast_ng50',            select_first([quast.stat_ng50,            "0"])],
    ['quast_misassemblies',   select_first([quast.stat_misassemblies,   "0"])],
    ['quast_genome_fraction', select_first([quast.stat_genome_fraction, "0"])],
    ['syri_hap1_inv',         select_first([syri_hap1.stat_inv_count,   "0"])],
    ['syri_hap1_trans',       select_first([syri_hap1.stat_trans_count, "0"])],
    ['syri_hap1_dup',         select_first([syri_hap1.stat_dup_count,   "0"])],
    ['syri_hap2_inv',         select_first([syri_hap2.stat_inv_count,   "0"])],
    ['syri_hap2_trans',       select_first([syri_hap2.stat_trans_count, "0"])],
    ['syri_hap2_dup',         select_first([syri_hap2.stat_dup_count,   "0"])]
  ]

  call Utilities.consolidate_stats {
    input:
      id                 = sample_id,
      stats              = stats,
      msg_array          = flatten([process_trgt_catalog.msg, upstream.msg]),
      runtime_attributes = default_runtime_attributes
  }

  output {
    # consolidated stats
    File stats_file = consolidate_stats.output_tsv
    File msg_file   = consolidate_stats.messages

    # bam stats
    File   bam_statistics                      = downstream.bam_statistics
    File   read_length_plot                    = downstream.read_length_plot
    File?  read_quality_plot                   = downstream.read_quality_plot
    File   mapq_distribution_plot              = downstream.mapq_distribution_plot
    File   mg_distribution_plot                = downstream.mg_distribution_plot
    String stat_read_count                     = downstream.stat_read_count
    String stat_read_length_mean               = downstream.stat_read_length_mean
    String stat_read_length_median             = downstream.stat_read_length_median
    String stat_read_length_n50                = downstream.stat_read_length_n50
    String stat_read_quality_mean              = downstream.stat_read_quality_mean
    String stat_read_quality_median            = downstream.stat_read_quality_median
    String stat_mapped_read_count              = downstream.stat_mapped_read_count
    String stat_mapped_read_percent            = downstream.stat_mapped_read_percent
    String stat_gap_compressed_identity_mean   = downstream.stat_gap_compressed_identity_mean
    String stat_gap_compressed_identity_median = downstream.stat_gap_compressed_identity_median

    # merged, haplotagged alignments
    File   merged_haplotagged_bam       = downstream.merged_haplotagged_bam
    File   merged_haplotagged_bam_index = downstream.merged_haplotagged_bam_index

    # mosdepth outputs
    File   mosdepth_summary                 = upstream.mosdepth_summary
    File   mosdepth_region_bed              = upstream.mosdepth_region_bed
    File   mosdepth_region_bed_index        = upstream.mosdepth_region_bed_index
    File   mosdepth_depth_distribution_plot = upstream.mosdepth_depth_distribution_plot
    String stat_depth_mean                  = upstream.stat_depth_mean
    String inferred_sex                     = upstream.inferred_sex

    # phasing stats
    File   phase_stats           = downstream.phase_stats
    File   phase_blocks          = downstream.phase_blocks
    File   phase_haplotags       = downstream.phase_haplotags
    String stat_phased_basepairs = downstream.stat_phased_basepairs
    String stat_phase_block_ng50 = downstream.stat_phase_block_ng50

    # methylation outputs and profile
    File?   cpg_combined_bed                = downstream.cpg_combined_bed
    File?   cpg_combined_bed_index          = downstream.cpg_combined_bed_index
    File?   cpg_hap1_bed                    = downstream.cpg_hap1_bed
    File?   cpg_hap1_bed_index              = downstream.cpg_hap1_bed_index
    File?   cpg_hap2_bed                    = downstream.cpg_hap2_bed
    File?   cpg_hap2_bed_index              = downstream.cpg_hap2_bed_index
    File?   cpg_combined_bw                 = downstream.cpg_combined_bw
    File?   cpg_hap1_bw                     = downstream.cpg_hap1_bw
    File?   cpg_hap2_bw                     = downstream.cpg_hap2_bw
    String  stat_cpg_hap1_count             = downstream.stat_hap1_cpg_count
    String  stat_cpg_hap2_count             = downstream.stat_hap2_cpg_count
    String  stat_cpg_combined_count         = downstream.stat_combined_cpg_count
    File?   methbat_profile                 = downstream.methbat_profile
    String stat_methbat_methylated_count   = downstream.stat_methbat_methylated_count
    String stat_methbat_unmethylated_count = downstream.stat_methbat_unmethylated_count
    String stat_methbat_asm_count          = downstream.stat_methbat_asm_count

    # sv outputs
    File phased_sv_vcf                 = downstream.phased_sv_vcf
    File phased_sv_vcf_index           = downstream.phased_sv_vcf_index
    File sv_supporting_reads           = select_first([upstream.sv_supporting_reads])
    File sv_copynum_bedgraph           = select_first([upstream.sv_copynum_bedgraph])
    File sv_depth_bw                   = select_first([upstream.sv_depth_bw])
    File sv_gc_bias_corrected_depth_bw = select_first([upstream.sv_gc_bias_corrected_depth_bw])
    File sv_maf_bw                     = select_first([upstream.sv_maf_bw])
    File sv_copynum_summary            = select_first([upstream.sv_copynum_summary])

    # sv stats
    String stat_sv_DUP_count  = downstream.stat_sv_DUP_count
    String stat_sv_DEL_count  = downstream.stat_sv_DEL_count
    String stat_sv_INS_count  = downstream.stat_sv_INS_count
    String stat_sv_INV_count  = downstream.stat_sv_INV_count
    String stat_sv_SWAP_count = downstream.stat_sv_SWAP_count
    String stat_sv_BND_count  = downstream.stat_sv_BND_count

    # small variant outputs
    File phased_small_variant_vcf       = downstream.phased_small_variant_vcf
    File phased_small_variant_vcf_index = downstream.phased_small_variant_vcf_index
    File small_variant_gvcf             = upstream.small_variant_gvcf
    File small_variant_gvcf_index       = upstream.small_variant_gvcf_index

    # small variant stats
    File   small_variant_stats             = downstream.small_variant_stats
    File   bcftools_roh_out                = downstream.bcftools_roh_out
    File   bcftools_roh_bed                = downstream.bcftools_roh_bed
    String stat_small_variant_SNV_count    = downstream.stat_SNV_count
    String stat_small_variant_INDEL_count  = downstream.stat_INDEL_count
    String stat_small_variant_TSTV_ratio   = downstream.stat_TSTV_ratio
    String stat_small_variant_HETHOM_ratio = downstream.stat_HETHOM_ratio
    File   snv_distribution_plot           = downstream.snv_distribution_plot
    File   indel_distribution_plot         = downstream.indel_distribution_plot

    # trgt outputs
    File   phased_trgt_vcf           = downstream.phased_trgt_vcf
    File   phased_trgt_vcf_index     = downstream.phased_trgt_vcf_index
    File   trgt_spanning_reads       = upstream.trgt_spanning_reads
    File   trgt_spanning_reads_index = upstream.trgt_spanning_reads_index
    File   trgt_coverage_dropouts    = downstream.trgt_coverage_dropouts
    String stat_trgt_genotyped_count = upstream.stat_trgt_genotyped_count
    String stat_trgt_uncalled_count  = upstream.stat_trgt_uncalled_count

    # paraphase outputs
    File? paraphase_output_json         = upstream.paraphase_output_json
    File? paraphase_realigned_bam       = upstream.paraphase_realigned_bam
    File? paraphase_realigned_bam_index = upstream.paraphase_realigned_bam_index
    File? paraphase_vcfs                = upstream.paraphase_vcfs

    # per sample mitorsaw outputs
    File mitorsaw_vcf       = upstream.mitorsaw_vcf
    File mitorsaw_vcf_index = upstream.mitorsaw_vcf_index
    File mitorsaw_hap_stats = upstream.mitorsaw_hap_stats

    # PGx outputs
    File? pbstarphase_json        = downstream.pbstarphase_json
    File? pharmcat_match_json     = downstream.pharmcat_match_json
    File? pharmcat_phenotype_json = downstream.pharmcat_phenotype_json
    File? pharmcat_report_html    = downstream.pharmcat_report_html
    File? pharmcat_report_json    = downstream.pharmcat_report_json

    # tertiary analysis outputs
    File? tertiary_small_variant_filtered_vcf           = tertiary_analysis.small_variant_filtered_vcf
    File? tertiary_small_variant_filtered_vcf_index     = tertiary_analysis.small_variant_filtered_vcf_index
    File? tertiary_small_variant_filtered_tsv           = tertiary_analysis.small_variant_filtered_tsv
    File? tertiary_small_variant_compound_het_vcf       = tertiary_analysis.small_variant_compound_het_vcf
    File? tertiary_small_variant_compound_het_vcf_index = tertiary_analysis.small_variant_compound_het_vcf_index
    File? tertiary_small_variant_compound_het_tsv       = tertiary_analysis.small_variant_compound_het_tsv
    File? tertiary_sv_filtered_vcf                      = tertiary_analysis.sv_filtered_vcf
    File? tertiary_sv_filtered_vcf_index                = tertiary_analysis.sv_filtered_vcf_index
    File? tertiary_sv_filtered_tsv                      = tertiary_analysis.sv_filtered_tsv

    # de novo assembly outputs (only present when run_assembly = true)
    File? assembly_primary_gfa   = hifiasm.primary_gfa
    File? assembly_hap1_gfa      = hifiasm.hap1_gfa
    File? assembly_hap2_gfa      = hifiasm.hap2_gfa
    File? assembly_primary_fasta = gfa_to_fasta_primary.fasta
    File? assembly_hap1_fasta    = gfa_to_fasta_hap1.fasta
    File? assembly_hap2_fasta    = gfa_to_fasta_hap2.fasta
    File? assembly_busco_primary = busco_primary.short_summary
    File? assembly_busco_hap1    = busco_hap1.short_summary
    File? assembly_busco_hap2    = busco_hap2.short_summary
    File? assembly_merqury_qv    = merqury.qv_file
    File? assembly_merqury_stats = merqury.completeness_stats

    # assembly stats (0 when run_assembly = false)
    String stat_assembly_primary_size  = select_first([gfa_to_fasta_primary.stat_assembly_size, "0"])
    String stat_assembly_primary_ctgs  = select_first([gfa_to_fasta_primary.stat_contig_count,  "0"])
    String stat_assembly_hap1_size     = select_first([gfa_to_fasta_hap1.stat_assembly_size,    "0"])
    String stat_assembly_hap1_ctgs     = select_first([gfa_to_fasta_hap1.stat_contig_count,     "0"])
    String stat_assembly_hap2_size     = select_first([gfa_to_fasta_hap2.stat_assembly_size,    "0"])
    String stat_assembly_hap2_ctgs     = select_first([gfa_to_fasta_hap2.stat_contig_count,     "0"])
    String stat_busco_primary_complete = select_first([busco_primary.stat_complete,             "0"])
    String stat_busco_primary_single   = select_first([busco_primary.stat_single,               "0"])
    String stat_busco_primary_missing  = select_first([busco_primary.stat_missing,              "0"])
    String stat_merqury_qv             = select_first([merqury.stat_qv,                         "0"])
    String stat_merqury_completeness   = select_first([merqury.stat_completeness,               "0"])

    # assembly-to-reference outputs (only present when run_assembly = true)
    File? assembly_minimap2_hap1_paf  = minimap2_hap1.paf
    File? assembly_minimap2_hap2_paf  = minimap2_hap2.paf
    File? assembly_syri_hap1_out      = syri_hap1.syri_out
    File? assembly_syri_hap1_summary  = syri_hap1.syri_summary
    File? assembly_syri_hap2_out      = syri_hap2.syri_out
    File? assembly_syri_hap2_summary  = syri_hap2.syri_summary
    File? assembly_quast_report       = quast.report_tsv
    File? assembly_quast_html         = quast.report_html

    # assembly-to-reference stats (0 when run_assembly = false)
    String stat_quast_n50              = select_first([quast.stat_n50,             "0"])
    String stat_quast_ng50             = select_first([quast.stat_ng50,            "0"])
    String stat_quast_misassemblies    = select_first([quast.stat_misassemblies,   "0"])
    String stat_quast_genome_fraction  = select_first([quast.stat_genome_fraction, "0"])
    String stat_syri_hap1_inv          = select_first([syri_hap1.stat_inv_count,   "0"])
    String stat_syri_hap1_trans        = select_first([syri_hap1.stat_trans_count, "0"])
    String stat_syri_hap1_dup          = select_first([syri_hap1.stat_dup_count,   "0"])
    String stat_syri_hap2_inv          = select_first([syri_hap2.stat_inv_count,   "0"])
    String stat_syri_hap2_trans        = select_first([syri_hap2.stat_trans_count, "0"])
    String stat_syri_hap2_dup          = select_first([syri_hap2.stat_dup_count,   "0"])

    # qc messages
    Array[String] msg = flatten(
      [
        process_trgt_catalog.msg,
        upstream.msg
      ]
    )

    # workflow metadata
    String workflow_name    = "humanwgs_singleton"
    String workflow_version = "v3.1.1" + if defined(debug_version) then "~{"-" + debug_version}" else ""
  }
}
