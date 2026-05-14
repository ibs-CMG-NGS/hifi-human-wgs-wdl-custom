version 1.0

import "../structs.wdl"

task hifiasm {
  meta {
    description: "De novo genome assembly from HiFi reads using hifiasm"
  }

  parameter_meta {
    sample_id:          { name: "Sample ID" }
    hifi_reads:         { name: "Array of HiFi reads in unaligned BAM format" }
    out_prefix:         { name: "Output prefix for hifiasm GFA files" }
    runtime_attributes: { name: "Runtime attribute structure" }
    primary_gfa:        { name: "Primary assembly GFA (bp.p_ctg.gfa)" }
    hap1_gfa:           { name: "Haplotype 1 assembly GFA (bp.hap1.p_ctg.gfa)" }
    hap2_gfa:           { name: "Haplotype 2 assembly GFA (bp.hap2.p_ctg.gfa)" }
  }

  input {
    String      sample_id
    Array[File] hifi_reads
    String      out_prefix = sample_id + ".hifiasm"

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 32
  Int mem_gb    = 64
  Int disk_size = ceil(size(hifi_reads, "GB") * 4 + 50)

  command <<<
    set -euo pipefail

    hifiasm --version

    hifiasm \
      -o ~{out_prefix} \
      -t ~{threads} \
      ~{sep=" " hifi_reads}
  >>>

  output {
    File primary_gfa = "~{out_prefix}.bp.p_ctg.gfa"
    File hap1_gfa    = "~{out_prefix}.bp.hap1.p_ctg.gfa"
    File hap2_gfa    = "~{out_prefix}.bp.hap2.p_ctg.gfa"
  }

  runtime {
    docker:                "quay.io/biocontainers/hifiasm:0.25.0--h5ca1c30_0"
    cpu:                   threads
    memory:                mem_gb + " GiB"
    disk:                  disk_size + " GB"
    disks:                 "local-disk " + disk_size + " HDD"
    preemptible:           runtime_attributes.preemptible_tries
    maxRetries:            runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones:                 runtime_attributes.zones
    cpuPlatform:           runtime_attributes.cpuPlatform
  }
}

task gfa_to_fasta {
  meta {
    description: "Convert GFA assembly graph to FASTA and compute basic assembly statistics"
  }

  parameter_meta {
    gfa:                { name: "Input GFA file" }
    out_prefix:         { name: "Output FASTA file prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    fasta:              { name: "Output FASTA file" }
    stat_assembly_size: { name: "Total assembly size in bp" }
    stat_contig_count:  { name: "Number of contigs" }
  }

  input {
    File   gfa
    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil(size(gfa, "GB") * 3 + 10)

  command <<<
    set -euo pipefail

    awk '/^S/{print ">"$2; print $3}' ~{gfa} > ~{out_prefix}.fa

    awk '/^>/{next} {total += length($0)} END {print total}' \
      ~{out_prefix}.fa > assembly_size.txt || echo "0" > assembly_size.txt

    grep -c "^>" ~{out_prefix}.fa > contig_count.txt || echo "0" > contig_count.txt
  >>>

  output {
    File   fasta              = "~{out_prefix}.fa"
    String stat_assembly_size = read_string("assembly_size.txt")
    String stat_contig_count  = read_string("contig_count.txt")
  }

  runtime {
    docker:                "quay.io/biocontainers/hifiasm:0.19.5--h5b5514e_1"
    cpu:                   threads
    memory:                mem_gb + " GiB"
    disk:                  disk_size + " GB"
    disks:                 "local-disk " + disk_size + " HDD"
    preemptible:           runtime_attributes.preemptible_tries
    maxRetries:            runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones:                 runtime_attributes.zones
    cpuPlatform:           runtime_attributes.cpuPlatform
  }
}

task bam_to_fastq {
  meta {
    description: "Convert lima-demultiplexed HiFi BAM to FASTQ for hifiasm compatibility"
  }

  parameter_meta {
    bam:                { name: "HiFi reads BAM (lima-demultiplexed)" }
    runtime_attributes: { name: "Runtime attribute structure" }
    fastq:              { name: "Output FASTQ file" }
  }

  input {
    File               bam
    RuntimeAttributes  runtime_attributes
  }

  Int threads   = 8
  Int mem_gb    = 8
  Int disk_size = ceil(size(bam, "GB") * 3 + 10)
  String fastq_name = basename(bam, ".bam") + ".fastq"

  command <<<
    set -euo pipefail
    samtools fastq -@ ~{threads} ~{bam} > ~{fastq_name}
  >>>

  output {
    File fastq = fastq_name
  }

  runtime {
    docker:                "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
    cpu:                   threads
    memory:                mem_gb + " GiB"
    disk:                  disk_size + " GB"
    disks:                 "local-disk " + disk_size + " HDD"
    preemptible:           runtime_attributes.preemptible_tries
    maxRetries:            runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones:                 runtime_attributes.zones
    cpuPlatform:           runtime_attributes.cpuPlatform
  }
}
