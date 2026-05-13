version 1.0

import "../structs.wdl"

task meryl_count {
  meta {
    description: "Build k-mer database from HiFi reads using meryl (prerequisite for Merqury)"
  }

  parameter_meta {
    hifi_reads:         { name: "Array of HiFi reads in unaligned BAM format" }
    out_prefix:         { name: "Output meryl database prefix" }
    kmer_size:          { name: "K-mer size (default: 21)" }
    runtime_attributes: { name: "Runtime attribute structure" }
    meryl_db:           { name: "Meryl k-mer database (tar.gz)" }
  }

  input {
    Array[File] hifi_reads
    String      out_prefix
    Int         kmer_size = 21

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 16
  Int mem_gb    = 32
  Int disk_size = ceil(size(hifi_reads, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    meryl --version 2>&1 | head -1

    meryl count \
      k=~{kmer_size} \
      threads=~{threads} \
      memory=~{mem_gb} \
      ~{sep=" " hifi_reads} \
      output ~{out_prefix}.meryl

    tar -czf ~{out_prefix}.meryl.tar.gz ~{out_prefix}.meryl
  >>>

  output {
    File meryl_db = "~{out_prefix}.meryl.tar.gz"
  }

  runtime {
    docker:                "quay.io/biocontainers/merqury:1.3--hdfd78af_4"
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

task merqury {
  meta {
    description: "Evaluate assembly accuracy and completeness using k-mer spectrum (Merqury)"
  }

  parameter_meta {
    meryl_db:           { name: "Meryl k-mer database (tar.gz)" }
    hap1_fasta:         { name: "Haplotype 1 assembly FASTA" }
    hap2_fasta:         { name: "Haplotype 2 assembly FASTA" }
    out_prefix:         { name: "Output prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    qv_file:            { name: "QV scores TSV" }
    completeness_stats: { name: "K-mer completeness statistics" }
    stat_qv:            { name: "Assembly QV score" }
    stat_completeness:  { name: "K-mer completeness (%)" }
  }

  input {
    File   meryl_db
    File   hap1_fasta
    File   hap2_fasta
    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 8
  Int mem_gb    = 16
  Int disk_size = ceil((size(hap1_fasta, "GB") + size(hap2_fasta, "GB")) * 3 + size(meryl_db, "GB") + 20)

  command <<<
    set -euo pipefail

    tar -xzf ~{meryl_db}
    DB_NAME=$(basename ~{meryl_db} .tar.gz)

    merqury.sh "${DB_NAME}" ~{hap1_fasta} ~{hap2_fasta} ~{out_prefix}

    mv ~{out_prefix}.qv ~{out_prefix}.merqury.qv 2>/dev/null || true
    mv ~{out_prefix}.completeness.stats ~{out_prefix}.merqury.completeness.stats 2>/dev/null || true

    python3 << 'EOF'
    import os

    qv_val   = "0.0"
    comp_val = "0.0"

    qv_file   = "~{out_prefix}.merqury.qv"
    comp_file = "~{out_prefix}.merqury.completeness.stats"

    if os.path.exists(qv_file):
        lines = [l for l in open(qv_file) if not l.startswith('#') and l.strip()]
        if lines:
            parts = lines[-1].split()
            qv_val = parts[3] if len(parts) > 3 else "0.0"

    if os.path.exists(comp_file):
        lines = [l for l in open(comp_file) if not l.startswith('#') and l.strip()]
        if lines:
            parts = lines[-1].split()
            comp_val = parts[2] if len(parts) > 2 else "0.0"

    with open("merqury_qv.txt", "w") as f:
        f.write(qv_val)
    with open("merqury_completeness.txt", "w") as f:
        f.write(comp_val)
    EOF
  >>>

  output {
    File   qv_file            = "~{out_prefix}.merqury.qv"
    File   completeness_stats = "~{out_prefix}.merqury.completeness.stats"
    String stat_qv            = read_string("merqury_qv.txt")
    String stat_completeness  = read_string("merqury_completeness.txt")
  }

  runtime {
    docker:                "quay.io/biocontainers/merqury:1.3--hdfd78af_4"
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
