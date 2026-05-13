version 1.0

import "../structs.wdl"

task minimap2_asm_to_ref {
  meta {
    description: "Align de novo assembly to reference genome using minimap2 (asm5 preset, PAF output for SyRI)"
  }

  parameter_meta {
    assembly_fasta:     { name: "Assembly FASTA (hap1 or hap2)" }
    ref_fasta:          { name: "Reference genome FASTA" }
    out_prefix:         { name: "Output file prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    paf:                { name: "PAF alignment file with --cs tag" }
  }

  input {
    File   assembly_fasta
    File   ref_fasta
    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 8
  Int mem_gb    = 16
  Int disk_size = ceil((size(assembly_fasta, "GB") + size(ref_fasta, "GB")) * 3 + 20)

  command <<<
    set -euo pipefail

    minimap2 --version

    minimap2 \
      -cx asm5 \
      --cs \
      -t ~{threads} \
      ~{ref_fasta} \
      ~{assembly_fasta} \
      > ~{out_prefix}.asm_to_ref.paf
  >>>

  output {
    File paf = "~{out_prefix}.asm_to_ref.paf"
  }

  runtime {
    docker:                "quay.io/biocontainers/minimap2:2.30--h577a1d6_0"
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

task syri {
  meta {
    description: "Detect synteny and rearrangements (inversions, translocations, duplications) by comparing assembly to reference"
  }

  parameter_meta {
    paf:                { name: "PAF alignment file with --cs tag (from minimap2 asm5)" }
    assembly_fasta:     { name: "Query assembly FASTA (same assembly used for PAF)" }
    ref_fasta:          { name: "Reference genome FASTA" }
    out_prefix:         { name: "Output file prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    syri_out:           { name: "SyRI variant output file" }
    syri_summary:       { name: "SyRI summary counts file" }
    stat_inv_count:     { name: "Number of inversions" }
    stat_trans_count:   { name: "Number of translocations" }
    stat_dup_count:     { name: "Number of duplications" }
  }

  input {
    File   paf
    File   assembly_fasta
    File   ref_fasta
    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 4
  Int mem_gb    = 8
  Int disk_size = ceil((size(paf, "GB") + size(assembly_fasta, "GB") + size(ref_fasta, "GB")) * 2 + 10)

  command <<<
    set -euo pipefail

    syri --version 2>&1 || true

    syri \
      -c ~{paf} \
      -r ~{ref_fasta} \
      -q ~{assembly_fasta} \
      -F P \
      --prefix ~{out_prefix}. \
      --nc ~{threads} \
      --lf ~{out_prefix}.syri.log

    python3 << 'EOF'
    import os

    inv_count   = "0"
    trans_count = "0"
    dup_count   = "0"

    summary_file = "~{out_prefix}.syri.summary"
    if os.path.exists(summary_file):
        with open(summary_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith("INV"):
                    parts = line.split()
                    inv_count = parts[1] if len(parts) > 1 else "0"
                elif line.startswith("TRANS"):
                    parts = line.split()
                    trans_count = parts[1] if len(parts) > 1 else "0"
                elif line.startswith("DUP"):
                    parts = line.split()
                    dup_count = parts[1] if len(parts) > 1 else "0"

    with open("syri_inv.txt", "w")   as f: f.write(inv_count)
    with open("syri_trans.txt", "w") as f: f.write(trans_count)
    with open("syri_dup.txt", "w")   as f: f.write(dup_count)
    EOF
  >>>

  output {
    File   syri_out         = "~{out_prefix}.syri.out"
    File   syri_summary     = "~{out_prefix}.syri.summary"
    String stat_inv_count   = read_string("syri_inv.txt")
    String stat_trans_count = read_string("syri_trans.txt")
    String stat_dup_count   = read_string("syri_dup.txt")
  }

  runtime {
    docker:                "quay.io/biocontainers/syri:1.7.1--py311hcf77733_1"
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

task quast {
  meta {
    description: "Evaluate assembly quality against reference genome using QUAST (N50, NG50, misassemblies, genome fraction)"
  }

  parameter_meta {
    assembly_fasta:       { name: "Primary assembly FASTA" }
    ref_fasta:            { name: "Reference genome FASTA" }
    out_prefix:           { name: "Output file prefix" }
    runtime_attributes:   { name: "Runtime attribute structure" }
    report_tsv:           { name: "QUAST tabular report" }
    report_html:          { name: "QUAST interactive HTML report" }
    stat_n50:             { name: "Assembly N50 (bp)" }
    stat_ng50:            { name: "Assembly NG50 relative to reference (bp)" }
    stat_misassemblies:   { name: "Number of misassemblies" }
    stat_genome_fraction: { name: "Genome fraction covered (%)" }
  }

  input {
    File   assembly_fasta
    File   ref_fasta
    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 8
  Int mem_gb    = 16
  Int disk_size = ceil((size(assembly_fasta, "GB") + size(ref_fasta, "GB")) * 3 + 20)

  command <<<
    set -euo pipefail

    quast.py --version

    quast.py \
      -r ~{ref_fasta} \
      -t ~{threads} \
      -o quast_out \
      --no-icarus \
      ~{assembly_fasta}

    cp quast_out/report.tsv  ~{out_prefix}.quast.report.tsv
    cp quast_out/report.html ~{out_prefix}.quast.report.html

    python3 << 'EOF'
    fields = {
        "N50":                 "quast_n50.txt",
        "NG50":                "quast_ng50.txt",
        "# misassemblies":     "quast_misassemblies.txt",
        "Genome fraction (%)": "quast_genome_fraction.txt",
    }
    values = {k: "0" for k in fields}

    with open("~{out_prefix}.quast.report.tsv") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] in values:
                values[parts[0]] = parts[1]

    for field, fname in fields.items():
        with open(fname, "w") as f:
            f.write(values[field])
    EOF
  >>>

  output {
    File   report_tsv           = "~{out_prefix}.quast.report.tsv"
    File   report_html          = "~{out_prefix}.quast.report.html"
    String stat_n50             = read_string("quast_n50.txt")
    String stat_ng50            = read_string("quast_ng50.txt")
    String stat_misassemblies   = read_string("quast_misassemblies.txt")
    String stat_genome_fraction = read_string("quast_genome_fraction.txt")
  }

  runtime {
    docker:                "quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2"
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
