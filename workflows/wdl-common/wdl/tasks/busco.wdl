version 1.0

import "../structs.wdl"

task busco {
  meta {
    description: "Assess assembly completeness using BUSCO conserved gene markers"
  }

  parameter_meta {
    assembly_fasta:     { name: "Assembly FASTA file" }
    out_prefix:         { name: "Output prefix (used in output file naming)" }
    lineage:            { name: "BUSCO lineage dataset (e.g. mammalia_odb10, eukaryota_odb10)" }
    runtime_attributes: { name: "Runtime attribute structure" }
    short_summary:      { name: "BUSCO short summary text file" }
    stat_complete:      { name: "% complete BUSCOs" }
    stat_single:        { name: "% complete and single-copy BUSCOs" }
    stat_duplicated:    { name: "% complete and duplicated BUSCOs" }
    stat_fragmented:    { name: "% fragmented BUSCOs" }
    stat_missing:       { name: "% missing BUSCOs" }
  }

  input {
    File   assembly_fasta
    String out_prefix
    String lineage = "mammalia_odb10"

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 8
  Int mem_gb    = 16
  Int disk_size = ceil(size(assembly_fasta, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    busco --version

    busco \
      -m genome \
      -i ~{assembly_fasta} \
      -o busco_out \
      -l ~{lineage} \
      -c ~{threads} \
      --offline 2>/dev/null || \
    busco \
      -m genome \
      -i ~{assembly_fasta} \
      -o busco_out \
      -l ~{lineage} \
      -c ~{threads}

    cp busco_out/short_summary*.txt ~{out_prefix}.busco.short_summary.txt

    python3 << 'EOF'
    import re

    with open("~{out_prefix}.busco.short_summary.txt") as f:
        text = f.read()

    def get_pct(pattern):
        m = re.search(pattern, text)
        return m.group(1) if m else "0.0"

    complete   = get_pct(r'C:(\d+\.\d+)%')
    single     = get_pct(r'S:(\d+\.\d+)%')
    duplicated = get_pct(r'D:(\d+\.\d+)%')
    fragmented = get_pct(r'F:(\d+\.\d+)%')
    missing    = get_pct(r'M:(\d+\.\d+)%')

    for name, val in [
        ("complete",   complete),
        ("single",     single),
        ("duplicated", duplicated),
        ("fragmented", fragmented),
        ("missing",    missing),
    ]:
        with open(f"busco_{name}.txt", "w") as f:
            f.write(val)
    EOF
  >>>

  output {
    File   short_summary   = "~{out_prefix}.busco.short_summary.txt"
    String stat_complete   = read_string("busco_complete.txt")
    String stat_single     = read_string("busco_single.txt")
    String stat_duplicated = read_string("busco_duplicated.txt")
    String stat_fragmented = read_string("busco_fragmented.txt")
    String stat_missing    = read_string("busco_missing.txt")
  }

  runtime {
    docker:                "quay.io/biocontainers/busco:5.7.1--pyhdfd78af_1"
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
