version 1.0

import "../structs.wdl"

task map_transgene_to_assembly {
  meta {
    description: "Map transgene FASTA to assembled contigs to identify integration-containing (chimeric) contigs"
  }

  parameter_meta {
    transgene_fasta:    { name: "Transgene sequence FASTA" }
    assembly_fasta:     { name: "Haplotype assembly FASTA (hap1 or hap2)" }
    out_prefix:         { name: "Output prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    paf:                { name: "PAF: transgene vs assembly contigs" }
  }

  input {
    File   transgene_fasta
    File   assembly_fasta
    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 8
  Int mem_gb    = 24
  Int disk_size = ceil(size(assembly_fasta, "GB") * 3 + 5)

  command <<<
    set -euo pipefail
    minimap2 --version

    # asm20: 더 발산된 서열도 포착 (일반 transgene 대비 host 삽입 변이 허용)
    minimap2 \
      -cx asm20 \
      --cs \
      -t ~{threads} \
      ~{assembly_fasta} \
      ~{transgene_fasta} \
      > ~{out_prefix}.tg_vs_asm.paf
  >>>

  output {
    File paf = "~{out_prefix}.tg_vs_asm.paf"
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

task extract_chimeric_contigs {
  meta {
    description: "Identify and extract contigs containing transgene sequence from assembly"
  }

  parameter_meta {
    tg_paf:             { name: "PAF from map_transgene_to_assembly" }
    assembly_fasta:     { name: "Haplotype assembly FASTA" }
    out_prefix:         { name: "Output prefix" }
    min_match_bp:       { name: "Minimum matching bases to call a chimeric contig (default: 500)" }
    runtime_attributes: { name: "Runtime attribute structure" }
    chimeric_fasta:     { name: "FASTA of contigs containing transgene" }
    chimeric_contig_list: { name: "List of chimeric contig names" }
    hit_count:          { name: "Number of chimeric contigs found" }
  }

  input {
    File   tg_paf
    File   assembly_fasta
    String out_prefix
    Int    min_match_bp = 500

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 8
  Int disk_size = ceil(size(assembly_fasta, "GB") * 2 + 5)

  command <<<
    set -euo pipefail

    # PAF 컬럼: query=transgene, target=contig
    # mapq 필터 없이 match 길이만으로 판단
    # (같은 contig에 다중 카피 → mapq 낮아지는 정상 현상)
    awk '$10 >= ~{min_match_bp} {print $6}' ~{tg_paf} \
      | sort -u > ~{out_prefix}.chimeric_contigs.txt

    COUNT=$(wc -l < ~{out_prefix}.chimeric_contigs.txt)
    echo "Chimeric contigs found: ${COUNT}"
    echo "${COUNT}" > hit_count.txt

    if [ "${COUNT}" -eq 0 ]; then
      echo "WARNING: No chimeric contigs found. Check transgene FASTA and assembly." >&2
      # 빈 FASTA 생성 (파이프라인 계속 진행)
      touch ~{out_prefix}.chimeric_contigs.fa
    else
      samtools faidx ~{assembly_fasta}
      xargs samtools faidx ~{assembly_fasta} \
        < ~{out_prefix}.chimeric_contigs.txt \
        > ~{out_prefix}.chimeric_contigs.fa
    fi

    # PAF 내용 요약 출력
    echo "=== Transgene hit summary ==="
    awk '$10 >= ~{min_match_bp} {
      tg_cov = ($4-$3)/$2*100
      printf "  contig=%s len=%d tg_pos=%d-%d/%d (%.1f%%) strand=%s matches=%d\n",
        $6, $7, $3, $4, $2, tg_cov, $5, $10
    }' ~{tg_paf} | sort -k2 -rn
  >>>

  output {
    File   chimeric_fasta       = "~{out_prefix}.chimeric_contigs.fa"
    File   chimeric_contig_list = "~{out_prefix}.chimeric_contigs.txt"
    String hit_count            = read_string("hit_count.txt")
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

task align_chimeric_to_ref {
  meta {
    description: "Align chimeric contigs to reference genome to identify integration breakpoints"
  }

  parameter_meta {
    chimeric_fasta:     { name: "FASTA of chimeric contigs" }
    ref_fasta:          { name: "Reference genome FASTA" }
    out_prefix:         { name: "Output prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    paf:                { name: "PAF: chimeric contigs vs reference (with --eqx for split alignment)" }
  }

  input {
    File   chimeric_fasta
    File   ref_fasta
    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 8
  Int mem_gb    = 16
  Int disk_size = ceil((size(chimeric_fasta, "GB") + size(ref_fasta, "GB")) * 3 + 10)

  command <<<
    set -euo pipefail

    if [ ! -s ~{chimeric_fasta} ]; then
      echo "No chimeric contigs to align." >&2
      touch ~{out_prefix}.chimeric_to_ref.paf
    else
      minimap2 --version
      minimap2 \
        -cx asm5 \
        --cs \
        --eqx \
        -t ~{threads} \
        ~{ref_fasta} \
        ~{chimeric_fasta} \
        > ~{out_prefix}.chimeric_to_ref.paf
    fi
  >>>

  output {
    File paf = "~{out_prefix}.chimeric_to_ref.paf"
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

task integration_report {
  meta {
    description: "Parse alignment results and generate transgene integration site report"
  }

  parameter_meta {
    sample_id:            { name: "Sample ID" }
    transgene_name:       { name: "Transgene construct name" }
    hap1_tg_paf:          { name: "PAF: transgene vs hap1 contigs" }
    hap2_tg_paf:          { name: "PAF: transgene vs hap2 contigs" }
    hap1_ref_paf:         { name: "PAF: hap1 chimeric contigs vs reference" }
    hap2_ref_paf:         { name: "PAF: hap2 chimeric contigs vs reference" }
    out_prefix:           { name: "Output prefix" }
    runtime_attributes:   { name: "Runtime attribute structure" }
    report_txt:           { name: "Integration site report (text)" }
    integration_tsv:      { name: "Integration sites table (TSV)" }
    stat_hap1_hit_count:  { name: "Number of hap1 chimeric contigs" }
    stat_hap2_hit_count:  { name: "Number of hap2 chimeric contigs" }
    stat_integration_chr: { name: "Integration chromosome" }
    stat_integration_pos: { name: "Integration position (estimated)" }
  }

  input {
    String sample_id
    String transgene_name
    File   hap1_tg_paf
    File   hap2_tg_paf
    File   hap1_ref_paf
    File   hap2_ref_paf
    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = 10

  command <<<
    set -euo pipefail

    python3 << 'PYEOF'
import os
from collections import defaultdict

sample_id      = "~{sample_id}"
transgene_name = "~{transgene_name}"
out_prefix     = "~{out_prefix}"

def parse_tg_paf(paf_file, min_match=500):
    """트랜스진 vs assembly PAF에서 hit contig 정보 추출"""
    hits = defaultdict(list)
    if not os.path.exists(paf_file) or os.path.getsize(paf_file) == 0:
        return hits
    with open(paf_file) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) < 12:
                continue
            tg_name, tg_len, tg_start, tg_end = p[0], int(p[1]), int(p[2]), int(p[3])
            strand = p[4]
            contig, contig_len = p[5], int(p[6])
            c_start, c_end = int(p[7]), int(p[8])
            matches, mapq = int(p[9]), int(p[11])
            if matches >= min_match:
                hits[contig].append({
                    "tg_start": tg_start, "tg_end": tg_end, "tg_len": tg_len,
                    "contig_len": contig_len, "c_start": c_start, "c_end": c_end,
                    "strand": strand, "matches": matches, "mapq": mapq,
                })
    return hits

def parse_ref_paf(paf_file, min_match=1000):
    """chimeric contig vs reference PAF에서 정렬 블록 추출"""
    alns = defaultdict(list)
    if not os.path.exists(paf_file) or os.path.getsize(paf_file) == 0:
        return alns
    with open(paf_file) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) < 12:
                continue
            contig, q_len = p[0], int(p[1])
            q_start, q_end = int(p[2]), int(p[3])
            strand = p[4]
            chrom = p[5]
            t_start, t_end = int(p[7]), int(p[8])
            matches, mapq = int(p[9]), int(p[11])
            if matches >= min_match:
                alns[contig].append({
                    "q_len": q_len, "q_start": q_start, "q_end": q_end,
                    "chrom": chrom, "t_start": t_start, "t_end": t_end,
                    "strand": strand, "matches": matches, "mapq": mapq,
                })
    return alns

def estimate_ref_pos(tg_center, ref_alns):
    """contig 내 트랜스진 중심 좌표를 reference 좌표로 선형 보간"""
    for a in sorted(ref_alns, key=lambda x: x["q_start"]):
        if a["q_start"] <= tg_center <= a["q_end"]:
            ratio = (tg_center - a["q_start"]) / max(1, a["q_end"] - a["q_start"])
            ref_pos = int(a["t_start"] + ratio * (a["t_end"] - a["t_start"]))
            return a["chrom"], ref_pos, a["mapq"]
    return "unknown", 0, 0

lines = []
lines.append("=" * 70)
lines.append(f"Transgene Integration Site Report")
lines.append(f"Sample:     {sample_id}")
lines.append(f"Transgene:  {transgene_name}")
lines.append("=" * 70)

integration_chr = "unknown"
integration_pos = "0"
hap1_count = 0
hap2_count = 0
all_positions = []

tsv_rows = []
tsv_rows.append("\t".join([
    "sample_id", "transgene", "haplotype", "contig", "contig_len",
    "copy_num", "copy_idx", "contig_pos_start", "contig_pos_end",
    "tg_coverage_pct", "strand", "ref_chrom", "ref_pos_est", "ref_mapq",
    "inter_copy_gap_bp",
]))

for hap, tg_paf, ref_paf in [
    ("hap1", "~{hap1_tg_paf}", "~{hap1_ref_paf}"),
    ("hap2", "~{hap2_tg_paf}", "~{hap2_ref_paf}"),
]:
    tg_hits = parse_tg_paf(tg_paf)
    ref_alns = parse_ref_paf(ref_paf)

    lines.append(f"\n{'─'*60}")
    lines.append(f"[ {hap.upper()} ]")

    if not tg_hits:
        lines.append("  트랜스진 hit 없음")
        continue

    count = len(tg_hits)
    if hap == "hap1":
        hap1_count = count
    else:
        hap2_count = count

    for contig, copies in tg_hits.items():
        copies_sorted = sorted(copies, key=lambda x: x["c_start"])
        copy_num = len(copies_sorted)
        lines.append(f"\n  Contig: {contig} ({copies_sorted[0]['contig_len']:,} bp)")
        lines.append(f"  트랜스진 카피 수: {copy_num}개")

        contig_ref_alns = ref_alns.get(contig, [])
        if contig_ref_alns:
            chroms = list(set(a["chrom"] for a in contig_ref_alns))
            lines.append(f"  Reference mapping: {', '.join(chroms)} ({len(contig_ref_alns)} block(s))")

        for i, cp in enumerate(copies_sorted, 1):
            tg_cov = (cp["tg_end"] - cp["tg_start"]) / cp["tg_len"] * 100
            tg_center = (cp["c_start"] + cp["c_end"]) // 2
            chrom, ref_pos, ref_mapq = estimate_ref_pos(tg_center, contig_ref_alns)

            lines.append(f"\n    Copy {i}:")
            lines.append(f"      Contig pos:   {cp['c_start']:,} - {cp['c_end']:,}")
            lines.append(f"      TG coverage:  {cp['tg_start']}-{cp['tg_end']}/{cp['tg_len']} ({tg_cov:.1f}%)")
            lines.append(f"      Direction:    {cp['strand']} ({'역방향' if cp['strand']=='-' else '정방향'})")
            lines.append(f"      Ref position: {chrom}:{ref_pos:,} (mapq={ref_mapq})")

            if chrom != "unknown" and ref_pos > 0:
                all_positions.append((chrom, ref_pos))
                if integration_chr == "unknown":
                    integration_chr = chrom
                    integration_pos = str(ref_pos)

            gap_bp = ""
            if i < copy_num:
                gap = copies_sorted[i]["c_start"] - cp["c_end"]
                gap_bp = str(gap)
                lines.append(f"      → 다음 카피까지 간격: {gap:,} bp")

            tsv_rows.append("\t".join([
                sample_id, transgene_name, hap, contig,
                str(copies_sorted[0]["contig_len"]),
                str(copy_num), str(i),
                str(cp["c_start"]), str(cp["c_end"]),
                f"{tg_cov:.1f}",
                cp["strand"],
                chrom, str(ref_pos), str(ref_mapq),
                gap_bp,
            ]))

lines.append(f"\n{'='*70}")
lines.append("통합 위치 요약")
lines.append(f"{'='*70}")
lines.append(f"  Chromosome:  {integration_chr}")
lines.append(f"  Position:    ~{int(integration_pos):,} (보간 추정값, ±수kb 오차)")
lines.append(f"  hap1 chimeric contigs: {hap1_count}")
lines.append(f"  hap2 chimeric contigs: {hap2_count}")
lines.append("")
lines.append("  ※ 정확한 breakpoint는 CIGAR 파싱 또는 IGV 시각화로 확인 권장")
lines.append(f"{'='*70}")

report_text = "\n".join(lines)
print(report_text)

with open(f"{out_prefix}.integration_report.txt", "w") as f:
    f.write(report_text + "\n")

with open(f"{out_prefix}.integration_sites.tsv", "w") as f:
    f.write("\n".join(tsv_rows) + "\n")

with open("stat_hap1_hit_count.txt", "w") as f:
    f.write(str(hap1_count))
with open("stat_hap2_hit_count.txt", "w") as f:
    f.write(str(hap2_count))
with open("stat_integration_chr.txt", "w") as f:
    f.write(integration_chr)
with open("stat_integration_pos.txt", "w") as f:
    f.write(integration_pos)
PYEOF
  >>>

  output {
    File   report_txt           = "~{out_prefix}.integration_report.txt"
    File   integration_tsv      = "~{out_prefix}.integration_sites.tsv"
    String stat_hap1_hit_count  = read_string("stat_hap1_hit_count.txt")
    String stat_hap2_hit_count  = read_string("stat_hap2_hit_count.txt")
    String stat_integration_chr = read_string("stat_integration_chr.txt")
    String stat_integration_pos = read_string("stat_integration_pos.txt")
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
