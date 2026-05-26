version 1.0

import "../structs.wdl"

# ── Step 3-1 ──────────────────────────────────────────────────────────────────

task map_transgene_to_assembly {
  meta { description: "Map transgene FASTA to assembled contigs (chimeric contig 탐색)" }

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
    minimap2 -cx asm20 --cs -t ~{threads} \
      ~{assembly_fasta} ~{transgene_fasta} \
      > ~{out_prefix}.tg_vs_asm.paf
  >>>

  output { File paf = "~{out_prefix}.tg_vs_asm.paf" }

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

# ── Step 3-2 ──────────────────────────────────────────────────────────────────

task extract_chimeric_contigs {
  meta { description: "Transgene 서열을 포함하는 chimeric contig 추출" }

  parameter_meta {
    tg_paf:               { name: "PAF from map_transgene_to_assembly" }
    assembly_fasta:       { name: "Haplotype assembly FASTA" }
    out_prefix:           { name: "Output prefix" }
    min_match_bp:         { name: "Minimum matching bases (default: 500)" }
    runtime_attributes:   { name: "Runtime attribute structure" }
    chimeric_fasta:       { name: "FASTA of chimeric contigs" }
    chimeric_contig_list: { name: "List of chimeric contig names" }
    hit_count:            { name: "Number of chimeric contigs found" }
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
    # mapq 필터 없이 match 길이만으로 판단
    # (다중 카피 → mapq 낮아지는 정상 현상)
    awk '$10 >= ~{min_match_bp} {print $6}' ~{tg_paf} \
      | sort -u > ~{out_prefix}.chimeric_contigs.txt

    COUNT=$(wc -l < ~{out_prefix}.chimeric_contigs.txt)
    echo "Chimeric contigs found: ${COUNT}"
    echo "${COUNT}" > hit_count.txt

    if [ "${COUNT}" -eq 0 ]; then
      echo "WARNING: No chimeric contigs found." >&2
      touch ~{out_prefix}.chimeric_contigs.fa
    else
      samtools faidx ~{assembly_fasta}
      xargs samtools faidx ~{assembly_fasta} \
        < ~{out_prefix}.chimeric_contigs.txt \
        > ~{out_prefix}.chimeric_contigs.fa
    fi

    echo "=== Transgene hit summary ==="
    awk '$10 >= ~{min_match_bp} {
      tg_cov = ($4-$3)/$2*100
      printf "  contig=%-20s len=%d tg_pos=%d-%d/%d (%.1f%%) strand=%s matches=%d\n",
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

# ── Step 4-1 ──────────────────────────────────────────────────────────────────

task align_chimeric_to_ref {
  meta { description: "Chimeric contig → reference 정렬 (split alignment로 통합 위치 파악)" }

  parameter_meta {
    chimeric_fasta:     { name: "FASTA of chimeric contigs" }
    ref_fasta:          { name: "Reference genome FASTA" }
    out_prefix:         { name: "Output prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    paf:                { name: "PAF: chimeric contigs vs reference" }
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
      minimap2 -cx asm5 --cs --eqx -t ~{threads} \
        ~{ref_fasta} ~{chimeric_fasta} \
        > ~{out_prefix}.chimeric_to_ref.paf
    fi
  >>>

  output { File paf = "~{out_prefix}.chimeric_to_ref.paf" }

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

# ── Step 4-2 (경량) ───────────────────────────────────────────────────────────

task extract_integration_coords {
  meta {
    description: "Chimeric-to-ref PAF에서 대략적 삽입 위치(chr/pos)를 빠르게 추출 — 이후 hybrid reference 분석의 시작점으로 사용"
  }

  parameter_meta {
    hap1_ref_paf:       { name: "PAF: hap1 chimeric contigs vs reference" }
    hap2_ref_paf:       { name: "PAF: hap2 chimeric contigs vs reference" }
    min_match_bp:       { name: "Minimum matching bases to consider (default: 1000)" }
    runtime_attributes: { name: "Runtime attribute structure" }
    integration_chr:    { name: "Estimated integration chromosome" }
    integration_pos:    { name: "Estimated integration position (bp)" }
  }

  input {
    File   hap1_ref_paf
    File   hap2_ref_paf
    Int    min_match_bp = 1000
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail
    echo "unknown" > integration_chr.txt
    echo "0"       > integration_pos.txt

    for PAF in "~{hap1_ref_paf}" "~{hap2_ref_paf}"; do
      if [ -s "$PAF" ]; then
        # PAF 블록을 qname, qstart 순으로 정렬한 뒤
        # 인접 블록 사이의 query gap이 TG insert 크기(500~50000 bp)인 지점을 찾는다.
        # 그 gap 직전 블록의 reference 좌표가 실제 삽입 위치.
        #   minus strand 블록: qend 위치는 ref tstart에 해당
        #   plus  strand 블록: qend 위치는 ref tend에 해당
        sort -k1,1 -k3,3n "$PAF" | awk -v min_gap=~{min_match_bp} '
          prev_qname != "" && $1 == prev_qname {
            gap = $3 - prev_qend
            if (gap >= min_gap && gap <= 50000) {
              if (prev_strand == "-") {
                chr = prev_tname; pos = prev_tstart
              } else {
                chr = prev_tname; pos = prev_tend
              }
              print chr > "integration_chr.txt"
              print pos > "integration_pos.txt"
              exit
            }
          }
          {
            prev_qname=$1; prev_qend=$4
            prev_tname=$6; prev_tstart=$8; prev_tend=$9; prev_strand=$5
          }
        ' -
        CHR=$(cat integration_chr.txt)
        if [ "$CHR" != "unknown" ]; then break; fi
      fi
    done

    echo "Integration coords: $(cat integration_chr.txt):$(cat integration_pos.txt)"
  >>>

  output {
    String integration_chr = read_string("integration_chr.txt")
    String integration_pos = read_string("integration_pos.txt")
  }

  runtime {
    docker:                "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
    cpu:                   1
    memory:                "2 GiB"
    disk:                  "5 GB"
    disks:                 "local-disk 5 HDD"
    preemptible:           runtime_attributes.preemptible_tries
    maxRetries:            runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones:                 runtime_attributes.zones
    cpuPlatform:           runtime_attributes.cpuPlatform
  }
}

# ── Step 5-1 ──────────────────────────────────────────────────────────────────

task build_hybrid_ref {
  meta {
    description: "Reference에서 삽입 염색체 서열 추출 후 transgene을 붙여 hybrid reference 생성"
  }

  parameter_meta {
    ref_fasta:          { name: "Reference genome FASTA (with .fai alongside)" }
    transgene_fasta:    { name: "Transgene sequence FASTA" }
    integration_chr:    { name: "Chromosome to extract (e.g. chr2)" }
    out_prefix:         { name: "Output prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    hybrid_ref_fa:      { name: "Hybrid FASTA: integration_chr + transgene" }
    hybrid_ref_fai:     { name: "FAI index of hybrid FASTA" }
    tg_seqname:         { name: "Transgene sequence name in the hybrid reference" }
  }

  input {
    File   ref_fasta
    File   transgene_fasta
    String integration_chr
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(ref_fasta, "GB") / 10 + 5)

  command <<<
    set -euo pipefail

    # FAI가 없으면 생성 (ref_fasta 경로 기준 — 이미 존재하는 게 정상)
    if [ ! -f "~{ref_fasta}.fai" ]; then
      echo "FAI not found, creating..." >&2
      samtools faidx "~{ref_fasta}"
    fi

    # 염색체 추출 + transgene 추가
    samtools faidx "~{ref_fasta}" "~{integration_chr}" > "~{out_prefix}.hybrid_ref.fa"
    awk '/^>/{gsub(/ .*/,""); print; next} {print}' "~{transgene_fasta}" >> "~{out_prefix}.hybrid_ref.fa"
    samtools faidx "~{out_prefix}.hybrid_ref.fa"

    # transgene sequence name 저장
    awk '/^>/{gsub(/^>/,""); print $1; exit}' "~{transgene_fasta}" > tg_seqname.txt

    echo "Hybrid ref sequences:"
    awk '{printf "  %s  %d bp\n", $1, $2}' "~{out_prefix}.hybrid_ref.fa.fai"
  >>>

  output {
    File   hybrid_ref_fa  = "~{out_prefix}.hybrid_ref.fa"
    File   hybrid_ref_fai = "~{out_prefix}.hybrid_ref.fa.fai"
    String tg_seqname     = read_string("tg_seqname.txt")
  }

  runtime {
    docker:                "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
    cpu:                   2
    memory:                "8 GiB"
    disk:                  disk_size + " GB"
    disks:                 "local-disk " + disk_size + " HDD"
    preemptible:           runtime_attributes.preemptible_tries
    maxRetries:            runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones:                 runtime_attributes.zones
    cpuPlatform:           runtime_attributes.cpuPlatform
  }
}

# ── Step 5-2 ──────────────────────────────────────────────────────────────────

task extract_region_reads {
  meta {
    description: "삽입 위치 주변 HiFi reads를 haplotagged BAM에서 추출"
  }

  parameter_meta {
    haplotagged_bam:    { name: "HiFi reads aligned to reference (haplotagged BAM)" }
    haplotagged_bai:    { name: "BAM index (.bai)" }
    integration_chr:    { name: "Chromosome of integration site" }
    integration_pos:    { name: "Approximate integration position (bp)" }
    window_bp:          { name: "Extraction window on each side of integration_pos (default: 300000)" }
    out_prefix:         { name: "Output prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    region_reads_fq:    { name: "Extracted reads in FASTQ.gz" }
    read_count:         { name: "Number of extracted reads" }
  }

  input {
    File   haplotagged_bam
    File   haplotagged_bai
    String integration_chr
    String integration_pos
    Int    window_bp = 300000
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(haplotagged_bam, "GB") / 50 + 5)

  command <<<
    set -euo pipefail

    # BAI를 BAM 이름 규칙에 맞게 symlink (samtools view 요구)
    BAM="~{haplotagged_bam}"
    BAI="~{haplotagged_bai}"
    ln -sf "$BAI" "${BAM}.bai" 2>/dev/null || true

    CENTER=~{integration_pos}
    WINDOW=~{window_bp}
    START=$(( CENTER > WINDOW ? CENTER - WINDOW : 0 ))
    END=$(( CENTER + WINDOW ))
    REGION="~{integration_chr}:${START}-${END}"
    echo "Extracting reads from: $REGION"

    # -F 2304: supplementary(2048) + secondary(256) 제외 → 중복 없이 primary reads만
    samtools view -b -F 2304 "$BAM" "$REGION" \
      | samtools fastq - \
      | gzip > ~{out_prefix}.region_reads.fastq.gz

    COUNT=$(zcat ~{out_prefix}.region_reads.fastq.gz | awk 'NR%4==1' | wc -l)
    echo "Extracted reads: $COUNT"
    echo "$COUNT" > read_count.txt
  >>>

  output {
    File   region_reads_fq = "~{out_prefix}.region_reads.fastq.gz"
    String read_count      = read_string("read_count.txt")
  }

  runtime {
    docker:                "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
    cpu:                   4
    memory:                "8 GiB"
    disk:                  disk_size + " GB"
    disks:                 "local-disk " + disk_size + " HDD"
    preemptible:           runtime_attributes.preemptible_tries
    maxRetries:            runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones:                 runtime_attributes.zones
    cpuPlatform:           runtime_attributes.cpuPlatform
  }
}

# ── Step 5-3 ──────────────────────────────────────────────────────────────────

task align_to_hybrid_ref {
  meta {
    description: "추출된 HiFi reads를 hybrid reference에 재정렬 (PAF 출력)"
  }

  parameter_meta {
    hybrid_ref_fa:      { name: "Hybrid reference FASTA" }
    hybrid_ref_fai:     { name: "Hybrid reference FAI" }
    region_reads_fq:    { name: "Region reads FASTQ.gz" }
    out_prefix:         { name: "Output prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    hybrid_paf:         { name: "PAF: region reads vs hybrid reference" }
  }

  input {
    File   hybrid_ref_fa
    File   hybrid_ref_fai
    File   region_reads_fq
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 16
  Int mem_gb    = 24
  Int disk_size = ceil(size(region_reads_fq, "GB") * 5 + size(hybrid_ref_fa, "GB") * 2 + 5)

  command <<<
    set -euo pipefail
    minimap2 --version
    minimap2 \
      -cx map-hifi \
      --cs \
      --secondary=no \
      -t ~{threads} \
      ~{hybrid_ref_fa} \
      ~{region_reads_fq} \
      > ~{out_prefix}.hybrid_aligned.paf

    echo "Total alignments: $(wc -l < ~{out_prefix}.hybrid_aligned.paf)"
  >>>

  output { File hybrid_paf = "~{out_prefix}.hybrid_aligned.paf" }

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

# ── Step 5-4 ──────────────────────────────────────────────────────────────────

task detect_chimeric_reads {
  meta {
    description: "Hybrid reference PAF에서 chimeric reads 검출 및 정확한 breakpoint 결정"
  }

  parameter_meta {
    hybrid_paf:         { name: "PAF: region reads vs hybrid reference" }
    tg_seqname:         { name: "Transgene sequence name in hybrid reference" }
    out_prefix:         { name: "Output prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    breakpoint_chr:     { name: "Breakpoint chromosome" }
    breakpoint_pos:     { name: "Breakpoint position (median of chimeric reads)" }
    breakpoint_range:   { name: "Breakpoint position range (min-max)" }
    chimeric_count:     { name: "Number of chimeric reads supporting the breakpoint" }
    chimeric_tsv:       { name: "Per-read chimeric alignment details (TSV)" }
  }

  input {
    File   hybrid_paf
    String tg_seqname
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail

    python3 << 'PYEOF'
import re, os
from collections import defaultdict

hybrid_paf = "~{hybrid_paf}"
tg_seqname = "~{tg_seqname}"
out_prefix = "~{out_prefix}"

# ── PAF 파싱: query별 정렬 블록 수집 ─────────────────────────────────────────
reads = defaultdict(list)
with open(hybrid_paf) as f:
    for line in f:
        p = line.strip().split('\t')
        if len(p) < 12:
            continue
        reads[p[0]].append({
            "q_start": int(p[2]), "q_end": int(p[3]),
            "strand":  p[4],
            "target":  p[5],
            "t_start": int(p[7]), "t_end": int(p[8]),
            "matches": int(p[9]), "mapq": int(p[11]),
        })

# ── Chimeric reads 탐색 ────────────────────────────────────────────────────────
# chimeric = 동일 read가 host chr AND transgene 모두에 정렬된 경우
chimeric = []
for qname, alns in reads.items():
    targets  = set(a["target"] for a in alns)
    host_set = {t for t in targets if t != tg_seqname}
    has_host = bool(host_set)
    has_tg   = tg_seqname in targets
    if has_host and has_tg:
        chimeric.append({
            "name":      qname,
            "host_alns": [a for a in alns if a["target"] != tg_seqname],
            "tg_alns":   [a for a in alns if a["target"] == tg_seqname],
        })

print(f"Total chimeric reads: {len(chimeric)}")

# ── Breakpoint 결정 ───────────────────────────────────────────────────────────
# Read 내 좌표 기준으로 host와 TG가 인접한 경우의 host reference 좌표를 breakpoint로 사용
# tolerance: 인접 판정 허용 거리 (bp)
TOLERANCE = 1000

breakpoints = []
host_chrs   = []
tsv_rows    = ["read_name\thost_chr\thost_q_start\thost_q_end\thost_t_start\thost_t_end\ttg_q_start\ttg_q_end\tbreakpoint_pos\tjunction_type"]

for c in chimeric:
    for ha in c["host_alns"]:
        for ta in c["tg_alns"]:
            # Case 1: host 정렬 끝 ≈ TG 정렬 시작 (read 상에서 host→TG 방향)
            if abs(ha["q_end"] - ta["q_start"]) <= TOLERANCE:
                bp = ha["t_end"]
                jtype = "host->TG"
                breakpoints.append(bp)
                host_chrs.append(ha["target"])
                tsv_rows.append(f'{c["name"]}\t{ha["target"]}\t{ha["q_start"]}\t{ha["q_end"]}\t{ha["t_start"]}\t{ha["t_end"]}\t{ta["q_start"]}\t{ta["q_end"]}\t{bp}\t{jtype}')
            # Case 2: TG 정렬 끝 ≈ host 정렬 시작 (read 상에서 TG→host 방향)
            elif abs(ta["q_end"] - ha["q_start"]) <= TOLERANCE:
                bp = ha["t_start"]
                jtype = "TG->host"
                breakpoints.append(bp)
                host_chrs.append(ha["target"])
                tsv_rows.append(f'{c["name"]}\t{ha["target"]}\t{ha["q_start"]}\t{ha["q_end"]}\t{ha["t_start"]}\t{ha["t_end"]}\t{ta["q_start"]}\t{ta["q_end"]}\t{bp}\t{jtype}')

bp_chr = max(set(host_chrs), key=host_chrs.count) if host_chrs else "unknown"

if breakpoints:
    bp_sorted  = sorted(breakpoints)
    bp_median  = bp_sorted[len(bp_sorted) // 2]
    bp_min     = min(breakpoints)
    bp_max     = max(breakpoints)
    bp_range   = f"{bp_chr}:{bp_min:,}-{bp_max:,}"
    print(f"Breakpoint chromosome: {bp_chr}")
    print(f"Breakpoint (median):   {bp_chr}:{bp_median:,}")
    print(f"Breakpoint range:      {bp_range}")
    print(f"Supporting reads:      {len(breakpoints)}")
else:
    bp_median = 0
    bp_range  = "N/A"
    print("WARNING: No clear breakpoint detected from chimeric reads.")

with open(f"{out_prefix}.chimeric_reads.tsv", "w") as f:
    f.write("\n".join(tsv_rows) + "\n")
with open("breakpoint_chr.txt",   "w") as f: f.write(bp_chr)
with open("breakpoint_pos.txt",   "w") as f: f.write(str(bp_median))
with open("breakpoint_range.txt", "w") as f: f.write(bp_range)
with open("chimeric_count.txt",   "w") as f: f.write(str(len(chimeric)))
PYEOF
  >>>

  output {
    File   chimeric_tsv     = "~{out_prefix}.chimeric_reads.tsv"
    String breakpoint_chr   = read_string("breakpoint_chr.txt")
    String breakpoint_pos   = read_string("breakpoint_pos.txt")
    String breakpoint_range = read_string("breakpoint_range.txt")
    String chimeric_count   = read_string("chimeric_count.txt")
  }

  runtime {
    docker:                "quay.io/biocontainers/busco:5.7.1--pyhdfd78af_1"
    cpu:                   2
    memory:                "8 GiB"
    disk:                  "10 GB"
    disks:                 "local-disk 10 HDD"
    preemptible:           runtime_attributes.preemptible_tries
    maxRetries:            runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones:                 runtime_attributes.zones
    cpuPlatform:           runtime_attributes.cpuPlatform
  }
}

# ── Step 5-5 (Optional) ───────────────────────────────────────────────────────

task query_gene_annotation {
  meta {
    description: "Breakpoint 주변 유전자 어노테이션 조회 (Gencode GTF tabix)"
  }

  parameter_meta {
    annotation_gtf_bgz:    { name: "Sorted, bgzipped Gencode GTF" }
    annotation_gtf_tbi:    { name: "Tabix index for the GTF" }
    breakpoint_chr:        { name: "Breakpoint chromosome" }
    breakpoint_pos:        { name: "Breakpoint position (bp)" }
    window_bp:             { name: "Query window on each side (default: 200000)" }
    out_prefix:            { name: "Output prefix" }
    runtime_attributes:    { name: "Runtime attribute structure" }
    nearest_gene_name:     { name: "Gene containing or nearest to the breakpoint" }
    nearest_gene_type:     { name: "Biotype of the nearest gene" }
    insertion_feature_type:{ name: "Feature type at breakpoint: exon/intron/intergenic" }
    nearest_exon_dist_bp:  { name: "Distance to nearest exon boundary (bp)" }
    annotation_tsv:        { name: "Genes near breakpoint (TSV)" }
  }

  input {
    File   annotation_gtf_bgz
    File   annotation_gtf_tbi
    String breakpoint_chr
    String breakpoint_pos
    Int    window_bp = 200000
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail

    python3 << 'PYEOF'
import subprocess, re, os

gtf_bgz    = "~{annotation_gtf_bgz}"
gtf_tbi    = "~{annotation_gtf_tbi}"
bp_chr     = "~{breakpoint_chr}"
bp_pos     = int("~{breakpoint_pos}") if "~{breakpoint_pos}" != "0" else 0
window     = ~{window_bp}
out_prefix = "~{out_prefix}"

def parse_attr(attr_str, key):
    m = re.search(rf'{key} "([^"]+)"', attr_str)
    return m.group(1) if m else ""

if bp_pos == 0:
    print("Breakpoint position is 0, skipping annotation query.")
    for fn, val in [
        ("nearest_gene_name.txt",      "N/A"),
        ("nearest_gene_type.txt",      "N/A"),
        ("insertion_feature_type.txt", "N/A"),
        ("nearest_exon_dist_bp.txt",   "N/A"),
    ]:
        open(fn, "w").write(val)
    open(f"{out_prefix}.annotation.tsv", "w").write(
        "gene_name\tgene_type\tgene_start\tgene_end\tstrand\tdistance_to_bp\n")
    exit(0)

query_start = max(0, bp_pos - window)
query_end   = bp_pos + window
region      = f"{bp_chr}:{query_start}-{query_end}"
print(f"Querying: {region}")

# GTF tabix를 실행할 때 tbi도 같은 이름으로 접근 가능한지 확인
tbi_expected = gtf_bgz + ".tbi"
if not os.path.exists(tbi_expected):
    import shutil
    shutil.copy(gtf_tbi, tbi_expected)

result = subprocess.run(
    ["tabix", gtf_bgz, region],
    capture_output=True, text=True
)
lines = result.stdout.strip().split("\n") if result.stdout.strip() else []

genes, exons = [], []
for line in lines:
    if not line or line.startswith("#"):
        continue
    p = line.split("\t")
    if len(p) < 9:
        continue
    feature = p[2]
    start, end = int(p[3]), int(p[4])
    strand = p[6]
    attrs  = p[8]
    gname  = parse_attr(attrs, "gene_name")
    gtype  = parse_attr(attrs, "gene_type")
    if feature == "gene":
        genes.append({"name": gname, "type": gtype, "start": start, "end": end, "strand": strand})
    elif feature == "exon":
        exons.append({"gene": gname, "start": start, "end": end})

# breakpoint를 포함하는 유전자 찾기
containing_genes = [g for g in genes if g["start"] <= bp_pos <= g["end"]]
if not containing_genes:
    # 가장 가까운 유전자
    containing_genes = sorted(genes, key=lambda g: min(abs(g["start"]-bp_pos), abs(g["end"]-bp_pos)))

nearest_gene = containing_genes[0] if containing_genes else None

# insertion feature type: exon / intron / intergenic
insertion_feature = "intergenic"
if nearest_gene:
    if nearest_gene["start"] <= bp_pos <= nearest_gene["end"]:
        gene_exons = [e for e in exons if e["gene"] == nearest_gene["name"]]
        in_exon = any(e["start"] <= bp_pos <= e["end"] for e in gene_exons)
        insertion_feature = "exon" if in_exon else "intron"
    else:
        insertion_feature = "intergenic"

# 가장 가까운 exon 경계까지 거리
gene_exons = [e for e in exons if nearest_gene and e["gene"] == nearest_gene["name"]]
exon_dists = []
for e in gene_exons:
    exon_dists.extend([abs(e["start"] - bp_pos), abs(e["end"] - bp_pos)])
nearest_exon_dist = min(exon_dists) if exon_dists else -1

gene_name = nearest_gene["name"] if nearest_gene else "intergenic"
gene_type = nearest_gene["type"] if nearest_gene else "N/A"

print(f"Nearest gene:      {gene_name} ({gene_type})")
print(f"Insertion feature: {insertion_feature}")
print(f"Nearest exon dist: {nearest_exon_dist:,} bp")

# TSV 출력
tsv_rows = ["gene_name\tgene_type\tgene_start\tgene_end\tstrand\tdistance_to_bp"]
for g in sorted(containing_genes[:10], key=lambda x: min(abs(x["start"]-bp_pos), abs(x["end"]-bp_pos))):
    dist = 0 if g["start"] <= bp_pos <= g["end"] else min(abs(g["start"]-bp_pos), abs(g["end"]-bp_pos))
    tsv_rows.append(f'{g["name"]}\t{g["type"]}\t{g["start"]}\t{g["end"]}\t{g["strand"]}\t{dist}')
with open(f"{out_prefix}.annotation.tsv", "w") as f:
    f.write("\n".join(tsv_rows) + "\n")

for fn, val in [
    ("nearest_gene_name.txt",      gene_name),
    ("nearest_gene_type.txt",      gene_type),
    ("insertion_feature_type.txt", insertion_feature),
    ("nearest_exon_dist_bp.txt",   str(nearest_exon_dist)),
]:
    with open(fn, "w") as f:
        f.write(val)
PYEOF
  >>>

  output {
    File   annotation_tsv         = "~{out_prefix}.annotation.tsv"
    String nearest_gene_name      = read_string("nearest_gene_name.txt")
    String nearest_gene_type      = read_string("nearest_gene_type.txt")
    String insertion_feature_type = read_string("insertion_feature_type.txt")
    String nearest_exon_dist_bp   = read_string("nearest_exon_dist_bp.txt")
  }

  runtime {
    docker:                "quay.io/biocontainers/busco:5.7.1--pyhdfd78af_1"
    cpu:                   2
    memory:                "4 GiB"
    disk:                  "10 GB"
    disks:                 "local-disk 10 HDD"
    preemptible:           runtime_attributes.preemptible_tries
    maxRetries:            runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones:                 runtime_attributes.zones
    cpuPlatform:           runtime_attributes.cpuPlatform
  }
}

# ── Step 6: 통합 리포트 ───────────────────────────────────────────────────────

task integration_report {
  meta {
    description: "De novo + hybrid reference 분석 결과를 통합한 리포트 생성 (txt / TSV / HTML)"
  }

  parameter_meta {
    sample_id:              { name: "Sample ID" }
    transgene_name:         { name: "Transgene construct name" }
    hap1_tg_paf:            { name: "PAF: transgene vs hap1 contigs" }
    hap2_tg_paf:            { name: "PAF: transgene vs hap2 contigs" }
    hap1_ref_paf:           { name: "PAF: hap1 chimeric contigs vs reference" }
    hap2_ref_paf:           { name: "PAF: hap2 chimeric contigs vs reference" }
    out_prefix:             { name: "Output prefix" }
    precise_breakpoint_chr: { name: "Precise breakpoint chromosome (from hybrid ref; 'N/A' if not run)" }
    precise_breakpoint_pos: { name: "Precise breakpoint position (from hybrid ref; '0' if not run)" }
    breakpoint_range:       { name: "Breakpoint position range string (from hybrid ref)" }
    chimeric_count:         { name: "Number of chimeric reads supporting the breakpoint" }
    nearest_gene_name:      { name: "Gene at/near the breakpoint (from annotation; 'N/A' if not run)" }
    nearest_gene_type:      { name: "Biotype of nearest gene" }
    insertion_feature_type: { name: "exon / intron / intergenic" }
    nearest_exon_dist_bp:   { name: "Distance to nearest exon boundary (bp)" }
    runtime_attributes:     { name: "Runtime attribute structure" }
    report_txt:             { name: "Integration site report (text)" }
    report_html:            { name: "Integration site report (HTML)" }
    integration_tsv:        { name: "Integration sites table (TSV)" }
    stat_hap1_hit_count:    { name: "Hap1 chimeric contig count" }
    stat_hap2_hit_count:    { name: "Hap2 chimeric contig count" }
    stat_integration_chr:   { name: "Integration chromosome (de novo estimate)" }
    stat_integration_pos:   { name: "Integration position (de novo estimate)" }
  }

  input {
    String sample_id
    String transgene_name
    File   hap1_tg_paf
    File   hap2_tg_paf
    File   hap1_ref_paf
    File   hap2_ref_paf
    String out_prefix

    # hybrid reference 분석 결과 (선택적 — 미실행 시 기본값)
    String precise_breakpoint_chr = "N/A"
    String precise_breakpoint_pos = "0"
    String breakpoint_range       = "N/A"
    String chimeric_count         = "0"

    # 유전자 어노테이션 결과 (선택적 — 미실행 시 기본값)
    String nearest_gene_name      = "N/A"
    String nearest_gene_type      = "N/A"
    String insertion_feature_type = "N/A"
    String nearest_exon_dist_bp   = "N/A"

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = 10

  command <<<
    set -euo pipefail

    python3 << 'PYEOF'
import os, re
from collections import defaultdict
from datetime import date

sample_id      = "~{sample_id}"
transgene_name = "~{transgene_name}"
out_prefix     = "~{out_prefix}"
today          = date.today().isoformat()

# hybrid ref 결과
prec_chr       = "~{precise_breakpoint_chr}"
prec_pos_str   = "~{precise_breakpoint_pos}"
prec_pos       = int(prec_pos_str) if prec_pos_str.isdigit() and prec_pos_str != "0" else 0
bp_range       = "~{breakpoint_range}"
chimeric_n     = "~{chimeric_count}"
has_hybrid     = prec_chr != "N/A" and prec_pos > 0

# 어노테이션 결과
gene_name      = "~{nearest_gene_name}"
gene_type      = "~{nearest_gene_type}"
feat_type      = "~{insertion_feature_type}"
exon_dist_str  = "~{nearest_exon_dist_bp}"
exon_dist      = int(exon_dist_str) if exon_dist_str.lstrip("-").isdigit() else -1
has_annot      = gene_name != "N/A"

# ── PAF 파싱 ───────────────────────────────────────────────────────────────────
def parse_tg_paf(paf_file, min_match=500):
    hits = defaultdict(list)
    if not os.path.exists(paf_file) or os.path.getsize(paf_file) == 0:
        return hits
    with open(paf_file) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) < 12 or int(p[9]) < min_match:
                continue
            hits[p[5]].append({
                "tg_start": int(p[2]), "tg_end": int(p[3]), "tg_len": int(p[1]),
                "contig_len": int(p[6]),
                "c_start": int(p[7]), "c_end": int(p[8]),
                "strand": p[4], "matches": int(p[9]), "mapq": int(p[11]),
            })
    return hits

def parse_ref_paf(paf_file, min_match=1000):
    alns = defaultdict(list)
    if not os.path.exists(paf_file) or os.path.getsize(paf_file) == 0:
        return alns
    with open(paf_file) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) < 12 or int(p[9]) < min_match:
                continue
            alns[p[0]].append({
                "q_len": int(p[1]), "q_start": int(p[2]), "q_end": int(p[3]),
                "chrom": p[5], "t_start": int(p[7]), "t_end": int(p[8]),
                "strand": p[4], "matches": int(p[9]), "mapq": int(p[11]),
            })
    return alns

def estimate_ref_pos(tg_center, ref_alns):
    for a in sorted(ref_alns, key=lambda x: x["q_start"]):
        if a["q_start"] <= tg_center <= a["q_end"]:
            ratio = (tg_center - a["q_start"]) / max(1, a["q_end"] - a["q_start"])
            return a["chrom"], int(a["t_start"] + ratio * (a["t_end"] - a["t_start"])), a["mapq"]
    return "unknown", 0, 0

hap_data        = {}
integration_chr = "unknown"
integration_pos = "0"

for hap, tg_paf, ref_paf in [
    ("hap1", "~{hap1_tg_paf}", "~{hap1_ref_paf}"),
    ("hap2", "~{hap2_tg_paf}", "~{hap2_ref_paf}"),
]:
    tg_hits  = parse_tg_paf(tg_paf)
    ref_alns = parse_ref_paf(ref_paf)
    contigs  = []
    for contig, copies in tg_hits.items():
        copies_s = sorted(copies, key=lambda x: x["c_start"])
        ref_aln  = ref_alns.get(contig, [])
        enriched = []
        for cp in copies_s:
            center = (cp["c_start"] + cp["c_end"]) // 2
            chrom, ref_pos, ref_mapq = estimate_ref_pos(center, ref_aln)
            tg_cov = (cp["tg_end"] - cp["tg_start"]) / cp["tg_len"] * 100
            enriched.append({**cp, "chrom": chrom, "ref_pos": ref_pos,
                              "ref_mapq": ref_mapq, "tg_cov": tg_cov})
        gaps = [enriched[i+1]["c_start"] - enriched[i]["c_end"]
                for i in range(len(enriched)-1)]
        contigs.append({
            "name": contig, "contig_len": copies_s[0]["contig_len"],
            "copies": enriched, "gaps": gaps, "ref_alns": ref_aln,
        })
    hap_data[hap] = contigs

for hap in ["hap1", "hap2"]:
    for ctg in hap_data.get(hap, []):
        for cp in ctg["copies"]:
            if cp["chrom"] != "unknown" and cp["ref_pos"] > 0:
                if integration_chr == "unknown":
                    integration_chr = cp["chrom"]
                    integration_pos = str(cp["ref_pos"])

# ── Text Report ────────────────────────────────────────────────────────────────
lines = []
lines.append("=" * 70)
lines.append("Transgene Integration Site Report")
lines.append(f"Sample:     {sample_id}")
lines.append(f"Transgene:  {transgene_name}")
lines.append(f"Generated:  {today}")
lines.append("=" * 70)

if has_hybrid:
    lines.append(f"\n[PRECISE BREAKPOINT — hybrid reference analysis]")
    lines.append(f"  Chromosome: {prec_chr}")
    lines.append(f"  Position:   {prec_pos:,} bp (±4 bp)")
    lines.append(f"  Range:      {bp_range}")
    lines.append(f"  Support:    {chimeric_n} chimeric reads")
if has_annot:
    lines.append(f"\n[GENE ANNOTATION — Gencode]")
    lines.append(f"  Gene:       {gene_name} ({gene_type})")
    lines.append(f"  Feature:    {feat_type}")
    if exon_dist >= 0:
        lines.append(f"  Nearest exon: {exon_dist:,} bp")

hap1_count = 0
hap2_count = 0
tsv_rows = ["\t".join([
    "sample_id", "transgene", "haplotype", "contig", "contig_len",
    "copy_num", "copy_idx", "contig_pos_start", "contig_pos_end",
    "tg_coverage_pct", "strand", "ref_chrom", "ref_pos_est", "ref_mapq",
    "inter_copy_gap_bp",
])]

for hap in ["hap1", "hap2"]:
    contigs = hap_data.get(hap, [])
    lines.append(f"\n{'─'*60}")
    lines.append(f"[ {hap.upper()} ]")
    if not contigs:
        lines.append("  트랜스진 hit 없음")
        continue
    count = len(contigs)
    if hap == "hap1": hap1_count = count
    else: hap2_count = count
    for c in contigs:
        copies_s = c["copies"]
        copy_num = len(copies_s)
        lines.append(f"\n  Contig: {c['name']} ({c['contig_len']:,} bp)")
        lines.append(f"  트랜스진 카피 수: {copy_num}개")
        if c["ref_alns"]:
            chroms = list(set(a["chrom"] for a in c["ref_alns"]))
            lines.append(f"  Reference mapping: {', '.join(chroms)} ({len(c['ref_alns'])} block(s))")
        for i, cp in enumerate(copies_s, 1):
            gap_bp = ""
            if i < copy_num:
                gap = copies_s[i]["c_start"] - cp["c_end"]
                gap_bp = str(gap)
                lines.append(f"      → 다음 카피까지 간격: {gap:,} bp")
            lines.append(f"\n    Copy {i}:")
            lines.append(f"      Contig pos:   {cp['c_start']:,} - {cp['c_end']:,}")
            lines.append(f"      TG coverage:  {cp['tg_start']}-{cp['tg_end']}/{cp['tg_len']} ({cp['tg_cov']:.1f}%)")
            lines.append(f"      Direction:    {cp['strand']} ({'역방향' if cp['strand']=='-' else '정방향'})")
            lines.append(f"      Ref position: {cp['chrom']}:{cp['ref_pos']:,} (mapq={cp['ref_mapq']})")
            tsv_rows.append("\t".join([
                sample_id, transgene_name, hap, c["name"],
                str(c["contig_len"]), str(copy_num), str(i),
                str(cp["c_start"]), str(cp["c_end"]),
                f"{cp['tg_cov']:.1f}", cp["strand"],
                cp["chrom"], str(cp["ref_pos"]), str(cp["ref_mapq"]),
                gap_bp,
            ]))

lines.append(f"\n{'='*70}")
lines.append("통합 위치 요약 (de novo estimate)")
lines.append(f"  Chromosome: {integration_chr}")
lines.append(f"  Position:   {int(integration_pos):,} (보간 추정, ±수kb)")
lines.append(f"  hap1 chimeric contigs: {hap1_count}")
lines.append(f"  hap2 chimeric contigs: {hap2_count}")
lines.append("  ※ 정확한 위치는 hybrid reference 분석 결과 참조")
lines.append(f"{'='*70}")

report_text = "\n".join(lines)
print(report_text)

with open(f"{out_prefix}.integration_report.txt", "w") as f:
    f.write(report_text + "\n")
with open(f"{out_prefix}.integration_sites.tsv", "w") as f:
    f.write("\n".join(tsv_rows) + "\n")

with open("stat_hap1_hit_count.txt", "w") as f: f.write(str(hap1_count))
with open("stat_hap2_hit_count.txt", "w") as f: f.write(str(hap2_count))
with open("stat_integration_chr.txt", "w") as f: f.write(integration_chr)
with open("stat_integration_pos.txt", "w") as f: f.write(integration_pos)

# ── HTML Report ────────────────────────────────────────────────────────────────
def pct_bar(v, color="#4CAF50"):
    return (f'<div style="background:#e9ecef;border-radius:4px;height:8px;width:110px;display:inline-block;vertical-align:middle">'
            f'<div style="background:{color};width:{min(v,100):.1f}%;height:100%;border-radius:4px"></div></div>'
            f' <span style="font-size:0.83em;color:#555">{v:.1f}%</span>')

def strand_badge(s):
    c,lbl = ("#e74c3c","&#8722; 역방향") if s=="-" else ("#27ae60","+ 정방향")
    return f'<span style="background:{c};color:#fff;padding:1px 7px;border-radius:10px;font-size:0.79em">{lbl}</span>'

def make_contig_svg(contig_len, copies):
    W, H, scale = 640, 52, 640/contig_len
    parts = [
        f'<rect x="0" y="18" width="{W}" height="16" rx="3" fill="#b0bec5"/>',
        f'<text x="2" y="14" font-size="11" fill="#546e7a">host</text>',
        f'<text x="{W-2}" y="14" font-size="11" fill="#546e7a" text-anchor="end">{contig_len:,} bp</text>',
    ]
    for i,cp in enumerate(copies):
        x1=cp["c_start"]*scale; w=max((cp["c_end"]-cp["c_start"])*scale,2)
        parts+=[
            f'<rect x="{x1:.1f}" y="12" width="{w:.1f}" height="28" rx="2" fill="#e74c3c" opacity="0.9"/>',
            f'<text x="{x1+w/2:.1f}" y="29" font-size="10" fill="white" text-anchor="middle" font-weight="bold">TG{i+1}</text>',
            f'<text x="{x1:.1f}" y="46" font-size="9" fill="#555">{cp["c_start"]:,}</text>',
        ]
    return f'<svg width="{W}" height="{H}" style="display:block;margin:5px 0">'+"".join(parts)+'</svg>'

def build_hap_html(hap, contigs):
    color="#1565c0" if hap=="hap1" else "#6a1b9a"
    label=hap.upper()
    if not contigs: return f'<p style="color:#888">No chimeric contigs in {label}.</p>'
    blocks=[]
    for c in contigs:
        copies,gaps,ref_alns=c["copies"],c["gaps"],c["ref_alns"]
        copy_rows=""
        for i,cp in enumerate(copies,1):
            copy_rows+=(
                f'<tr><td style="text-align:center;font-weight:600">{i}</td>'
                f'<td style="font-family:monospace;white-space:nowrap">{cp["c_start"]:,}&#8211;{cp["c_end"]:,}</td>'
                f'<td>{pct_bar(cp["tg_cov"])}</td><td>{strand_badge(cp["strand"])}</td>'
                f'<td style="font-family:monospace">{cp["chrom"]}:{cp["ref_pos"]:,}</td>'
                f'<td style="text-align:center">{cp["ref_mapq"]}</td></tr>'
            )
            if i<=len(gaps):
                copy_rows+=f'<tr><td colspan="6" style="text-align:center;color:#888;font-size:0.82em;padding:1px">&#8595; gap {gaps[i-1]:,} bp</td></tr>'
        ref_rows="".join(
            f'<tr><td style="font-family:monospace;white-space:nowrap">{a["chrom"]}:{a["t_start"]:,}&#8211;{a["t_end"]:,}</td>'
            f'<td style="font-family:monospace;white-space:nowrap">{a["q_start"]:,}&#8211;{a["q_end"]:,}</td>'
            f'<td style="text-align:center">{a["strand"]}</td><td style="text-align:center">{a["mapq"]}</td></tr>'
            for a in ref_alns
        )
        n=len(ref_alns)
        ref_note=('<p style="color:#27ae60;font-size:0.82em;margin-top:5px">&#10003; 단일 블록 &#8212; 단순 삽입</p>'
                  if n<=1 else
                  f'<p style="color:#e67e22;font-size:0.82em;margin-top:5px">&#9888; {n}개 블록 &#8594; 추가 구조변이 가능</p>')
        blocks.append(
            f'<div style="background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:16px;margin-bottom:12px">'
            f'<div style="display:flex;align-items:center;gap:10px;margin-bottom:8px">'
            f'<span style="background:{color};color:#fff;padding:2px 10px;border-radius:12px;font-size:0.82em">{label}</span>'
            f'<code style="font-size:0.93em;font-weight:600">{c["name"]}</code>'
            f'<span style="color:#888;font-size:0.85em">{c["contig_len"]:,} bp &middot; {len(copies)}카피</span></div>'
            f'{make_contig_svg(c["contig_len"],copies)}'
            f'<div style="display:grid;grid-template-columns:1fr 1fr;gap:12px;margin-top:10px">'
            f'<div><div style="font-size:0.76em;color:#666;font-weight:600;margin-bottom:3px">트랜스진 카피 상세</div>'
            f'<table style="width:100%;border-collapse:collapse;font-size:0.84em"><thead><tr style="background:#f8f9fa">'
            f'<th style="padding:3px 6px;border-bottom:2px solid #dee2e6">#</th>'
            f'<th style="padding:3px 6px;text-align:left;border-bottom:2px solid #dee2e6">Contig 위치</th>'
            f'<th style="padding:3px 6px;text-align:left;border-bottom:2px solid #dee2e6">TG 커버리지</th>'
            f'<th style="padding:3px 6px;text-align:left;border-bottom:2px solid #dee2e6">방향</th>'
            f'<th style="padding:3px 6px;text-align:left;border-bottom:2px solid #dee2e6">Ref 위치 (추정)</th>'
            f'<th style="padding:3px 6px;border-bottom:2px solid #dee2e6">mapq</th>'
            f'</tr></thead><tbody>{copy_rows}</tbody></table></div>'
            f'<div><div style="font-size:0.76em;color:#666;font-weight:600;margin-bottom:3px">Reference 정렬 블록</div>'
            f'<table style="width:100%;border-collapse:collapse;font-size:0.84em"><thead><tr style="background:#f8f9fa">'
            f'<th style="padding:3px 6px;text-align:left;border-bottom:2px solid #dee2e6">Ref 좌표</th>'
            f'<th style="padding:3px 6px;text-align:left;border-bottom:2px solid #dee2e6">Contig 좌표</th>'
            f'<th style="padding:3px 6px;border-bottom:2px solid #dee2e6">Strand</th>'
            f'<th style="padding:3px 6px;border-bottom:2px solid #dee2e6">mapq</th>'
            f'</tr></thead><tbody>{ref_rows}</tbody></table>{ref_note}</div>'
            f'</div></div>'
        )
    return "\n".join(blocks)

hap1_html=build_hap_html("hap1",hap_data.get("hap1",[]))
hap2_html=build_hap_html("hap2",hap_data.get("hap2",[]))

def hap_val(hap,fn,default="N/A"):
    d=hap_data.get(hap,[])
    return fn(d[0]) if d else default

def ref_badge(hap):
    d=hap_data.get(hap,[])
    if not d: return "N/A"
    n=len(d[0]["ref_alns"])
    if n<=1: return f'<span style="background:#d4edda;color:#155724;padding:1px 6px;border-radius:8px;font-size:0.79em">{n} 단순</span>'
    return f'<span style="background:#fff3cd;color:#856404;padding:1px 6px;border-radius:8px;font-size:0.79em">{n} 복잡</span>'

cmp_rows="".join(
    f'<tr><td style="color:#555;padding:4px 7px">{r}</td>'
    f'<td style="text-align:center;padding:4px 7px">{v1}</td>'
    f'<td style="text-align:center;padding:4px 7px">{v2}</td></tr>'
    for r,v1,v2 in [
        ("Contig",     hap_val("hap1",lambda c:f'<code>{c["name"]}</code>'),hap_val("hap2",lambda c:f'<code>{c["name"]}</code>')),
        ("Contig 크기",hap_val("hap1",lambda c:f'{c["contig_len"]:,} bp'),hap_val("hap2",lambda c:f'{c["contig_len"]:,} bp')),
        ("카피 수",    hap_val("hap1",lambda c:str(len(c["copies"])),"0"),hap_val("hap2",lambda c:str(len(c["copies"])),"0")),
        ("카피 간격",  hap_val("hap1",lambda c:f'{c["gaps"][0]:,} bp' if c["gaps"] else "&#8212;"),
                       hap_val("hap2",lambda c:f'{c["gaps"][0]:,} bp' if c["gaps"] else "&#8212;")),
        ("방향",       hap_val("hap1",lambda c:"역방향 (&#8722;)" if c["copies"] and c["copies"][0]["strand"]=="-" else "정방향 (+)"),
                       hap_val("hap2",lambda c:"역방향 (&#8722;)" if c["copies"] and c["copies"][0]["strand"]=="-" else "정방향 (+)")),
        ("Ref 정렬 블록",ref_badge("hap1"),ref_badge("hap2")),
    ]
)

# 요약 카드 데이터
all_copies  = [cp for hap in ["hap1","hap2"] for c in hap_data.get(hap,[]) for cp in c["copies"]]
max_copies  = max((len(c["copies"]) for hap in ["hap1","hap2"] for c in hap_data.get(hap,[])), default=0)
all_strands = [cp["strand"] for cp in all_copies]
rep_strand  = max(set(all_strands),key=all_strands.count) if all_strands else "?"
strand_lbl  = "역방향" if rep_strand=="-" else "정방향"

# 표시할 대표 위치
if has_hybrid:
    disp_chr   = prec_chr
    disp_pos   = f"{prec_pos:,} bp"
    disp_prec  = "&#177; 4 bp"
    disp_color = "#b71c1c"
else:
    disp_chr   = integration_chr
    disp_pos   = f"{int(integration_pos)/1e6:.2f} Mb"
    disp_prec  = "&#177; 수 kb"
    disp_color = "#283593"

# ── 정확한 breakpoint 섹션 (hybrid ref 결과 있을 때만) ─────────────────────
hybrid_section = ""
if has_hybrid:
    hybrid_section = (
        '<div style="background:#fff;border-radius:8px;padding:16px;margin-bottom:14px;border-left:4px solid #b71c1c;border:1px solid #dee2e6">'
        '<h2 style="margin-top:0;color:#b71c1c;font-size:1.05em">&#127919; 정확한 Breakpoint (Hybrid Reference 분석)</h2>'
        '<div style="display:grid;grid-template-columns:1fr 1fr;gap:14px">'
        '<div>'
        f'<table style="width:100%;border-collapse:collapse;font-size:0.86em">'
        f'<tr><td style="color:#666;padding:4px 6px;width:45%">Chromosome</td><td style="font-family:monospace;font-weight:600;padding:4px 6px">{prec_chr}</td></tr>'
        f'<tr><td style="color:#666;padding:4px 6px">Breakpoint</td><td style="font-family:monospace;font-weight:600;color:#b71c1c;padding:4px 6px">{prec_pos:,} bp</td></tr>'
        f'<tr><td style="color:#666;padding:4px 6px">Range</td><td style="font-family:monospace;padding:4px 6px">{bp_range}</td></tr>'
        f'<tr><td style="color:#666;padding:4px 6px">Supporting reads</td><td style="padding:4px 6px">{chimeric_n}개 chimeric reads</td></tr>'
        f'</table>'
        '</div>'
        '<div>'
        '<table style="border-collapse:collapse;font-size:0.85em;width:100%"><thead><tr style="background:#f8f9fa">'
        '<th style="padding:4px 7px;text-align:left;border-bottom:2px solid #dee2e6">방법</th>'
        '<th style="padding:4px 7px;text-align:left;border-bottom:2px solid #dee2e6">결과</th>'
        '<th style="padding:4px 7px;border-bottom:2px solid #dee2e6">정밀도</th>'
        '</tr></thead><tbody>'
        f'<tr><td style="padding:4px 7px">De novo (보간)</td><td style="font-family:monospace;padding:4px 7px">{disp_chr}:{int(integration_pos)/1e6:.2f} Mb</td>'
        f'<td style="padding:4px 7px"><span style="background:#fff3cd;color:#856404;padding:1px 6px;border-radius:8px;font-size:0.8em">&#177; 수 kb</span></td></tr>'
        f'<tr style="background:#f0fff0"><td style="padding:4px 7px"><strong>Hybrid ref</strong></td>'
        f'<td style="font-family:monospace;font-weight:600;color:#b71c1c;padding:4px 7px">{prec_chr}:{prec_pos:,}</td>'
        f'<td style="padding:4px 7px"><span style="background:#d4edda;color:#155724;padding:1px 6px;border-radius:8px;font-size:0.8em">&#177; 4 bp</span></td></tr>'
        '</tbody></table>'
        '</div>'
        '</div></div>'
    )

# ── 유전자 어노테이션 섹션 (annotation 결과 있을 때만) ─────────────────────
annot_section = ""
if has_annot:
    feat_color = {"exon":"#c62828","intron":"#2e7d32","intergenic":"#616161"}.get(feat_type,"#616161")
    exon_dist_txt = f"{exon_dist:,} bp" if exon_dist >= 0 else "N/A"
    annot_section = (
        '<div style="background:#fff;border-radius:8px;padding:16px;margin-bottom:14px;border-left:4px solid #2e7d32;border:1px solid #dee2e6">'
        '<h2 style="margin-top:0;color:#2e7d32;font-size:1.05em">&#127809; 유전자 어노테이션</h2>'
        '<div style="display:grid;grid-template-columns:auto 1fr;gap:14px;align-items:start">'
        f'<table style="border-collapse:collapse;font-size:0.86em;white-space:nowrap">'
        f'<tr><td style="color:#666;padding:4px 6px">유전자</td><td style="padding:4px 6px"><em><strong>{gene_name}</strong></em>'
        f' <span style="background:#d4edda;color:#155724;padding:1px 6px;border-radius:8px;font-size:0.79em">{gene_type}</span></td></tr>'
        f'<tr><td style="color:#666;padding:4px 6px">삽입 위치</td><td style="padding:4px 6px">'
        f'<span style="background:{feat_color};color:#fff;padding:1px 8px;border-radius:8px;font-size:0.82em">{feat_type}</span></td></tr>'
        f'<tr><td style="color:#666;padding:4px 6px">최근접 exon</td><td style="padding:4px 6px;font-family:monospace">{exon_dist_txt}</td></tr>'
        f'</table>'
        '<div style="background:#fff8e1;border-radius:6px;padding:10px;font-size:0.84em">'
        '<div style="font-weight:600;color:#f57f17;margin-bottom:5px">&#9888; 모니터링 권장</div>'
        f'<ul style="list-style:disc;padding-left:13px;color:#555;line-height:1.8">'
        f'<li><em>{gene_name}</em> mRNA 발현량 qRT-PCR 확인</li>'
        f'<li>삽입 부위({feat_type}) 관련 전사체 이상 여부 확인</li>'
        f'<li>표적 조직에서의 <em>{gene_name}</em> 발현 패턴 비교</li>'
        f'</ul></div>'
        '</div></div>'
    )

pos_label = f"{int(integration_pos)/1e6:.2f} Mb" if not has_hybrid else f"{prec_pos:,} bp"

html_out = (
    '<!DOCTYPE html>\n<html lang="ko">\n<head>\n<meta charset="UTF-8">\n'
    f'<title>Transgene Integration &#8211; {sample_id}</title>\n'
    '<style>\n'
    '* { box-sizing:border-box; margin:0; padding:0 }\n'
    'body { font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; background:#f0f2f5; color:#212529; line-height:1.5 }\n'
    '.wrap { max-width:1120px; margin:0 auto; padding:22px }\n'
    'h2 { font-size:1.08em; font-weight:600; margin:18px 0 8px }\n'
    'table td,table th { padding:5px 8px; border-bottom:1px solid #f0f0f0 }\n'
    'table tbody tr:hover { background:#fafafa }\n'
    'code { background:#f4f4f4; padding:1px 5px; border-radius:3px; font-size:0.89em }\n'
    '</style>\n</head>\n<body>\n<div class="wrap">\n\n'

    f'<div style="background:linear-gradient(135deg,#1a237e,#283593);color:#fff;border-radius:10px;padding:20px 24px;margin-bottom:16px">\n'
    f'  <div style="font-size:0.81em;opacity:.7;margin-bottom:3px">De Novo Assembly + Hybrid Reference &#8211; Transgene Integration Analysis</div>\n'
    f'  <div style="font-size:1.55em;font-weight:700">{sample_id} Integration Report</div>\n'
    f'  <div style="margin-top:7px;font-size:0.87em;opacity:.85">Transgene: <strong>{transgene_name}</strong> &nbsp;&middot;&nbsp; Generated: {today}</div>\n'
    '</div>\n\n'

    f'<div style="display:grid;grid-template-columns:repeat({"5" if has_annot else "4"},1fr);gap:9px;margin-bottom:14px">\n'
    f'  <div style="background:#fff;border-radius:8px;padding:13px;border-left:4px solid #1565c0"><div style="font-size:0.72em;color:#888;font-weight:600;text-transform:uppercase">통합 염색체</div><div style="font-size:1.35em;font-weight:700;color:#1565c0;margin-top:2px">{disp_chr}</div><div style="font-size:0.74em;color:#666;margin-top:1px">reference</div></div>\n'
    f'  <div style="background:#fff;border-radius:8px;padding:13px;border-left:4px solid {disp_color}"><div style="font-size:0.72em;color:#888;font-weight:600;text-transform:uppercase">통합 위치</div><div style="font-size:1.2em;font-weight:700;color:{disp_color};margin-top:2px;font-family:monospace">{pos_label}</div><div style="font-size:0.74em;color:#666;margin-top:1px">{disp_prec}</div></div>\n'
    f'  <div style="background:#fff;border-radius:8px;padding:13px;border-left:4px solid #e74c3c"><div style="font-size:0.72em;color:#888;font-weight:600;text-transform:uppercase">탠덤 카피 수</div><div style="font-size:1.35em;font-weight:700;color:#e74c3c;margin-top:2px">{max_copies}카피</div><div style="font-size:0.74em;color:#666;margin-top:1px">hap1 &amp; hap2</div></div>\n'
    f'  <div style="background:#fff;border-radius:8px;padding:13px;border-left:4px solid #6a1b9a"><div style="font-size:0.72em;color:#888;font-weight:600;text-transform:uppercase">삽입 방향</div><div style="font-size:1.15em;font-weight:700;color:#6a1b9a;margin-top:2px">{strand_lbl}</div><div style="font-size:0.74em;color:#666;margin-top:1px">Reverse complement</div></div>\n'
    + (f'  <div style="background:#fff;border-radius:8px;padding:13px;border-left:4px solid #2e7d32"><div style="font-size:0.72em;color:#888;font-weight:600;text-transform:uppercase">삽입 유전자</div><div style="font-size:1.2em;font-weight:700;color:#2e7d32;margin-top:2px"><em>{gene_name}</em></div><div style="font-size:0.74em;color:#666;margin-top:1px">{feat_type}</div></div>\n' if has_annot else '')
    + '</div>\n\n'

    + hybrid_section
    + annot_section

    + '<div style="background:#fff;border-radius:8px;padding:16px;margin-bottom:14px;border:1px solid #dee2e6">\n'
    '<h2 style="margin-top:0">&#129516; De Novo Assembly 분석</h2>\n'
    '<div style="display:grid;grid-template-columns:auto 1fr;gap:14px;align-items:start;margin-bottom:14px">\n'
    f'<table style="border-collapse:collapse;font-size:0.85em;white-space:nowrap"><thead><tr style="background:#f8f9fa"><th style="padding:4px 7px;text-align:left;border-bottom:2px solid #dee2e6">항목</th><th style="padding:4px 7px;text-align:center;border-bottom:2px solid #dee2e6">hap1</th><th style="padding:4px 7px;text-align:center;border-bottom:2px solid #dee2e6">hap2</th></tr></thead><tbody>{cmp_rows}</tbody></table>\n'
    '<div style="background:#fff8e1;border-radius:6px;padding:11px;font-size:0.84em"><div style="font-weight:600;color:#f57f17;margin-bottom:5px">&#9888; 주의사항</div><ul style="list-style:disc;padding-left:13px;color:#555;line-height:1.8"><li>De novo ref 위치는 선형 보간 추정값 &#8212; hybrid ref breakpoint가 더 정확</li><li>mapq 0/1 : 다중 카피로 인한 정상적 낮은 mapq</li><li>hap2 복수 블록 → 추가 구조변이 가능</li></ul></div>\n'
    '</div>\n\n'

    '<h3 style="font-size:0.95em;font-weight:600;margin-bottom:8px">&#128204; Hap1</h3>\n'
    f'{hap1_html}\n'
    '<h3 style="font-size:0.95em;font-weight:600;margin:14px 0 8px">&#128204; Hap2</h3>\n'
    f'{hap2_html}\n'
    '</div>\n\n'

    f'<div style="text-align:center;color:#aaa;font-size:0.76em;padding:12px">Generated by transgene_integration WDL pipeline &middot; {today}</div>\n'
    '</div>\n</body>\n</html>\n'
)

with open(f"{out_prefix}.integration_report.html", "w") as f:
    f.write(html_out)

print(f"Reports: {out_prefix}.integration_report.txt / .tsv / .html")
PYEOF
  >>>

  output {
    File   report_txt           = "~{out_prefix}.integration_report.txt"
    File   report_html          = "~{out_prefix}.integration_report.html"
    File   integration_tsv      = "~{out_prefix}.integration_sites.tsv"
    String stat_hap1_hit_count  = read_string("stat_hap1_hit_count.txt")
    String stat_hap2_hit_count  = read_string("stat_hap2_hit_count.txt")
    String stat_integration_chr = read_string("stat_integration_chr.txt")
    String stat_integration_pos = read_string("stat_integration_pos.txt")
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
