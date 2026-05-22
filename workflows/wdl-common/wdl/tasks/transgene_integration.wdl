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
    description: "Parse alignment results and generate transgene integration site report (text, TSV, HTML)"
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
    report_html:          { name: "Integration site report (HTML)" }
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
from datetime import date

sample_id      = "~{sample_id}"
transgene_name = "~{transgene_name}"
out_prefix     = "~{out_prefix}"
today          = date.today().isoformat()

def parse_tg_paf(paf_file, min_match=500):
    hits = defaultdict(list)
    if not os.path.exists(paf_file) or os.path.getsize(paf_file) == 0:
        return hits
    with open(paf_file) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) < 12:
                continue
            matches = int(p[9])
            if matches < min_match:
                continue
            hits[p[5]].append({
                "tg_name": p[0], "tg_start": int(p[2]), "tg_end": int(p[3]),
                "tg_len": int(p[1]), "contig_len": int(p[6]),
                "c_start": int(p[7]), "c_end": int(p[8]),
                "strand": p[4], "matches": matches, "mapq": int(p[11]),
            })
    return hits

def parse_ref_paf(paf_file, min_match=1000):
    alns = defaultdict(list)
    if not os.path.exists(paf_file) or os.path.getsize(paf_file) == 0:
        return alns
    with open(paf_file) as f:
        for line in f:
            p = line.strip().split('\t')
            if len(p) < 12:
                continue
            if int(p[9]) < min_match:
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

# ── 데이터 파싱 (text/TSV/HTML 공용) ──────────────────────────────────────────
hap_data = {}
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
lines.append("=" * 70)

hap1_count = 0
hap2_count = 0
tsv_rows   = ["\t".join([
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
    if hap == "hap1":
        hap1_count = count
    else:
        hap2_count = count
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

# ── HTML Report ────────────────────────────────────────────────────────────────
def pct_bar(value, color="#4CAF50"):
    w = min(value, 100)
    return (f'<div style="background:#e9ecef;border-radius:4px;height:8px;width:120px;'
            f'display:inline-block;vertical-align:middle">'
            f'<div style="background:{color};width:{w:.1f}%;height:100%;border-radius:4px">'
            f'</div></div> <span style="font-size:0.84em;color:#555">{value:.1f}%</span>')

def strand_badge(s):
    if s == "-":
        return ('<span style="background:#e74c3c;color:#fff;padding:1px 7px;'
                'border-radius:10px;font-size:0.8em">&#8722; 역방향</span>')
    return ('<span style="background:#27ae60;color:#fff;padding:1px 7px;'
            'border-radius:10px;font-size:0.8em">+ 정방향</span>')

def make_contig_svg(contig_len, copies):
    W, H, scale = 680, 56, 680 / contig_len
    parts = [
        f'<rect x="0" y="20" width="{W}" height="16" rx="3" fill="#b0bec5"/>',
        f'<text x="2" y="16" font-size="11" fill="#546e7a">host</text>',
        f'<text x="{W-2}" y="16" font-size="11" fill="#546e7a" text-anchor="end">{contig_len:,} bp</text>',
    ]
    for i, cp in enumerate(copies):
        x1 = cp["c_start"] * scale
        w  = max((cp["c_end"] - cp["c_start"]) * scale, 2)
        parts += [
            f'<rect x="{x1:.1f}" y="14" width="{w:.1f}" height="28" rx="2" fill="#e74c3c" opacity="0.9"/>',
            f'<text x="{x1+w/2:.1f}" y="31" font-size="10" fill="white" text-anchor="middle" font-weight="bold">TG {i+1}</text>',
            f'<text x="{x1:.1f}" y="50" font-size="9" fill="#666">{cp["c_start"]:,}</text>',
        ]
    return (f'<svg width="{W}" height="{H}" style="display:block;margin:6px 0">'
            + "".join(parts) + '</svg>')

def build_hap_section(hap, contigs):
    color = "#1565c0" if hap == "hap1" else "#6a1b9a"
    label = hap.upper()
    if not contigs:
        return f'<p style="color:#888">No chimeric contigs detected in {label}.</p>'
    blocks = []
    for c in contigs:
        copies, gaps, ref_alns = c["copies"], c["gaps"], c["ref_alns"]
        copy_rows = ""
        for i, cp in enumerate(copies, 1):
            copy_rows += (
                f'<tr>'
                f'<td style="text-align:center;font-weight:600">{i}</td>'
                f'<td style="font-family:monospace;white-space:nowrap">{cp["c_start"]:,}&#8211;{cp["c_end"]:,}</td>'
                f'<td>{pct_bar(cp["tg_cov"])}</td>'
                f'<td>{strand_badge(cp["strand"])}</td>'
                f'<td style="font-family:monospace">{cp["chrom"]}:{cp["ref_pos"]:,}</td>'
                f'<td style="text-align:center">{cp["ref_mapq"]}</td>'
                f'</tr>'
            )
            if i <= len(gaps):
                copy_rows += (f'<tr><td colspan="6" style="text-align:center;color:#888;'
                              f'font-size:0.84em;padding:2px">&#8595; gap {gaps[i-1]:,} bp</td></tr>')
        ref_rows = "".join(
            f'<tr>'
            f'<td style="font-family:monospace;white-space:nowrap">{a["chrom"]}:{a["t_start"]:,}&#8211;{a["t_end"]:,}</td>'
            f'<td style="font-family:monospace;white-space:nowrap">{a["q_start"]:,}&#8211;{a["q_end"]:,}</td>'
            f'<td style="text-align:center">{a["strand"]}</td>'
            f'<td style="text-align:center">{a["mapq"]}</td>'
            f'</tr>'
            for a in ref_alns
        )
        n_ref = len(ref_alns)
        ref_note = (
            '<p style="color:#27ae60;font-size:0.84em;margin-top:6px">&#10003; 단일 블록 (단순 삽입)</p>'
            if n_ref <= 1 else
            f'<p style="color:#e67e22;font-size:0.84em;margin-top:6px">&#9888; {n_ref}개 블록 &#8594; 추가 구조변이 가능</p>'
        )
        blocks.append(
            f'<div style="background:#fff;border:1px solid #dee2e6;border-radius:8px;padding:18px;margin-bottom:14px">'
            f'<div style="display:flex;align-items:center;gap:10px;margin-bottom:10px">'
            f'<span style="background:{color};color:#fff;padding:2px 10px;border-radius:12px;font-size:0.84em">{label}</span>'
            f'<code style="font-size:0.95em;font-weight:600">{c["name"]}</code>'
            f'<span style="color:#888;font-size:0.88em">{c["contig_len"]:,} bp &middot; {len(copies)}카피 탠덤</span>'
            f'</div>'
            f'{make_contig_svg(c["contig_len"], copies)}'
            f'<div style="display:grid;grid-template-columns:1fr 1fr;gap:14px;margin-top:12px">'
            f'<div>'
            f'<div style="font-size:0.78em;color:#666;font-weight:600;margin-bottom:4px">트랜스진 카피 상세</div>'
            f'<table style="width:100%;border-collapse:collapse;font-size:0.86em">'
            f'<thead><tr style="background:#f8f9fa">'
            f'<th style="padding:4px 7px;text-align:center;border-bottom:2px solid #dee2e6">#</th>'
            f'<th style="padding:4px 7px;text-align:left;border-bottom:2px solid #dee2e6">Contig 위치</th>'
            f'<th style="padding:4px 7px;text-align:left;border-bottom:2px solid #dee2e6">TG 커버리지</th>'
            f'<th style="padding:4px 7px;text-align:left;border-bottom:2px solid #dee2e6">방향</th>'
            f'<th style="padding:4px 7px;text-align:left;border-bottom:2px solid #dee2e6">Ref 위치 (추정)</th>'
            f'<th style="padding:4px 7px;text-align:center;border-bottom:2px solid #dee2e6">mapq</th>'
            f'</tr></thead>'
            f'<tbody>{copy_rows}</tbody>'
            f'</table>'
            f'</div>'
            f'<div>'
            f'<div style="font-size:0.78em;color:#666;font-weight:600;margin-bottom:4px">Reference 정렬 블록</div>'
            f'<table style="width:100%;border-collapse:collapse;font-size:0.86em">'
            f'<thead><tr style="background:#f8f9fa">'
            f'<th style="padding:4px 7px;text-align:left;border-bottom:2px solid #dee2e6">Ref 좌표</th>'
            f'<th style="padding:4px 7px;text-align:left;border-bottom:2px solid #dee2e6">Contig 좌표</th>'
            f'<th style="padding:4px 7px;text-align:center;border-bottom:2px solid #dee2e6">Strand</th>'
            f'<th style="padding:4px 7px;text-align:center;border-bottom:2px solid #dee2e6">mapq</th>'
            f'</tr></thead>'
            f'<tbody>{ref_rows}</tbody>'
            f'</table>'
            f'{ref_note}'
            f'</div>'
            f'</div>'
            f'</div>'
        )
    return "\n".join(blocks)

hap1_html_content = build_hap_section("hap1", hap_data.get("hap1", []))
hap2_html_content = build_hap_section("hap2", hap_data.get("hap2", []))

rep_chrom   = integration_chr
rep_pos_int = int(integration_pos) if integration_pos != "0" else 0
pos_label   = f"~{rep_pos_int/1e6:.2f} Mb" if rep_pos_int else "N/A"

all_copies  = [cp for hap in ["hap1","hap2"] for c in hap_data.get(hap,[]) for cp in c["copies"]]
max_copies  = max((len(c["copies"]) for hap in ["hap1","hap2"] for c in hap_data.get(hap,[])), default=0)
all_strands = [cp["strand"] for cp in all_copies]
rep_strand  = max(set(all_strands), key=all_strands.count) if all_strands else "?"
strand_lbl  = "역방향" if rep_strand == "-" else "정방향"
strand_desc = "Reverse complement" if rep_strand == "-" else "Forward"

def hap_val(hap, fn, default="N/A"):
    d = hap_data.get(hap, [])
    return fn(d[0]) if d else default

def ref_badge(hap):
    d = hap_data.get(hap, [])
    if not d: return "N/A"
    n = len(d[0]["ref_alns"])
    if n <= 1:
        return f'<span style="background:#d4edda;color:#155724;padding:2px 7px;border-radius:10px;font-size:0.8em">{n} (단순)</span>'
    return f'<span style="background:#fff3cd;color:#856404;padding:2px 7px;border-radius:10px;font-size:0.8em">{n} (복잡)</span>'

cmp_rows = "".join(
    f'<tr>'
    f'<td style="color:#555;padding:5px 8px">{r}</td>'
    f'<td style="text-align:center;padding:5px 8px">{v1}</td>'
    f'<td style="text-align:center;padding:5px 8px">{v2}</td>'
    f'</tr>'
    for r, v1, v2 in [
        ("Contig",
         hap_val("hap1", lambda c: c["name"]),
         hap_val("hap2", lambda c: c["name"])),
        ("Contig 크기",
         hap_val("hap1", lambda c: f'{c["contig_len"]:,} bp'),
         hap_val("hap2", lambda c: f'{c["contig_len"]:,} bp')),
        ("카피 수",
         hap_val("hap1", lambda c: str(len(c["copies"])), "0"),
         hap_val("hap2", lambda c: str(len(c["copies"])), "0")),
        ("삽입 방향",
         hap_val("hap1", lambda c: "역방향 (&#8722;)" if c["copies"] and c["copies"][0]["strand"]=="-" else "정방향 (+)"),
         hap_val("hap2", lambda c: "역방향 (&#8722;)" if c["copies"] and c["copies"][0]["strand"]=="-" else "정방향 (+)")),
        ("Ref 정렬 블록", ref_badge("hap1"), ref_badge("hap2")),
    ]
)

html_out = (
    '<!DOCTYPE html>\n'
    '<html lang="ko">\n'
    '<head>\n'
    '<meta charset="UTF-8">\n'
    f'<title>Transgene Integration Report &#8211; {sample_id}</title>\n'
    '<style>\n'
    '* { box-sizing:border-box; margin:0; padding:0 }\n'
    'body { font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;\n'
    '       background:#f0f2f5; color:#212529; line-height:1.5 }\n'
    '.wrap { max-width:1100px; margin:0 auto; padding:24px }\n'
    'h2 { font-size:1.1em; font-weight:600; margin:20px 0 10px }\n'
    'table td,table th { padding:6px 10px; border-bottom:1px solid #f0f0f0 }\n'
    'table tbody tr:hover { background:#fafafa }\n'
    'code { background:#f4f4f4; padding:1px 5px; border-radius:3px; font-size:0.9em }\n'
    '</style>\n'
    '</head>\n'
    '<body>\n'
    '<div class="wrap">\n'
    '\n'
    '<div style="background:#1a237e;color:#fff;border-radius:10px;padding:22px 26px;margin-bottom:18px">\n'
    '  <div style="font-size:0.83em;opacity:.7;margin-bottom:3px">De Novo Assembly &#8211; Transgene Integration Analysis</div>\n'
    f'  <div style="font-size:1.55em;font-weight:700">{sample_id} Integration Report</div>\n'
    f'  <div style="margin-top:8px;font-size:0.88em;opacity:.85">Transgene: <strong>{transgene_name}</strong> &nbsp;&middot;&nbsp; Generated: {today}</div>\n'
    '</div>\n'
    '\n'
    '<div style="display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin-bottom:18px">\n'
    f'  <div style="background:#fff;border-radius:8px;padding:14px;border-left:4px solid #1565c0"><div style="font-size:0.75em;color:#888;font-weight:600;text-transform:uppercase">통합 염색체</div><div style="font-size:1.35em;font-weight:700;color:#1565c0;margin-top:3px">{rep_chrom}</div><div style="font-size:0.78em;color:#666;margin-top:1px">GRCm39</div></div>\n'
    f'  <div style="background:#fff;border-radius:8px;padding:14px;border-left:4px solid #283593"><div style="font-size:0.75em;color:#888;font-weight:600;text-transform:uppercase">통합 위치 (추정)</div><div style="font-size:1.35em;font-weight:700;color:#283593;margin-top:3px">{pos_label}</div><div style="font-size:0.78em;color:#666;margin-top:1px">&#177; 수 kb 오차</div></div>\n'
    f'  <div style="background:#fff;border-radius:8px;padding:14px;border-left:4px solid #e74c3c"><div style="font-size:0.75em;color:#888;font-weight:600;text-transform:uppercase">탠덤 카피 수</div><div style="font-size:1.35em;font-weight:700;color:#e74c3c;margin-top:3px">{max_copies}카피</div><div style="font-size:0.78em;color:#666;margin-top:1px">hap1 &amp; hap2 공통</div></div>\n'
    f'  <div style="background:#fff;border-radius:8px;padding:14px;border-left:4px solid #6a1b9a"><div style="font-size:0.75em;color:#888;font-weight:600;text-transform:uppercase">삽입 방향</div><div style="font-size:1.35em;font-weight:700;color:#6a1b9a;margin-top:3px">{strand_lbl}</div><div style="font-size:0.78em;color:#666;margin-top:1px">{strand_desc}</div></div>\n'
    '</div>\n'
    '\n'
    '<div style="background:#fff;border-radius:8px;padding:18px;margin-bottom:18px;border:1px solid #dee2e6">\n'
    '  <h2 style="margin-top:0">hap1 vs hap2 비교</h2>\n'
    '  <div style="display:grid;grid-template-columns:auto 1fr;gap:16px;align-items:start">\n'
    '    <table style="border-collapse:collapse;font-size:0.87em;white-space:nowrap">\n'
    '      <thead><tr style="background:#f8f9fa">\n'
    '        <th style="padding:5px 8px;text-align:left;border-bottom:2px solid #dee2e6">항목</th>\n'
    '        <th style="padding:5px 8px;text-align:center;border-bottom:2px solid #dee2e6">hap1</th>\n'
    '        <th style="padding:5px 8px;text-align:center;border-bottom:2px solid #dee2e6">hap2</th>\n'
    '      </tr></thead>\n'
    f'      <tbody>{cmp_rows}</tbody>\n'
    '    </table>\n'
    '    <div style="background:#fff8e1;border-radius:6px;padding:12px;font-size:0.86em">\n'
    '      <div style="font-weight:600;margin-bottom:6px;color:#f57f17">&#9888; 해석 시 주의사항</div>\n'
    '      <ul style="list-style:disc;padding-left:14px;color:#555;line-height:1.8">\n'
    '        <li>Reference 위치는 선형 보간 추정값 (&#177; 수 kb 오차)</li>\n'
    '        <li>정확한 breakpoint는 CIGAR 파싱 또는 IGV 확인 권장</li>\n'
    '        <li>mapq 0/1 : 다중 카피로 인한 정상적인 낮은 mapq</li>\n'
    '        <li>삽입 부위 주변 유전자 어노테이션 별도 확인 필요</li>\n'
    '      </ul>\n'
    '    </div>\n'
    '  </div>\n'
    '</div>\n'
    '\n'
    '<h2>&#128204; Hap1 상세</h2>\n'
    f'{hap1_html_content}\n'
    '\n'
    '<h2>&#128204; Hap2 상세</h2>\n'
    f'{hap2_html_content}\n'
    '\n'
    f'<div style="text-align:center;color:#aaa;font-size:0.78em;padding:14px">Generated by transgene_integration WDL pipeline &middot; {today}</div>\n'
    '</div>\n'
    '</body>\n'
    '</html>\n'
)

with open(f"{out_prefix}.integration_report.html", "w") as fh:
    fh.write(html_out)

print(f"HTML: {out_prefix}.integration_report.html")
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
