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
    # -p 0.05 -N 50: asm20 기본값(-p 0.8)은 부분 카피·짧은 카피 정렬을 secondary로 버림 → tandem array의 partial copy 누락
    minimap2 -cx asm20 --cs -p 0.05 -N 50 -t ~{threads} \
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
    description: "Chimeric contig PAF에서 TG 블록 양쪽의 host 접합 후보를 모두 수집 — 후보 클러스터링, 다중 위치·염색체 전환 경고"
  }

  parameter_meta {
    hap1_ref_paf:       { name: "PAF: hap1 chimeric contigs vs reference" }
    hap2_ref_paf:       { name: "PAF: hap2 chimeric contigs vs reference" }
    hap1_tg_paf:        { name: "PAF: transgene vs hap1 assembly" }
    hap2_tg_paf:        { name: "PAF: transgene vs hap2 assembly" }
    min_match_bp:       { name: "Minimum matching bases for a TG hit (default: 500)" }
    min_flank_bp:       { name: "Minimum host alignment length to call a flank (default: 2000)" }
    transgene_fasta:    { name: "Transgene FASTA (checked for duplicated / circular ends)" }
    bac_region:         { name: "[BAC] Reference interval of the BAC insert (chr:start-end); alignments inside it are construct, not host" }
    runtime_attributes: { name: "Runtime attribute structure" }
    integration_chr:    { name: "Primary junction cluster chromosome" }
    integration_pos:    { name: "Primary junction cluster position (bp, 1-based last host base)" }
    junctions_tsv:      { name: "All host/TG junction candidates (hap, contig, TG block, side, ref position, gap, mapq, cluster)" }
    qc_warnings:        { name: "Semicolon-separated warnings (MULTI_SITE, SPLIT_JUNCTION, UNALIGNED_GAP, ONE_SIDED, TG_HOST_HOMOLOGY, DUPLICATED_ENDS)" }
    loci_txt:           { name: "One candidate locus per line (chr<TAB>pos): every junction cluster, for read extraction" }
  }

  input {
    File   hap1_ref_paf
    File   hap2_ref_paf
    File   hap1_tg_paf
    File   hap2_tg_paf
    Int    min_match_bp = 500
    Int    min_flank_bp = 2000
    File?  transgene_fasta
    String construct_type = "plasmid"
    String bac_region = ""
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail
    python3 << 'PYEOF'
from collections import defaultdict

BAC = "~{bac_region}".replace(",", "")
bac = None
if BAC:
    c, se = BAC.split(":")
    bac = (c, int(se.split("-")[0]), int(se.split("-")[1]))

def in_bac(chrom, s, e):
    return bac is not None and chrom == bac[0] and min(e, bac[2]) - max(s, bac[1]) > 0

MIN_MATCH  = ~{min_match_bp}
MIN_FLANK  = ~{min_flank_bp}
BLOCK_GAP  = 5000      # TG hits closer than this on a contig = one TG block (tandem array)
FLANK_TOL  = 2000      # host alignment may overlap the TG block edge by this much
MAX_GAP    = 20000     # a host alignment further than this from the TG block is not its flank
CLUSTER_BP = 1000      # junction candidates within this distance = one junction
SITE_BP    = 1000000   # junction clusters within this distance = one integration site

SPLIT_INDEL = 1000     # contig alignments are split at indels this large (asm5 often keeps a whole insertion inside one alignment)

def read_paf(path):
    rows = []
    try:
        fh = open(path)
    except OSError:
        return rows
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 12:
            continue
        cs = next((x[5:] for x in p[12:] if x.startswith("cs:Z:")), "")
        rows.append(dict(q=p[0], qlen=int(p[1]), qs=int(p[2]), qe=int(p[3]), strand=p[4], t=p[5],
                         ts=int(p[7]), te=int(p[8]), match=int(p[9]), mapq=int(p[11]), cs=cs))
    return rows

def split_at_indels(a):
    """Split one alignment into gap-free pieces at insertions/deletions >= SPLIT_INDEL (needs the cs tag).
    cs runs along the reference; on '-' strand the query coordinate runs backwards."""
    if not a["cs"]:
        return [a]
    import re as _re
    pieces, t, q = [], a["ts"], (a["qs"] if a["strand"] == "+" else a["qe"])
    t0, q0 = t, q
    def close(t1, q1):
        if t1 > t0:
            qs_, qe_ = (q0, q1) if a["strand"] == "+" else (q1, q0)
            pieces.append(dict(a, ts=t0, te=t1, qs=qs_, qe=qe_, cs=""))
    step = 1 if a["strand"] == "+" else -1
    for op, val in _re.findall(r"([:*+\-~])([0-9]+|[a-z]+)", a["cs"]):
        if op == ":":
            n = int(val); t += n; q += step * n
        elif op == "*":
            t += 1; q += step
        elif op == "+":
            if len(val) >= SPLIT_INDEL:
                close(t, q); q += step * len(val); t0, q0 = t, q
            else:
                q += step * len(val)
        elif op == "-":
            if len(val) >= SPLIT_INDEL:
                close(t, q); t += len(val); t0, q0 = t, q
            else:
                t += len(val)
    close(t, q)
    return pieces or [a]

def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))

junctions, warnings = [], []
if "~{construct_type}".upper() == "BAC" and not bac:
    warnings.append("BAC_REGION_MISSING construct_type=BAC but bac_region is empty: the BAC's native locus will be "
                    "mistaken for the integration site. Give the BAC interval (chr:start-end) and re-run")

# construct FASTA check: a sequence whose start re-occurs inside it was cut from a circular map or has duplicated ends
tg_fa = "~{default="" transgene_fasta}"
if tg_fa:
    name, seq = None, []
    for line in open(tg_fa):
        if line.startswith(">"):
            name = line[1:].split()[0]
        else:
            seq.append(line.strip().upper())
    seq = "".join(seq)
    k = 100
    if len(seq) > 4 * k and seq.find(seq[:k], k) > 0:
        warnings.append(f"DUPLICATED_ENDS {name}: the first {k} bp re-occur at {seq.find(seq[:k], k) + 1} "
                        "(circular map opened twice or duplicated backbone; copy counts and copy boundaries will be wrong)")

for hap, tg_paf, ref_paf in (("hap1", "~{hap1_tg_paf}", "~{hap1_ref_paf}"),
                             ("hap2", "~{hap2_tg_paf}", "~{hap2_ref_paf}")):
    hits = defaultdict(list)                       # contig -> [(start, end)] of TG hits (PAF target = contig)
    for r in read_paf(tg_paf):
        if r["match"] >= MIN_MATCH:
            hits[r["t"]].append((r["ts"], r["te"]))
    refs = defaultdict(list)                       # contig -> ref alignments (PAF query = contig)
    for r in read_paf(ref_paf):
        refs[r["q"]].extend(split_at_indels(r))
    for ctg, iv in hits.items():
        iv.sort()
        blocks = [list(iv[0])]
        for s, e in iv[1:]:
            if s - blocks[-1][1] <= BLOCK_GAP:
                blocks[-1][1] = max(blocks[-1][1], e)
            else:
                blocks.append([s, e])
        alns = [a for a in refs.get(ctg, []) if a["mapq"] >= 5 and a["qe"] - a["qs"] >= MIN_FLANK]
        if bac:
            # BAC mode: contig pieces that align inside the BAC interval are construct, never host flank
            n0 = len(alns)
            alns = [a for a in alns if not in_bac(a["t"], a["ts"], a["te"])]
            if len(alns) < n0:
                warnings.append(f"BAC_REGION {hap}:{ctg}: {n0 - len(alns)} alignment(s) inside {BAC} treated as construct")
        for b0, b1 in blocks:
            bid = f"{hap}:{ctg}:{b0 + 1}-{b1}"
            # construct sequence that itself aligns to the reference (BAC insert, genomic promoter/UTR pieces)
            homol = [a for a in alns if overlap(a["qs"], a["qe"], b0, b1) >= 0.5 * (b1 - b0)]
            if homol:
                # the TG hit is host sequence in its native place (e.g. a polyA/promoter/UTR piece of the construct,
                # or a BAC insert): report it, but do not call integration junctions from it
                for a in homol:
                    # project the block onto the reference (ignores indels; approximate)
                    rs = a["ts"] + (b0 - a["qs"]) if a["strand"] == "+" else a["te"] - (b1 - a["qs"])
                    rs = max(rs, a["ts"])
                    warnings.append(f"TG_HOST_HOMOLOGY {bid} ({b1 - b0} bp) sits in contiguous host sequence at approx. {a['t']}:{rs + 1}-{rs + b1 - b0} "
                                    "(construct carries host-derived sequence; not used for junctions; if the construct is a BAC, "
                                    "its native locus is NOT the integration site)")
                continue
            host = [a for a in alns if overlap(a["qs"], a["qe"], b0, b1) < 0.5 * (a["qe"] - a["qs"])]
            left = [a for a in host if a["qe"] <= b0 + FLANK_TOL and b0 - a["qe"] <= MAX_GAP]
            right = [a for a in host if a["qs"] >= b1 - FLANK_TOL and a["qs"] - b1 <= MAX_GAP]
            sides = {}
            if left:
                a = max(left, key=lambda x: x["qe"])
                # last host base before the TG block, in reference coordinates (1-based)
                pos = a["te"] if a["strand"] == "+" else a["ts"] + 1
                sides["left"] = (a, pos, b0 - a["qe"])
            if right:
                a = min(right, key=lambda x: x["qs"])
                pos = a["ts"] + 1 if a["strand"] == "+" else a["te"]
                sides["right"] = (a, pos, a["qs"] - b1)
            for side, (a, pos, gap) in sides.items():
                junctions.append(dict(hap=hap, contig=ctg, tg_block=bid, block_len=b1 - b0, side=side, chrom=a["t"], pos=pos,
                                      gap_bp=gap, mapq=a["mapq"], aln_len=a["qe"] - a["qs"], strand=a["strand"]))
                if gap > 500:
                    warnings.append(f"UNALIGNED_GAP {bid} {side}: {gap} bp between TG and host alignment (linker / unmapped sequence)")
            if len(sides) == 1:
                warnings.append(f"ONE_SIDED {bid}: host flank found only on the {next(iter(sides))} side")
            if len(sides) == 2:
                (al, pl, _), (ar, pr, _) = sides["left"], sides["right"]
                cl, cr = al["t"], ar["t"]
                if cl != cr or abs(pl - pr) > SITE_BP:
                    warnings.append(f"SPLIT_JUNCTION {bid}: left flank {cl}:{pl} vs right flank {cr}:{pr} "
                                    "(chromosome switch or large rearrangement at the integration)")

# cluster junction candidates (both haplotypes) -> junctions -> sites
junctions.sort(key=lambda j: (j["chrom"], j["pos"]))
clusters = []
for j in junctions:
    if clusters and clusters[-1]["chrom"] == j["chrom"] and j["pos"] - clusters[-1]["end"] <= CLUSTER_BP:
        c = clusters[-1]; c["members"].append(j); c["end"] = j["pos"]
    else:
        clusters.append(dict(chrom=j["chrom"], start=j["pos"], end=j["pos"], members=[j]))
for i, c in enumerate(clusters):
    c["id"] = f"J{i + 1}"
    for j in c["members"]:
        j["cluster"] = c["id"]
sites = []
for c in clusters:
    if sites and sites[-1]["chrom"] == c["chrom"] and c["start"] - sites[-1]["end"] <= SITE_BP:
        sites[-1]["clusters"].append(c); sites[-1]["end"] = c["end"]
    else:
        sites.append(dict(chrom=c["chrom"], start=c["start"], end=c["end"], clusters=[c]))
if len(sites) > 1:
    desc = ", ".join(f"{s['chrom']}:{s['start']}-{s['end']} ({sum(len(c['members']) for c in s['clusters'])} candidates)" for s in sites)
    warnings.append(f"MULTI_SITE {len(sites)} candidate integration sites: {desc}")

chrom, pos = "unknown", "0"
if clusters:
    # most haplotypes, most candidates, largest TG block (the real array), left side first
    best = max(clusters, key=lambda c: (len({m['hap'] for m in c['members']}), len(c["members"]),
                                         max(m["block_len"] for m in c["members"]),
                                         any(m["side"] == "left" for m in c["members"]),
                                         max(m["mapq"] for m in c["members"])))
    chrom = best["chrom"]
    pos = str(sorted(m["pos"] for m in best["members"])[len(best["members"]) // 2])

cols = ["cluster", "hap", "contig", "tg_block", "block_len", "side", "chrom", "pos", "gap_bp", "mapq", "aln_len", "strand"]
with open("junction_candidates.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for j in junctions:
        f.write("\t".join(str(j[c]) for c in cols) + "\n")
with open("loci.txt", "w") as f:
    for c in clusters:
        f.write(f"{c['chrom']}\t{c['start']}\n")
    if bac:
        f.write(f"{bac[0]}\t{bac[1]}\n{bac[0]}\t{bac[2]}\n")
open("integration_chr.txt", "w").write(chrom)
open("integration_pos.txt", "w").write(pos)
open("qc_warnings.txt", "w").write("; ".join(dict.fromkeys(warnings)) if warnings else "none")
print(f"Junction candidates: {len(junctions)} in {len(clusters)} clusters / {len(sites)} sites")
for c in clusters:
    print(f"  {c['id']} {c['chrom']}:{c['start']}-{c['end']}  n={len(c['members'])}  "
          f"sides={sorted({m['side'] for m in c['members']})}  haps={sorted({m['hap'] for m in c['members']})}")
print(f"Primary: {chrom}:{pos}")
for w in dict.fromkeys(warnings):
    print("WARNING:", w)
PYEOF
  >>>

  output {
    String integration_chr = read_string("integration_chr.txt")
    String integration_pos = read_string("integration_pos.txt")
    File   junctions_tsv   = "junction_candidates.tsv"
    String qc_warnings     = read_string("qc_warnings.txt")
    File   loci_txt        = "loci.txt"
  }

  runtime {
    docker:                "quay.io/biocontainers/busco:5.7.1--pyhdfd78af_1"
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
    description: "Whole reference genome + transgene = hybrid reference (reads may join the construct to any chromosome)"
  }

  parameter_meta {
    ref_fasta:          { name: "Reference genome FASTA" }
    transgene_fasta:    { name: "Transgene sequence FASTA" }
    out_prefix:         { name: "Output prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    hybrid_ref_fa:      { name: "Hybrid FASTA: all reference sequences + transgene" }
    hybrid_ref_fai:     { name: "FAI index of hybrid FASTA" }
    tg_seqname:         { name: "Transgene sequence name in the hybrid reference" }
  }

  input {
    File   ref_fasta
    File   transgene_fasta
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  # 이전 버전은 추정 염색체 하나만 넣어서, 다른 염색체의 접합(TG-6102 chr8, TG-6903 J2)을 볼 수 없었음
  Int disk_size = ceil(size(ref_fasta, "GB") * 3 + 5)

  command <<<
    set -euo pipefail
    cat "~{ref_fasta}" > "~{out_prefix}.hybrid_ref.fa"
    awk '/^>/{gsub(/ .*/,""); print; next} {print}' "~{transgene_fasta}" >> "~{out_prefix}.hybrid_ref.fa"
    samtools faidx "~{out_prefix}.hybrid_ref.fa"
    awk '/^>/{gsub(/^>/,""); print $1; exit}' "~{transgene_fasta}" > tg_seqname.txt
    echo "Hybrid ref: $(wc -l < ~{out_prefix}.hybrid_ref.fa.fai) sequences (+ $(cat tg_seqname.txt))"
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
    loci_txt:           { name: "Candidate loci (chr<TAB>pos per line) from extract_integration_coords" }
    window_bp:          { name: "Extraction window on each side of every locus (default: 300000)" }
    out_prefix:         { name: "Output prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    region_reads_fq:    { name: "Extracted reads in FASTQ.gz" }
    read_count:         { name: "Number of extracted reads" }
  }

  input {
    File   haplotagged_bam
    File   haplotagged_bai
    File   loci_txt
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

    # 모든 후보 위치(접합 클러스터, BAC 양 끝) ± window. 이전 버전은 첫 후보 하나만 사용
    REGIONS=$(awk -v w=~{window_bp} '{s=$2-w; if(s<1)s=1; print $1":"s"-"$2+w}' "~{loci_txt}" | sort -u | tr '\n' ' ')
    if [ -z "$REGIONS" ]; then echo "No candidate loci; nothing to extract." >&2; fi
    echo "Extracting reads from: $REGIONS"

    # -F 2304: supplementary(2048) + secondary(256) 제외, -M: 겹치는 region에서 같은 read 중복 출력 방지
    samtools view -b -F 2304 -M "$BAM" $REGIONS \
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
    description: "추출된 HiFi reads를 hybrid reference에 재정렬 (PAF + SAM 동시 출력)"
  }

  parameter_meta {
    hybrid_ref_fa:      { name: "Hybrid reference FASTA" }
    hybrid_ref_fai:     { name: "Hybrid reference FAI" }
    region_reads_fq:    { name: "Region reads FASTQ.gz" }
    out_prefix:         { name: "Output prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    hybrid_paf:         { name: "PAF: region reads vs hybrid reference (detect_chimeric_reads 입력용)" }
    hybrid_sam:         { name: "SAM: region reads vs hybrid reference (BAM 변환 → IGV 시각화용)" }
  }

  input {
    File   hybrid_ref_fa
    File   hybrid_ref_fai
    File   region_reads_fq
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 16
  Int mem_gb    = 32
  # SAM이 PAF보다 크므로 디스크 여유 증가; whole-genome index(~8 GB) 포함
  Int disk_size = ceil(size(region_reads_fq, "GB") * 15 + size(hybrid_ref_fa, "GB") * 6 + 5)

  command <<<
    set -euo pipefail
    minimap2 --version
    # whole-genome hybrid ref → index 한 번만 생성
    minimap2 -x map-hifi -t ~{threads} -d hybrid.mmi ~{hybrid_ref_fa}

    # PAF 출력 — detect_chimeric_reads 파서에서 직접 사용
    minimap2 \
      -cx map-hifi \
      --cs \
      --secondary=no \
      -t ~{threads} \
      hybrid.mmi \
      ~{region_reads_fq} \
      > ~{out_prefix}.hybrid_aligned.paf
    echo "PAF alignments: $(wc -l < ~{out_prefix}.hybrid_aligned.paf)"

    # SAM 출력 — sort_index_bam 태스크에서 sorted BAM+BAI로 변환 (IGV 용)
    minimap2 \
      -ax map-hifi \
      --secondary=no \
      -t ~{threads} \
      hybrid.mmi \
      ~{region_reads_fq} \
      > ~{out_prefix}.hybrid_aligned.sam
    echo "SAM records (non-header): $(grep -cv '^@' ~{out_prefix}.hybrid_aligned.sam || true)"
  >>>

  output {
    File hybrid_paf = "~{out_prefix}.hybrid_aligned.paf"
    File hybrid_sam = "~{out_prefix}.hybrid_aligned.sam"
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

# ── Step 5-3b ─────────────────────────────────────────────────────────────────

task sort_index_bam {
  meta {
    description: "Hybrid reference 정렬 SAM → 정렬된 BAM + BAI 생성 (IGV 시각화용)"
  }

  parameter_meta {
    in_sam:             { name: "Input SAM file (from align_to_hybrid_ref)" }
    out_prefix:         { name: "Output prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
    bam:                { name: "Sorted BAM: region reads vs hybrid reference" }
    bai:                { name: "BAM index (.bai)" }
  }

  input {
    File   in_sam
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 4
  Int mem_gb    = 8
  Int disk_size = ceil(size(in_sam, "GB") * 3 + 5)

  command <<<
    set -euo pipefail
    samtools sort -@ ~{threads} -o ~{out_prefix}.hybrid_aligned.bam ~{in_sam}
    samtools index ~{out_prefix}.hybrid_aligned.bam
    echo "BAM total records: $(samtools view -c ~{out_prefix}.hybrid_aligned.bam)"
  >>>

  output {
    File bam = "~{out_prefix}.hybrid_aligned.bam"
    File bai = "~{out_prefix}.hybrid_aligned.bam.bai"
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
    breakpoint_pos:     { name: "Left junction (last host base before the insert, 1-based); right junction if no left" }
    breakpoint_range:   { name: "chr:left_junction-right_junction (strand-aware read clusters)" }
    chimeric_count:     { name: "Number of chimeric reads supporting the breakpoint" }
    chimeric_tsv:       { name: "Per-read chimeric alignment details (TSV)" }
    tg_strand:          { name: "TG strand relative to reference genome (+/-/N/A)" }
    bac_region:         { name: "[BAC] chr:start-end of the BAC insert; alignments inside it count as construct" }
    junctions_tsv:      { name: "Every junction cluster (side, chrom, pos, reads), all chromosomes" }
  }

  input {
    File   hybrid_paf
    String bac_region = ""
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
# chimeric = 동일 read가 host AND construct 모두에 정렬된 경우.
# construct = transgene 서열, BAC 모드에서는 BAC 구간(bac_region)에 정렬된 부분도 포함
# (BAC 안의 Cre-Aldh1l1 같은 접합은 construct 내부 접합이므로 integration이 아님)
BAC = "~{bac_region}".replace(",", "")
bac = None
if BAC:
    _c, _se = BAC.split(":")
    bac = (_c, int(_se.split("-")[0]), int(_se.split("-")[1]))

def is_construct(a):
    if a["target"] == tg_seqname:
        return True
    return bac is not None and a["target"] == bac[0] and min(a["t_end"], bac[2]) - max(a["t_start"], bac[1]) > 0.5 * (a["t_end"] - a["t_start"])

chimeric = []
for qname, alns in reads.items():
    host = [a for a in alns if not is_construct(a)]
    cons = [a for a in alns if is_construct(a)]
    if host and cons:
        chimeric.append({"name": qname, "host_alns": host, "tg_alns": cons})

print(f"Total chimeric reads: {len(chimeric)}")

# ── Breakpoint 결정 ───────────────────────────────────────────────────────────
# Read 내 좌표에서 host와 TG 정렬이 맞닿은 지점을 host reference 좌표로 환산.
# host 정렬의 strand를 반드시 반영해야 함 ('-' 정렬은 read 끝 = ref t_start).
#   left  junction: 삽입 바로 앞의 마지막 host 염기 (삽입은 ref 좌표가 큰 쪽)
#   right junction: 삽입 바로 뒤의 첫 host 염기
# 좌표는 1-based. left/right를 따로 모아 각각의 최빈 클러스터를 접합부로 보고,
# 두 접합부 사이(right - left - 1)를 host 결실 길이로 보고.
TOLERANCE  = 1000   # read 상에서 host/TG 정렬 끝이 이만큼 떨어져도 인접으로 판정
CLUSTER_BP = 20     # 같은 접합부로 묶는 ref 좌표 허용 오차

junction_calls = []   # (side, chrom, pos)
host_chrs      = []
tsv_rows       = ["read_name\thost_chr\thost_strand\thost_q_start\thost_q_end\thost_t_start\thost_t_end\ttg_q_start\ttg_q_end\tbreakpoint_pos\tjunction_side\tread_order"]

for c in chimeric:
    for ha in c["host_alns"]:
        for ta in c["tg_alns"]:
            if abs(ha["q_end"] - ta["q_start"]) <= TOLERANCE:          # read: host → TG
                order = "host->TG"
                if ha["strand"] == "+":
                    bp, side = ha["t_end"], "left"
                else:
                    bp, side = ha["t_start"] + 1, "right"
            elif abs(ta["q_end"] - ha["q_start"]) <= TOLERANCE:        # read: TG → host
                order = "TG->host"
                if ha["strand"] == "+":
                    bp, side = ha["t_start"] + 1, "right"
                else:
                    bp, side = ha["t_end"], "left"
            else:
                continue
            junction_calls.append((side, ha["target"], bp))
            host_chrs.append(ha["target"])
            tsv_rows.append(f'{c["name"]}\t{ha["target"]}\t{ha["strand"]}\t{ha["q_start"]}\t{ha["q_end"]}\t{ha["t_start"]}\t{ha["t_end"]}\t{ta["q_start"]}\t{ta["q_end"]}\t{bp}\t{side}\t{order}')

bp_chr = max(set(host_chrs), key=host_chrs.count) if host_chrs else "unknown"

def best_cluster(positions):
    """Largest group of positions within CLUSTER_BP; returns (mode position, size)."""
    positions = sorted(positions)
    best, i = (None, 0), 0
    for k in range(len(positions)):
        while positions[k] - positions[i] > CLUSTER_BP:
            i += 1
        if k - i + 1 > best[1]:
            grp = positions[i:k + 1]
            best = (max(set(grp), key=grp.count), k - i + 1)
    return best

left_pos,  n_left  = best_cluster([p for sd, ch, p in junction_calls if sd == "left"  and ch == bp_chr])
right_pos, n_right = best_cluster([p for sd, ch, p in junction_calls if sd == "right" and ch == bp_chr])

# every junction cluster on every chromosome (>= 2 reads), for multi-site / rearranged integrations
all_clusters = []
for side in ("left", "right"):
    for chrom in sorted({ch for sd, ch, p in junction_calls if sd == side}):
        pos = sorted(p for sd, ch, p in junction_calls if sd == side and ch == chrom)
        grp = [pos[0]]
        for p in pos[1:] + [None]:
            if p is not None and p - grp[-1] <= CLUSTER_BP:
                grp.append(p); continue
            if len(grp) >= 2:
                all_clusters.append((side, chrom, max(set(grp), key=grp.count), len(grp)))
            if p is not None:
                grp = [p]
all_clusters.sort(key=lambda x: -x[3])
with open(f"{out_prefix}.junctions.tsv", "w") as f:
    f.write("side\tchrom\tpos\treads\n")
    for side, chrom, p, n in all_clusters:
        f.write(f"{side}\t{chrom}\t{p}\t{n}\n")
print("Junction clusters (>= 2 reads):")
for side, chrom, p, n in all_clusters:
    print(f"  {side:5s} {chrom}:{p}  {n} reads")
other = [c for c in all_clusters if c[1] != bp_chr]

if left_pos or right_pos:
    bp_median  = left_pos or right_pos
    bp_support = f"{n_left} left-junction / {n_right} right-junction reads"
    if left_pos and right_pos:
        bp_range = f"{bp_chr}:{left_pos:,}-{right_pos:,}"
        gap = right_pos - left_pos - 1
        if abs(gap) > 1000000:
            bp_support += f" (junctions {abs(gap):,} bp apart: rearrangement at the integration, not a simple host deletion)"
        elif gap >= 0:
            bp_support += f" (host deletion {gap:,} bp between junctions)"
        else:
            bp_support += f" (junctions overlap by {-gap:,} bp: target-site duplication or rearrangement)"
    else:
        bp_range = f"{bp_chr}:{bp_median:,} ({'left' if left_pos else 'right'} junction only)"
    if other:
        bp_support += f"; +{len(other)} junction cluster(s) on other chromosomes: " + ", ".join(f"{s} {c}:{p} ({n})" for s, c, p, n in other[:4])
    print(f"Breakpoint chromosome: {bp_chr}")
    print(f"Left junction:  {bp_chr}:{left_pos if left_pos else 'not found'} ({n_left} reads)")
    print(f"Right junction: {bp_chr}:{right_pos if right_pos else 'not found'} ({n_right} reads)")
    print(f"Breakpoint range: {bp_range}")
    print(f"Total chimeric reads:  {len(chimeric)}")
else:
    bp_median  = 0
    bp_range   = "N/A"
    bp_support = "N/A"
    print("WARNING: No clear breakpoint detected from chimeric reads.")

with open(f"{out_prefix}.chimeric_reads.tsv", "w") as f:
    f.write("\n".join(tsv_rows) + "\n")
with open("breakpoint_chr.txt",   "w") as f: f.write(bp_chr)
with open("breakpoint_pos.txt",   "w") as f: f.write(str(bp_median))
with open("breakpoint_range.txt", "w") as f: f.write(bp_range)
with open("chimeric_count.txt",   "w") as f: f.write(str(len(chimeric)))
with open("bp_support.txt",       "w") as f: f.write(bp_support)

# ── TG strand vs reference ────────────────────────────────────────────────
# host_strand == tg_strand → 동일 방향 read 정렬 → TG는 ref 기준 정방향(+)
# host_strand != tg_strand → 반대 방향 → TG는 ref 기준 역방향(-)
tg_strand_votes = []
for c in chimeric:
    for ha in c["host_alns"]:
        for ta in c["tg_alns"]:
            if (abs(ha["q_end"] - ta["q_start"]) <= TOLERANCE or
                abs(ta["q_end"] - ha["q_start"]) <= TOLERANCE):
                tg_strand_votes.append("+" if ha["strand"] == ta["strand"] else "-")

if tg_strand_votes:
    plus  = tg_strand_votes.count("+")
    minus = tg_strand_votes.count("-")
    tg_strand = "+" if plus >= minus else "-"
    print(f"TG strand vs reference: {tg_strand} ({plus}+/{minus}- from {len(tg_strand_votes)} junctions)")
else:
    tg_strand = "N/A"
    print("TG strand vs reference: N/A (no chimeric junctions)")

with open("tg_strand.txt", "w") as f: f.write(tg_strand)
PYEOF
  >>>

  output {
    File   chimeric_tsv     = "~{out_prefix}.chimeric_reads.tsv"
    File   junctions_tsv    = "~{out_prefix}.junctions.tsv"
    String breakpoint_chr   = read_string("breakpoint_chr.txt")
    String breakpoint_pos   = read_string("breakpoint_pos.txt")
    String breakpoint_range = read_string("breakpoint_range.txt")
    String chimeric_count   = read_string("chimeric_count.txt")
    String bp_support       = read_string("bp_support.txt")
    String tg_strand        = read_string("tg_strand.txt")
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

# ── Step 5-4b (BAC) ───────────────────────────────────────────────────────────

task bac_depth_qc {
  meta {
    description: "BAC 모드: BAC 구간 depth / 양쪽 flank depth. 1.3배 이상이면 BAC 사본이 native locus에 겹쳐 정렬된 것"
  }

  parameter_meta {
    haplotagged_bam:    { name: "HiFi reads aligned to reference" }
    haplotagged_bai:    { name: "BAM index" }
    bac_region:         { name: "chr:start-end of the BAC insert" }
    flank_bp:           { name: "Flank length on each side (default: 200000)" }
    runtime_attributes: { name: "Runtime attribute structure" }
    qc:                 { name: "QC string for the report" }
  }

  input {
    File   haplotagged_bam
    File   haplotagged_bai
    String bac_region
    Int    flank_bp = 200000
    RuntimeAttributes runtime_attributes
  }

  command <<<
    set -euo pipefail
    BAM="~{haplotagged_bam}"
    ln -sf "~{haplotagged_bai}" "${BAM}.bai" 2>/dev/null || true
    R=$(echo "~{bac_region}" | tr -d ',')
    C=${R%%:*}; S=$(echo ${R#*:} | cut -d- -f1); E=$(echo ${R#*:} | cut -d- -f2)
    mean() { samtools depth -a -Q 20 -r "$1" "$BAM" | awk '{s+=$3; n++} END{printf "%.1f", (n ? s/n : 0)}'; }
    IN=$(mean "$C:$S-$E")
    L=$(mean "$C:$(( S > ~{flank_bp} ? S - ~{flank_bp} : 1 ))-$(( S - 1 ))")
    RR=$(mean "$C:$(( E + 1 ))-$(( E + ~{flank_bp} ))")
    awk -v i="$IN" -v l="$L" -v r="$RR" -v reg="$R" 'BEGIN{
      f=(l+r)/2; ratio=(f>0 ? i/f : 0);
      msg=sprintf("BAC_DEPTH %s mean %.1fx vs flanks %.1fx/%.1fx (ratio %.2f)", reg, i, l, r, ratio);
      if (ratio >= 1.3) msg=msg " -> extra BAC copies map onto the native locus; the integration site is elsewhere (see junction clusters on other chromosomes)";
      print msg > "qc.txt"; print msg }'
  >>>

  output { String qc = read_string("qc.txt") }

  runtime {
    docker:                "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
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
    chimeric_count:         { name: "Total number of chimeric reads detected" }
    bp_support:             { name: "Breakpoint support string (e.g. '5 host->TG reads')" }
    tg_strand_vs_ref:       { name: "TG strand relative to reference genome from hybrid ref chimeric reads (+/-/N/A)" }
    nearest_gene_name:      { name: "Gene at/near the breakpoint (from annotation; 'N/A' if not run)" }
    nearest_gene_type:      { name: "Biotype of nearest gene" }
    insertion_feature_type: { name: "exon / intron / intergenic" }
    nearest_exon_dist_bp:   { name: "Distance to nearest exon boundary (bp)" }
    qc_warnings:            { name: "Semicolon-separated QC warnings from extract_integration_coords ('none' if clean)" }
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
    String bp_support             = "N/A"
    String tg_strand_vs_ref       = "N/A"

    # 유전자 어노테이션 결과 (선택적 — 미실행 시 기본값)
    String nearest_gene_name      = "N/A"
    String nearest_gene_type      = "N/A"
    String insertion_feature_type = "N/A"
    String nearest_exon_dist_bp   = "N/A"

    # extract_integration_coords 경고 (다중 위치, 염색체 전환, host 유래 construct 서열 등)
    String qc_warnings            = "none"

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
bp_support     = "~{bp_support}"
has_hybrid     = prec_chr != "N/A" and prec_pos > 0
tg_strand_ref  = "~{tg_strand_vs_ref}"   # hybrid ref chimeric reads 기반 TG 절대 방향
has_tg_strand  = tg_strand_ref not in ("N/A", "?")

# 어노테이션 결과
gene_name      = "~{nearest_gene_name}"
gene_type      = "~{nearest_gene_type}"
feat_type      = "~{insertion_feature_type}"
exon_dist_str  = "~{nearest_exon_dist_bp}"
exon_dist      = int(exon_dist_str) if exon_dist_str.lstrip("-").isdigit() else -1
has_annot      = gene_name != "N/A"

# QC 경고 — 하나라도 있으면 리포트 상단에 표시 (단일 위치 결과를 그대로 믿지 말 것)
import html as _html
qc_raw         = r'''~{qc_warnings}'''.strip()
qc_warns       = [w.strip() for w in qc_raw.split(";") if w.strip() and w.strip() not in ("none", "N/A")]
warn_html      = ("" if not qc_warns else
                  '<div style="background:#fff4e5;border-left:4px solid #e67e22;border-radius:8px;padding:12px 16px;margin-bottom:14px">'
                  '<div style="font-weight:700;color:#a04000;margin-bottom:6px">QC warnings &#8212; review before reporting a single integration site</div><ul style="margin:0;padding-left:18px">'
                  + "".join(f"<li style=\"font-size:0.86em;margin:2px 0\">{_html.escape(w)}</li>" for w in qc_warns) + "</ul></div>\n")

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
            # de:f: tag = divergence → identity = (1 - de) * 100
            de = 0.0
            for tag in p[12:]:
                if tag.startswith("de:f:"):
                    try: de = float(tag.split(":")[2])
                    except: pass
                    break
            q_start, q_end = int(p[2]), int(p[3])
            t_start, t_end = int(p[7]), int(p[8])
            alns[p[0]].append({
                "q_len": int(p[1]), "q_start": q_start, "q_end": q_end,
                "chrom": p[5], "t_start": t_start, "t_end": t_end,
                "strand": p[4], "matches": int(p[9]), "mapq": int(p[11]),
                "identity": round((1.0 - de) * 100, 2),
                "aln_len":  q_end - q_start,
                "ref_span": t_end - t_start,
            })
    return alns

def estimate_ref_pos(tg_center, ref_alns):
    for a in sorted(ref_alns, key=lambda x: x["q_start"]):
        if a["q_start"] <= tg_center <= a["q_end"]:
            ratio = (tg_center - a["q_start"]) / max(1, a["q_end"] - a["q_start"])
            return a["chrom"], int(a["t_start"] + ratio * (a["t_end"] - a["t_start"])), a["mapq"], a["strand"]
    return "unknown", 0, 0, "?"

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
            chrom, ref_pos, ref_mapq, contig_ref_strand = estimate_ref_pos(center, ref_aln)
            tg_cov = (cp["tg_end"] - cp["tg_start"]) / cp["tg_len"] * 100
            # ── 교정된 TG strand (genome 기준) ──────────────────────────────
            # 우선순위 1: mapq>0 contig → contig 방향 이용해 TG 절대 방향 계산
            # 우선순위 2: hybrid ref chimeric reads 결과
            # 우선순위 3: 불명확
            if ref_mapq > 0 and contig_ref_strand != "?":
                corrected_strand = "+" if contig_ref_strand == cp["strand"] else "-"
            elif has_tg_strand:
                corrected_strand = tg_strand_ref
            else:
                corrected_strand = "?"
            enriched.append({**cp, "chrom": chrom, "ref_pos": ref_pos,
                              "ref_mapq": ref_mapq, "tg_cov": tg_cov,
                              "corrected_strand": corrected_strand,
                              "contig_ref_strand": contig_ref_strand})
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

if qc_warns:
    lines.append("\n[QC WARNINGS — review before reporting a single integration site]")
    for w in qc_warns:
        lines.append(f"  ! {w}")

if has_hybrid:
    lines.append(f"\n[PRECISE BREAKPOINT — hybrid reference analysis]")
    lines.append(f"  Chromosome: {prec_chr}")
    lines.append(f"  Position:   {prec_pos:,} bp")
    lines.append(f"  Support:    {bp_support} / {chimeric_n} total chimeric reads")
    lines.append(f"  Range:      {bp_range}")
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
            cs = cp["corrected_strand"]
            cs_lbl = {"+" : "정방향", "-": "역방향"}.get(cs, "불명확")
            raw_s  = cp["strand"]
            raw_lbl = "역방향" if raw_s == "-" else "정방향"
            strand_note = "" if cs == raw_s else f" (de novo contig 기준: {raw_s} {raw_lbl})"
            lines.append(f"      Direction:    {cs} ({cs_lbl}){strand_note}  ← genome 기준")
            lines.append(f"      Ref position: {cp['chrom']}:{cp['ref_pos']:,} (mapq={cp['ref_mapq']})")
            tsv_rows.append("\t".join([
                sample_id, transgene_name, hap, c["name"],
                str(c["contig_len"]), str(copy_num), str(i),
                str(cp["c_start"]), str(cp["c_end"]),
                f"{cp['tg_cov']:.1f}", cp["corrected_strand"],
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

def strand_badge(s, source=""):
    if s == "-":   c, lbl = "#e74c3c", "&#8722; 역방향"
    elif s == "+": c, lbl = "#27ae60", "+ 정방향"
    else:          c, lbl = "#999",    "? 불명확"
    badge = f'<span style="background:{c};color:#fff;padding:1px 7px;border-radius:10px;font-size:0.79em">{lbl}</span>'
    if source:
        badge += f' <span style="font-size:0.74em;color:#888">{source}</span>'
    return badge

def bp_fmt(n):
    if n >= 1_000_000: return f"{n/1_000_000:.2f} Mb"
    if n >= 1_000:     return f"{n/1_000:.1f} kb"
    return f"{n} bp"

def collinearity_class(alns):
    """블록 순서 분석: 각 블록에 색상 힌트 부여."""
    if not alns: return []
    sorted_alns = sorted(enumerate(alns), key=lambda x: x[1]["q_start"])
    hints = [""] * len(alns)
    prev_chrom, prev_t_end, prev_strand = None, None, None
    for rank, (orig_idx, a) in enumerate(sorted_alns):
        chrom = a["chrom"]; strand = a["strand"]
        if prev_chrom is None:
            hints[orig_idx] = "syntenic"
        elif chrom != prev_chrom:
            hints[orig_idx] = "transloc"
        elif strand != prev_strand:
            hints[orig_idx] = "inversion"
        elif (strand == "+" and a["t_start"] < prev_t_end) or \
             (strand == "-" and a["t_end"]   > prev_t_end):
            hints[orig_idx] = "overlap"
        else:
            hints[orig_idx] = "syntenic"
        prev_chrom = chrom
        prev_t_end = a["t_end"] if strand == "+" else a["t_start"]
        prev_strand = strand
    return hints

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
        hints = collinearity_class(ref_alns)
        hint_bg = {"syntenic":"#f0fff0","inversion":"#fff0f0","transloc":"#fffbe6","overlap":"#f0f0ff","":""}
        hint_icon = {"syntenic":"","inversion":"&#8635; ","transloc":"&#8646; ","overlap":"&#8651; ","":""}
        hint_tip  = {"syntenic":"Syntenic (공선성)","inversion":"Inversion (역위)","transloc":"Translocation (전좌)","overlap":"Overlap/Inversion (중복/역위)","":""}
        copy_rows=""
        for i,cp in enumerate(copies,1):
            cs = cp["corrected_strand"]
            raw_s = cp["strand"]
            src = ""
            if cp["ref_mapq"] > 0 and cp["contig_ref_strand"] != "?":
                src = "contig→ref"
            elif has_tg_strand:
                src = "hybrid ref"
            mismatch_warn = ""
            if cs != "?" and cs != raw_s:
                mismatch_warn = f' <span title="de novo contig 방향({raw_s})과 다름" style="color:#e67e22;font-size:0.75em">&#9888;</span>'
            copy_rows+=(
                f'<tr><td style="text-align:center;font-weight:600">{i}</td>'
                f'<td style="font-family:monospace;white-space:nowrap">{cp["c_start"]:,}&#8211;{cp["c_end"]:,}</td>'
                f'<td>{pct_bar(cp["tg_cov"])}</td>'
                f'<td>{strand_badge(cs, src)}{mismatch_warn}</td>'
                f'<td style="font-family:monospace">{cp["chrom"]}:{cp["ref_pos"]:,}</td>'
                f'<td style="text-align:center">{cp["ref_mapq"]}</td></tr>'
            )
            if i<=len(gaps):
                copy_rows+=f'<tr><td colspan="6" style="text-align:center;color:#888;font-size:0.82em;padding:1px">&#8595; gap {gaps[i-1]:,} bp</td></tr>'
        ref_rows=""
        for idx,a in enumerate(ref_alns):
            h = hints[idx] if idx < len(hints) else ""
            bg = hint_bg.get(h,"")
            icon = hint_icon.get(h,"")
            tip  = hint_tip.get(h,"")
            ref_span_fmt = bp_fmt(a["ref_span"])
            aln_len_fmt  = bp_fmt(a["aln_len"])
            ident_str    = f'{a["identity"]:.2f}%' if "identity" in a else "N/A"
            ref_rows+=(
                f'<tr style="background:{bg}" title="{tip}">'
                f'<td style="font-family:monospace;white-space:nowrap;font-size:0.82em">{icon}{a["chrom"]}:{a["t_start"]:,}&#8211;{a["t_end"]:,}</td>'
                f'<td style="font-family:monospace;white-space:nowrap;font-size:0.82em">{a["q_start"]:,}&#8211;{a["q_end"]:,}</td>'
                f'<td style="text-align:center">{a["strand"]}</td>'
                f'<td style="text-align:right;font-family:monospace;font-size:0.82em">{ref_span_fmt}</td>'
                f'<td style="text-align:right;font-family:monospace;font-size:0.82em">{aln_len_fmt}</td>'
                f'<td style="text-align:right;font-size:0.82em">{ident_str}</td>'
                f'<td style="text-align:center">{a["mapq"]}</td></tr>'
            )
        n=len(ref_alns)
        n_inversion = hints.count("inversion")
        n_transloc  = hints.count("transloc")
        if n <= 1:
            ref_note = '<p style="color:#27ae60;font-size:0.82em;margin-top:5px">&#10003; 단일 블록 &#8212; 단순 삽입</p>'
        else:
            warns = []
            if n_inversion: warns.append(f'Inversion (역위) {n_inversion}개')
            if n_transloc:  warns.append(f'Translocation (전좌) {n_transloc}개')
            warn_txt = " / ".join(warns) if warns else ""
            extra = f" ({warn_txt})" if warn_txt else ""
            ref_note = f'<p style="color:#e67e22;font-size:0.82em;margin-top:5px">&#9888; {n}개 블록{extra} &#8594; 추가 구조변이 가능</p>'
        # 범례
        legend = ('<div style="font-size:0.74em;color:#888;margin-top:4px">'
                  '<span style="background:#f0fff0;padding:1px 5px;border-radius:3px">Syntenic (공선성)</span> '
                  '<span style="background:#fff0f0;padding:1px 5px;border-radius:3px">&#8635; Inversion (역위)</span> '
                  '<span style="background:#fffbe6;padding:1px 5px;border-radius:3px">&#8646; Translocation (전좌)</span>'
                  '</div>')
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
            f'<th style="padding:3px 6px;text-align:left;border-bottom:2px solid #dee2e6">방향 (genome 기준)</th>'
            f'<th style="padding:3px 6px;text-align:left;border-bottom:2px solid #dee2e6">Ref 위치 (추정)</th>'
            f'<th style="padding:3px 6px;border-bottom:2px solid #dee2e6">mapq</th>'
            f'</tr></thead><tbody>{copy_rows}</tbody></table></div>'
            f'<div><div style="font-size:0.76em;color:#666;font-weight:600;margin-bottom:3px">Reference 정렬 블록</div>'
            f'<table style="width:100%;border-collapse:collapse;font-size:0.84em"><thead><tr style="background:#f8f9fa">'
            f'<th style="padding:3px 6px;text-align:left;border-bottom:2px solid #dee2e6">Ref 좌표</th>'
            f'<th style="padding:3px 6px;text-align:left;border-bottom:2px solid #dee2e6">Contig 좌표</th>'
            f'<th style="padding:3px 6px;border-bottom:2px solid #dee2e6">Strand</th>'
            f'<th style="padding:3px 6px;text-align:right;border-bottom:2px solid #dee2e6">Ref 길이</th>'
            f'<th style="padding:3px 6px;text-align:right;border-bottom:2px solid #dee2e6">Aln 길이</th>'
            f'<th style="padding:3px 6px;text-align:right;border-bottom:2px solid #dee2e6">Identity</th>'
            f'<th style="padding:3px 6px;border-bottom:2px solid #dee2e6">mapq</th>'
            f'</tr></thead><tbody>{ref_rows}</tbody></table>{ref_note}{legend}</div>'
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
        ("방향 (genome)",
                       hap_val("hap1",lambda c:("역방향 (&#8722;)" if c["copies"][0]["corrected_strand"]=="-" else "정방향 (+)" if c["copies"][0]["corrected_strand"]=="+" else "불명확 (?)") if c["copies"] else "N/A"),
                       hap_val("hap2",lambda c:("역방향 (&#8722;)" if c["copies"][0]["corrected_strand"]=="-" else "정방향 (+)" if c["copies"][0]["corrected_strand"]=="+" else "불명확 (?)") if c["copies"] else "N/A")),
        ("Ref 정렬 블록",ref_badge("hap1"),ref_badge("hap2")),
    ]
)

# 요약 카드 데이터
all_copies  = [cp for hap in ["hap1","hap2"] for c in hap_data.get(hap,[]) for cp in c["copies"]]
max_copies  = max((len(c["copies"]) for hap in ["hap1","hap2"] for c in hap_data.get(hap,[])), default=0)
# 방향: hybrid ref 결과 우선, 없으면 교정된 de novo 다수결
if has_tg_strand:
    rep_strand = tg_strand_ref
else:
    corr_strands = [cp["corrected_strand"] for cp in all_copies if cp["corrected_strand"] != "?"]
    rep_strand   = max(set(corr_strands), key=corr_strands.count) if corr_strands else "?"
strand_lbl  = {"+" : "정방향", "-": "역방향"}.get(rep_strand, "불명확")
strand_src  = "hybrid ref" if has_tg_strand else "de novo (교정)"

# 표시할 대표 위치
if has_hybrid:
    disp_chr   = prec_chr
    disp_pos   = f"{prec_pos:,} bp"
    disp_prec  = bp_support
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
        f'<tr><td style="color:#666;padding:4px 6px">Support</td><td style="padding:4px 6px">{bp_support} / {chimeric_n}개 total</td></tr>'
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
        f'<td style="padding:4px 7px"><span style="background:#d4edda;color:#155724;padding:1px 6px;border-radius:8px;font-size:0.8em">{bp_support}</span></td></tr>'
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
    f'{warn_html}'

    f'<div style="display:grid;grid-template-columns:repeat({"5" if has_annot else "4"},1fr);gap:9px;margin-bottom:14px">\n'
    f'  <div style="background:#fff;border-radius:8px;padding:13px;border-left:4px solid #1565c0"><div style="font-size:0.72em;color:#888;font-weight:600;text-transform:uppercase">통합 염색체</div><div style="font-size:1.35em;font-weight:700;color:#1565c0;margin-top:2px">{disp_chr}</div><div style="font-size:0.74em;color:#666;margin-top:1px">reference</div></div>\n'
    f'  <div style="background:#fff;border-radius:8px;padding:13px;border-left:4px solid {disp_color}"><div style="font-size:0.72em;color:#888;font-weight:600;text-transform:uppercase">통합 위치</div><div style="font-size:1.2em;font-weight:700;color:{disp_color};margin-top:2px;font-family:monospace">{pos_label}</div><div style="font-size:0.74em;color:#666;margin-top:1px">{disp_prec}</div></div>\n'
    f'  <div style="background:#fff;border-radius:8px;padding:13px;border-left:4px solid #e74c3c"><div style="font-size:0.72em;color:#888;font-weight:600;text-transform:uppercase">탠덤 카피 수</div><div style="font-size:1.35em;font-weight:700;color:#e74c3c;margin-top:2px">{max_copies}카피</div><div style="font-size:0.74em;color:#666;margin-top:1px">hap1 &amp; hap2</div></div>\n'
    f'  <div style="background:#fff;border-radius:8px;padding:13px;border-left:4px solid #6a1b9a"><div style="font-size:0.72em;color:#888;font-weight:600;text-transform:uppercase">삽입 방향 (genome)</div><div style="font-size:1.15em;font-weight:700;color:#6a1b9a;margin-top:2px">{rep_strand} {strand_lbl}</div><div style="font-size:0.74em;color:#666;margin-top:1px">{strand_src}</div></div>\n'
    + (f'  <div style="background:#fff;border-radius:8px;padding:13px;border-left:4px solid #2e7d32"><div style="font-size:0.72em;color:#888;font-weight:600;text-transform:uppercase">삽입 유전자</div><div style="font-size:1.2em;font-weight:700;color:#2e7d32;margin-top:2px"><em>{gene_name}</em></div><div style="font-size:0.74em;color:#666;margin-top:1px">{feat_type}</div></div>\n' if has_annot else '')
    + '</div>\n\n'

    + hybrid_section
    + annot_section

    + '<div style="background:#fff;border-radius:8px;padding:16px;margin-bottom:14px;border:1px solid #dee2e6">\n'
    '<h2 style="margin-top:0">&#129516; De Novo Assembly 분석</h2>\n'
    '<div style="display:grid;grid-template-columns:auto 1fr;gap:14px;align-items:start;margin-bottom:14px">\n'
    f'<table style="border-collapse:collapse;font-size:0.85em;white-space:nowrap"><thead><tr style="background:#f8f9fa"><th style="padding:4px 7px;text-align:left;border-bottom:2px solid #dee2e6">항목</th><th style="padding:4px 7px;text-align:center;border-bottom:2px solid #dee2e6">hap1</th><th style="padding:4px 7px;text-align:center;border-bottom:2px solid #dee2e6">hap2</th></tr></thead><tbody>{cmp_rows}</tbody></table>\n'
    '<div style="background:#fff8e1;border-radius:6px;padding:11px;font-size:0.84em"><div style="font-weight:600;color:#f57f17;margin-bottom:5px">&#9888; 주의사항</div><ul style="list-style:disc;padding-left:13px;color:#555;line-height:1.8"><li>De novo ref 위치는 선형 보간 추정값 &#8212; hybrid ref breakpoint가 더 정확</li><li>방향(genome 기준): mapq&gt;0 contig는 ref 정렬로 교정, mapq=0은 hybrid ref chimeric reads 사용</li><li>mapq 0/1 : 다중 카피로 인한 정상적 낮은 mapq &#8212; 방향은 별도 소스로 교정됨</li><li>복수 블록 &#8594; Syntenic(공선성)/Inversion(역위)/Translocation(전좌) 여부는 Reference 정렬 블록 색상 참조</li></ul></div>\n'
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
