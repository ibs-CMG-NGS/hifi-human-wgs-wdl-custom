version 1.0

import "../structs.wdl"

# ── Task 0: HiFi kinetics → by-strand ip/pw 변환 (ipdSummary 입력용) ──────────
task ccs_kinetics_bystrandify {
  meta { description: "HiFi CCS kinetics(fi/fp/ri/rp)를 strand별 ip/pw 태그 pseudo-subread BAM으로 변환 — ipdSummary가 읽을 수 있는 형식. 정렬 전 raw CCS BAM에 적용해야 함" }

  parameter_meta {
    raw_hifi_bam:       { name: "Unaligned raw HiFi BAM with fi/fp/ri/rp kinetics" }
    sample_id:          { name: "Sample identifier" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    File   raw_hifi_bam
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 8
  Int mem_gb    = 16
  Int disk_size = ceil(size(raw_hifi_bam, "GB") * 4 + 20)

  command <<<
    set -euo pipefail
    ccs-kinetics-bystrandify --version

    # fi/fp/ri/rp → strand별 ip/pw (각 CCS read가 forward/reverse 2개 pseudo-read로 분리)
    ccs-kinetics-bystrandify \
      --num-threads ~{threads} \
      "~{raw_hifi_bam}" \
      "~{sample_id}.bystrand.bam"
  >>>

  output {
    File   bystrand_bam = "~{sample_id}.bystrand.bam"
  }

  runtime {
    docker:                "quay.io/biocontainers/pbtk:3.5.0--h9ee0642_0"
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

# ── Task 1: Kinetics 보존 정렬 ────────────────────────────────────────────────
task pbmm2_align_kinetics {
  meta { description: "ip/pw kinetics + MM/ML base-mod tags 모두 보존하며 GRCm39에 정렬" }

  parameter_meta {
    raw_hifi_bam:       { name: "Raw HiFi BAM (kinetics ip/pw tags 포함)" }
    ref_fasta:          { name: "Reference genome FASTA (GRCm39)" }
    ref_fai:            { name: "FAI index for ref_fasta" }
    sample_id:          { name: "Sample ID (read group SM tag에 사용)" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    File   raw_hifi_bam
    File   ref_fasta
    File   ref_fai
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 32
  Int mem_gb    = 80
  Int disk_size = ceil(size(raw_hifi_bam, "GB") * 4 + size(ref_fasta, "GB") * 2 + 20)

  command <<<
    set -euo pipefail
    pbmm2 --version

    # strip_kinetics 플래그 제거 → ip/pw 태그 유지
    # BAM 입력은 이미 @RG 헤더를 갖고 있어 --rg로 덮어쓸 수 없음 → --sample만 지정
    # 기본 --sort-memory(768M)*--sort-threads(4)=3G 버퍼는 692개 spill chunk를 만들어
    # merge 단계가 디스크 I/O에 묶이는 원인이었음 → 8 threads * 4G = 32G 버퍼로 chunk 수 최소화
    pbmm2 align \
      --preset HIFI \
      --sort \
      --sort-memory 4G \
      --sort-threads 8 \
      --sample "~{sample_id}" \
      --num-threads ~{threads} \
      --log-level INFO \
      "~{ref_fasta}" \
      "~{raw_hifi_bam}" \
      "~{sample_id}.aligned.bam"

    # pbmm2 --sort가 이미 BAI를 생성함 (samtools index 재실행은 323GB 파일을 다시
    # 읽는 불필요한 중복 단계였음)
    echo "aligned_reads=$(samtools idxstats "~{sample_id}.aligned.bam" | awk '{sum+=$3+$4} END{print sum}')"
    ip_count=$(samtools view "~{sample_id}.aligned.bam" | head -1000 | awk '{for(i=12;i<=NF;i++) if($i~/^ip:/) c++} END{print c+0}' || true)
    echo "ip_tag_check_first1000=${ip_count}"
  >>>

  output {
    File   aligned_bam = "~{sample_id}.aligned.bam"
    File   aligned_bai = "~{sample_id}.aligned.bam.bai"
  }

  runtime {
    docker:                "~{runtime_attributes.container_registry}/pbmm2@sha256:5f3f4d1f5dbea5cd4c388ee26b2fecbbb7dbcef449343633e039dca3d3725859"
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

# ── Task 2a: jasmine ML 변형 calling ──────────────────────────────────────────
task jasmine_modification_calling {
  meta { description: "ip/pw kinetics → ML 모델로 5mC / 5hmC / 6mA 확률 산출 후 MM/ML 태그 업데이트" }

  parameter_meta {
    aligned_bam:        { name: "pbmm2 정렬 BAM (ip/pw 태그 포함)" }
    aligned_bai:        { name: "BAI index" }
    sample_id:          { name: "Sample ID" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    File   aligned_bam
    File   aligned_bai
    String sample_id
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 32
  Int mem_gb    = 48
  Int disk_size = ceil(size(aligned_bam, "GB") * 3 + 10)

  command <<<
    set -euo pipefail
    jasmine --version 2>/dev/null || true

    jasmine \
      --log-level INFO \
      --num-threads ~{threads} \
      "~{aligned_bam}" \
      "~{sample_id}.modified.bam"
  >>>

  output {
    File   modified_bam = "~{sample_id}.modified.bam"
  }

  runtime {
    docker:                "quay.io/biocontainers/pbjasmine:26.1.3--hd63eeec_0"
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

# ── Task 2b: 그룹 BAM 병합 (~90x 커버리지 확보) ──────────────────────────────
task merge_group_bams {
  meta { description: "동일 그룹 3개 replicate BAM 병합 → ~90x pooled BAM" }

  parameter_meta {
    input_bams:         { name: "Array of jasmine-modified BAMs (3 replicates)" }
    out_prefix:         { name: "Output prefix (e.g. control, drug_a)" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    Array[File] input_bams
    String      out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 8
  Int mem_gb    = 16
  Int disk_size = ceil(size(input_bams, "GB") * 3 + 10)

  command <<<
    set -euo pipefail

    BAM_LIST=(~{sep=" " input_bams})
    echo "Input BAMs: ${#BAM_LIST[@]}"

    if [ ${#BAM_LIST[@]} -eq 1 ]; then
      cp "${BAM_LIST[0]}" "~{out_prefix}.merged.bam"
    else
      samtools merge -f -@ ~{threads} \
        "~{out_prefix}.merged.bam" \
        "${BAM_LIST[@]}"
    fi

    samtools sort -@ ~{threads} -o "~{out_prefix}.merged.sorted.bam" "~{out_prefix}.merged.bam"
    mv "~{out_prefix}.merged.sorted.bam" "~{out_prefix}.merged.bam"
    samtools index "~{out_prefix}.merged.bam"

    echo "Total merged reads: $(samtools view -c ~{out_prefix}.merged.bam)"
    samtools coverage "~{out_prefix}.merged.bam" | awk 'NR>1{s+=$6; n++} END{printf "mean_coverage=%.1f\n", s/n}' || true
  >>>

  output {
    File   merged_bam = "~{out_prefix}.merged.bam"
    File   merged_bai = "~{out_prefix}.merged.bam.bai"
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

# ── Task 2c: MM 태그 read 필터 (modkit pileup 전처리) ───────────────────────
task filter_modbam {
  meta { description: "MM 태그 없는 primary read 제거 — modkit 0.6.3은 태그 없는 read를 만나면 worker deadlock에 빠짐" }

  parameter_meta {
    merged_bam:         { name: "Pooled group BAM (sorted, indexed)" }
    merged_bai:         { name: "BAI index" }
    group_name:         { name: "Group label (e.g. control, drug_a)" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    File   merged_bam
    File   merged_bai
    String group_name
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 8
  Int mem_gb    = 16
  Int disk_size = ceil(size(merged_bam, "GB") * 2 + 10)

  command <<<
    set -euo pipefail

    # -d MM : MM 태그를 가진 read만 유지 (태그 없는 read가 modkit worker를 deadlock시킴)
    samtools view -b -d MM -@ ~{threads} "~{merged_bam}" -o "~{group_name}.mm_tagged.bam"
    samtools index -@ ~{threads} "~{group_name}.mm_tagged.bam"

    kept=$(samtools view -c "~{group_name}.mm_tagged.bam")
    echo "MM-tagged reads kept: ${kept}"
  >>>

  output {
    File   filtered_bam = "~{group_name}.mm_tagged.bam"
    File   filtered_bai = "~{group_name}.mm_tagged.bam.bai"
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

# ── Task 3a: ipdSummary — unknown 포함 전체 변형 two-sample 탐지 ──────────────
task ipd_summary {
  meta { description: "treated vs control per-base IPD ratio → unknown modification 탐지 (ipdSummary)" }

  parameter_meta {
    treated_bam:        { name: "Pooled treated group BAM (e.g. drug_a.merged.bam)" }
    treated_bai:        { name: "BAI for treated_bam" }
    control_bam:        { name: "Pooled control group BAM" }
    control_bai:        { name: "BAI for control_bam" }
    ref_fasta:          { name: "GRCm39 reference FASTA" }
    comparison_name:    { name: "Comparison label (e.g. drug_a_vs_ctrl)" }
    ipd_pvalue:         { name: "P-value cutoff for modification calling (default 0.001)" }
    min_coverage:       { name: "Minimum read coverage per site (default 25)" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    File   treated_bam
    File   treated_bai
    File   control_bam
    File   control_bai
    File   ref_fasta
    String comparison_name
    Float  ipd_pvalue   = 0.001
    Int    min_coverage = 25
    # Revio HiFi에는 전용 IPD 모델이 없어 Sequel 세대 모델을 강제 사용.
    # two-sample(--control) 차등 분석에서는 모델 편향이 case/control 간 상쇄됨.
    String ipd_model    = "/opt/kt37/lib/python3.7/site-packages/kineticsTools/resources/SP3-C3.npz.gz"
    # 선택적 reference window 제한 (예: "chr19" 또는 "chr1:1-1000000"). 빈 문자열이면 전체 게놈.
    String reference_window = ""
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 16
  Int mem_gb    = 48
  Int disk_size = ceil((size(treated_bam, "GB") + size(control_bam, "GB")) * 2 + size(ref_fasta, "GB") + 20)

  command <<<
    set -euo pipefail
    ipdSummary --version 2>/dev/null || python3 -c "import kineticsTools; print(kineticsTools.__version__)" || true

    # reference_window 지정 시 -w 옵션 추가 (전체 게놈 스캔 방지)
    WIN_ARG=""
    if [ -n "~{reference_window}" ]; then
      WIN_ARG="-w ~{reference_window}"
    fi

    ipdSummary \
      "~{treated_bam}" \
      --reference "~{ref_fasta}" \
      --control "~{control_bam}" \
      --gff "~{comparison_name}.modifications.gff" \
      --csv "~{comparison_name}.modifications.csv" \
      --pvalue ~{ipd_pvalue} \
      --minCoverage ~{min_coverage} \
      --numWorkers ~{threads} \
      --ipdModel "~{ipd_model}" \
      ${WIN_ARG} \
      --log-level INFO

    wc -l "~{comparison_name}.modifications.csv" | awk '{print "modification_sites="$1-1}'
    echo "GFF records: $(grep -vc '^#' ~{comparison_name}.modifications.gff || true)"
  >>>

  output {
    File   modifications_gff = "~{comparison_name}.modifications.gff"
    File   modifications_csv = "~{comparison_name}.modifications.csv"
  }

  runtime {
    docker:                "local/kineticstools:0.8.0"
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

# ── Task 3b: modkit pileup — 5mC / 5hmC / 6mA per-site 집계 ─────────────────
task jasmine_pileup {
  meta { description: "jasmine MM/ML 태그 BAM → modkit pileup으로 5mC / 5hmC / 6mA per-site 확률 BED 생성" }

  parameter_meta {
    merged_bam:         { name: "Pooled jasmine-modified BAM" }
    merged_bai:         { name: "BAI index" }
    ref_fasta:          { name: "GRCm39 reference FASTA" }
    ref_fai:            { name: "FAI index" }
    group_name:         { name: "Group label (e.g. control, drug_a)" }
    min_coverage:       { name: "Minimum coverage for pileup (default 5)" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    File   merged_bam
    File   merged_bai
    File   ref_fasta
    File   ref_fai
    String group_name
    Int    min_coverage = 5
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 8
  Int mem_gb    = 24
  Int disk_size = ceil(size(merged_bam, "GB") * 4 + size(ref_fasta, "GB") + 10)

  command <<<
    set -euo pipefail
    modkit --version

    # 포유류 전장유전체 per-site pileup은 비압축 시 그룹당 수백 GB(특히 6mA).
    # modkit stdout("-") → coverage 필터 → mod code별 분리 → bgzip 스트리밍 압축.
    # 거대한 비압축 중간파일(raw/all_mods)을 만들지 않아 디스크를 절약한다.
    # bedMethyl 컬럼 4 = modification_code (m=5mC, h=5hmC, a=6mA), 컬럼 10 = N_valid_cov
    modkit pileup \
      "~{merged_bam}" \
      - \
      --reference "~{ref_fasta}" \
      --threads ~{threads} \
      --log-filepath "~{group_name}.modkit.log" \
    | awk -v mc=~{min_coverage} \
          -v m5="~{group_name}.5mC.bed.gz" \
          -v h5="~{group_name}.5hmC.bed.gz" \
          -v a6="~{group_name}.6mA.bed.gz" '
        $10 >= mc {
          if      ($4=="m") print | ("bgzip -c > " m5)
          else if ($4=="h") print | ("bgzip -c > " h5)
          else if ($4=="a") print | ("bgzip -c > " a6)
        }'

    # 빈 파일 방지 (해당 변형이 없을 때 빈 gz 생성)
    for f in "~{group_name}.5mC.bed.gz" "~{group_name}.5hmC.bed.gz" "~{group_name}.6mA.bed.gz"; do
      [ -f "${f}" ] || : | bgzip -c > "${f}"
    done

    echo "5mC sites:  $(bgzip -dc ~{group_name}.5mC.bed.gz | wc -l)"
    echo "5hmC sites: $(bgzip -dc ~{group_name}.5hmC.bed.gz | wc -l)"
    echo "6mA sites:  $(bgzip -dc ~{group_name}.6mA.bed.gz | wc -l)"
  >>>

  output {
    File   fivemC_bed    = "~{group_name}.5mC.bed.gz"
    File   fivehmC_bed   = "~{group_name}.5hmC.bed.gz"
    File   sixmA_bed     = "~{group_name}.6mA.bed.gz"
  }

  runtime {
    docker:                "ontresearch/modkit:latest"
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

# ── Task 4a: unknown modification 후보 필터링 및 분류 ─────────────────────────
task filter_and_classify {
  meta { description: "ipdSummary CSV에서 known 변형 제외 후 unknown 후보 추출 및 8-oxoG context 표시" }

  parameter_meta {
    modifications_csv:  { name: "ipdSummary modifications.csv" }
    ref_fasta:          { name: "GRCm39 reference FASTA (context 추출용)" }
    group_label:        { name: "Comparison label (e.g. drug_a_vs_ctrl)" }
    min_ipd_ratio:      { name: "Minimum IPD ratio cutoff (default 2.0)" }
    min_mod_score:      { name: "Minimum modificationScore (default 30)" }
    min_coverage:       { name: "Minimum coverage (default 25)" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    File   modifications_csv
    File   ref_fasta
    String group_label
    Float  min_ipd_ratio  = 2.0
    Int    min_mod_score  = 30
    Int    min_coverage   = 25
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 8
  Int disk_size = ceil(size(modifications_csv, "GB") * 3 + size(ref_fasta, "GB") + 5)

  command <<<
    set -euo pipefail

    python3 << 'PYEOF'
import csv
import sys
import re

csv_path       = "~{modifications_csv}"
ref_fasta_path = "~{ref_fasta}"
group_label    = "~{group_label}"
min_ipd_ratio  = float("~{min_ipd_ratio}")
min_mod_score  = int("~{min_mod_score}")
min_coverage   = int("~{min_coverage}")

# 레퍼런스 FASTA 로드 (context 추출용)
ref = {}
current_chrom = None
buf = []
with open(ref_fasta_path) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if current_chrom:
                ref[current_chrom] = "".join(buf).upper()
            current_chrom = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
if current_chrom:
    ref[current_chrom] = "".join(buf).upper()

def get_context(chrom, pos, width=3):
    seq = ref.get(chrom, "")
    s = max(0, pos - width)
    e = min(len(seq), pos + width + 1)
    return seq[s:e]

unknown_out   = open(f"{group_label}.unknown_candidates.tsv", "w")
known_out     = open(f"{group_label}.known_mod_summary.tsv", "w")
ggg_out       = open(f"{group_label}.GGG_flagged.tsv", "w")

header = ["chrom","pos","strand","coverage","ipdRatio","modificationScore","context","is_CpG","context_type","GGG_flag"]
for fh in [unknown_out, known_out, ggg_out]:
    fh.write("\t".join(header) + "\n")

with open(csv_path) as f:
    reader = csv.DictReader(f)
    n_total = n_unknown = n_known = n_ggg = 0
    for row in reader:
        try:
            chrom    = row.get("refName","") or row.get("seqname","") or row.get("chrom","")
            pos      = int(row.get("tpl","0") or row.get("position","0") or "0")
            strand   = row.get("strand",".")
            cov      = int(float(row.get("coverage","0") or "0"))
            ipd_r    = float(row.get("ipdRatio","0") or "0")
            mod_sc   = float(row.get("modificationScore","0") or "0")
        except (ValueError, KeyError):
            continue

        n_total += 1
        if cov < min_coverage:
            continue

        ctx = get_context(chrom, pos)
        is_cpg = "CG" in ctx[len(ctx)//2:len(ctx)//2+2]
        ggg_flag = "GGG" in ctx.upper()

        # known 변형 분류
        # CpG context + high ipdRatio → likely 5mC/5hmC (handled by jasmine_pileup)
        # canonical 6mA context: GAATTC (EcoRI) → in mouse genome, m6A not canonical
        is_cpg_context = is_cpg
        context_type = "CpG" if is_cpg_context else ("G-rich" if ggg_flag else "other")

        line_vals = [chrom, str(pos), strand, str(cov), f"{ipd_r:.4f}", f"{mod_sc:.1f}",
                     ctx, str(is_cpg_context), context_type, str(ggg_flag)]

        if ipd_r >= min_ipd_ratio and mod_sc >= min_mod_score:
            if is_cpg_context:
                # CpG → 5mC/5hmC 기여 가능성, known_mod_summary로
                known_out.write("\t".join(line_vals) + "\n")
                n_known += 1
            else:
                unknown_out.write("\t".join(line_vals) + "\n")
                n_unknown += 1
                if ggg_flag:
                    ggg_out.write("\t".join(line_vals) + "\n")
                    n_ggg += 1

for fh in [unknown_out, known_out, ggg_out]:
    fh.close()

print(f"total_sites_read={n_total}")
print(f"unknown_candidates={n_unknown}")
print(f"known_cpg_sites={n_known}")
print(f"GGG_flagged_8oxoG_candidates={n_ggg}")
PYEOF
  >>>

  output {
    File   unknown_candidates_tsv = "~{group_label}.unknown_candidates.tsv"
    File   known_mod_summary_tsv  = "~{group_label}.known_mod_summary.tsv"
    File   ggg_flagged_tsv        = "~{group_label}.GGG_flagged.tsv"
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

# ── Task 4b: 4-group known modification 비교 ─────────────────────────────────
task compare_known_mods {
  meta { description: "4개 그룹 5mC / 5hmC / 6mA per-site 수준 비교 → differential BED 및 5hmC/(5mC+5hmC) 비율" }

  parameter_meta {
    group_names:        { name: "Array of 4 group labels [control, drug_a, drug_b, drug_ab]" }
    fivemC_beds:        { name: "Array of 5mC BED files (same order as group_names)" }
    fivehmC_beds:       { name: "Array of 5hmC BED files" }
    sixmA_beds:         { name: "Array of 6mA BED files" }
    out_prefix:         { name: "Output prefix" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    Array[String] group_names
    Array[File]   fivemC_beds
    Array[File]   fivehmC_beds
    Array[File]   sixmA_beds
    String        out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int threads      = 4
  Int mem_gb       = 12
  Int sort_mem_gb  = 4
  # 6mA beds cover every valid adenine genome-wide (not just CpGs); measured ~10.4x gzip
  # ratio and ~93% of rows at 0% modification, so the *2 (vs *5 for mc/hmc) plus filtering
  # below keeps sorted temp space well under the *5 naive estimate.
  Int disk_size    = ceil(size(fivemC_beds, "GB") * 5 + size(fivehmC_beds, "GB") * 5 + size(sixmA_beds, "GB") * 2 + 30)

  # NOTE: pileup bed.gz shards are not globally sorted (region-parallel output simply
  # concatenated), so each is sorted on disk first, then a streaming k-way merge builds
  # the diff table without ever holding all sites in memory. 6mA beds are genome-wide
  # (every A, not just CpGs) and blew past both the old per-task memory budget (loaded as
  # Python dicts) and available disk (~694GB decompressed across 4 groups vs ~400GB free)
  # -- ~93% of 6mA rows are 0% modified background, so those are dropped before sorting;
  # 5mC/5hmC are small enough to keep in full.
  command <<<
    set -euo pipefail
    export LC_ALL=C

    group_names=(~{sep=' ' group_names})
    mc_beds=(~{sep=' ' fivemC_beds})
    hmc_beds=(~{sep=' ' fivehmC_beds})
    ma_beds=(~{sep=' ' sixmA_beds})

    mkdir -p sorted
    sort_one() {
      zcat "$1" | sort -s -t $'\t' -T sorted -S ~{sort_mem_gb}G --parallel=~{threads} -k1,1 -k2,2n -k6,6 | gzip -1 > "$2"
    }
    sort_one_nonzero() {
      zcat "$1" | awk -F'\t' '$11+0 > 0' | sort -s -t $'\t' -T sorted -S ~{sort_mem_gb}G --parallel=~{threads} -k1,1 -k2,2n -k6,6 | gzip -1 > "$2"
    }

    for i in "${!group_names[@]}"; do
      sort_one "${mc_beds[$i]}"  "sorted/g${i}_mc.bed.gz"
      sort_one "${hmc_beds[$i]}" "sorted/g${i}_hmc.bed.gz"
      sort_one_nonzero "${ma_beds[$i]}"  "sorted/g${i}_ma.bed.gz"
    done

    python3 << 'PYEOF'
import heapq, gzip, os

group_names = "~{sep=',' group_names}".split(",")
out_prefix  = "~{out_prefix}"
n = len(group_names)

def line_iter(path):
    if not os.path.exists(path):
        return
    with gzip.open(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 11:
                continue
            try:
                chrom  = cols[0]
                pos    = int(cols[1])
                strand = cols[5]
                frac   = float(cols[10])   # percent_modified
                cov    = int(cols[9])
            except (ValueError, IndexError):
                continue
            yield (chrom, pos, strand), (frac, cov)

# k-way merge over the 3n sorted bed streams: only a handful of "current front" rows
# per file live in memory at once, never the whole genome.
iters = {}
heap = []
for i in range(n):
    for kind in ("mc", "hmc", "ma"):
        name = (kind, i)
        it = line_iter(f"sorted/g{i}_{kind}.bed.gz")
        iters[name] = it
        try:
            k, v = next(it)
            heapq.heappush(heap, (k, name, v))
        except StopIteration:
            pass

ctrl_idx = 0  # control is first group
total_sites = 0

with open(f"{out_prefix}.differential_known_mods.tsv", "w") as out:
    header = ["chrom", "pos", "strand"]
    for g in group_names:
        header += [f"{g}_5mC", f"{g}_5hmC", f"{g}_6mA",
                   f"{g}_5hmC_ratio", f"{g}_5mC_cov", f"{g}_5hmC_cov"]
    header += ["5mC_delta_drugA_vs_ctrl", "5hmC_delta_drugA_vs_ctrl",
               "5hmC_ratio_delta_drugA_vs_ctrl",
               "5mC_delta_rescue", "5hmC_delta_rescue"]
    out.write("\t".join(header) + "\n")

    def write_row(site, buf):
        chrom, pos, strand = site
        row = [chrom, str(pos), strand]
        frac_mc, frac_hmc = [], []
        for i in range(n):
            mc_v,  mc_c  = buf.get(("mc",  i), (0.0, 0))
            hmc_v, hmc_c = buf.get(("hmc", i), (0.0, 0))
            ma_v,  _     = buf.get(("ma",  i), (0.0, 0))
            ratio = hmc_v / (mc_v + hmc_v) if (mc_v + hmc_v) > 0 else 0.0
            row += [f"{mc_v:.4f}", f"{hmc_v:.4f}", f"{ma_v:.4f}",
                    f"{ratio:.4f}", str(mc_c), str(hmc_c)]
            frac_mc.append(mc_v)
            frac_hmc.append(hmc_v)

        # drug_a is index 1, drug_ab is index 3
        da_mc,  da_hmc  = frac_mc[1],  frac_hmc[1]
        ctrl_mc_v       = frac_mc[ctrl_idx]
        ctrl_hmc_v      = frac_hmc[ctrl_idx]
        ab_mc,  ab_hmc  = frac_mc[3],  frac_hmc[3]

        ctrl_mc_hmc  = ctrl_mc_v + ctrl_hmc_v
        da_mc_hmc    = da_mc + da_hmc
        ctrl_ratio   = ctrl_hmc_v / ctrl_mc_hmc if ctrl_mc_hmc > 0 else 0.0
        da_ratio     = da_hmc    / da_mc_hmc    if da_mc_hmc > 0   else 0.0

        row += [
            f"{da_mc - ctrl_mc_v:.4f}",    # 5mC delta A-ctrl
            f"{da_hmc - ctrl_hmc_v:.4f}",  # 5hmC delta A-ctrl
            f"{da_ratio - ctrl_ratio:.4f}",# ratio delta
            f"{ab_mc - da_mc:.4f}",        # rescue 5mC (AB-A)
            f"{ab_hmc - da_hmc:.4f}",      # rescue 5hmC (AB-A)
        ]
        out.write("\t".join(row) + "\n")

    current_key = None
    buffer = {}
    while heap:
        k, name, v = heapq.heappop(heap)
        if current_key is not None and k != current_key:
            write_row(current_key, buffer)
            total_sites += 1
            buffer = {}
        current_key = k
        buffer[name] = v
        it = iters[name]
        try:
            nk, nv = next(it)
            heapq.heappush(heap, (nk, name, nv))
        except StopIteration:
            pass
    if current_key is not None:
        write_row(current_key, buffer)
        total_sites += 1

print(f"total_sites_compared={total_sites}")
print("Output: {}.differential_known_mods.tsv".format(out_prefix))
PYEOF
  >>>

  output {
    File   differential_mods_tsv = "~{out_prefix}.differential_known_mods.tsv"
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

# ── Task 4b-2: differential_known_mods TSV → Parquet 변환 ────────────────────
task known_mods_to_parquet {
  meta { description: "compare_known_mods TSV(대부분 배경 0값, 파일럿 기준 수억 행)를 DuckDB로 컬럼 압축 Parquet 변환 -- 이후 필터링/시각화를 위한 다운스트림 분석 입력" }

  parameter_meta {
    differential_mods_tsv: { name: "compare_known_mods 결과 TSV" }
    out_prefix:             { name: "Output prefix" }
    runtime_attributes:     { name: "Runtime attribute structure" }
  }

  input {
    File   differential_mods_tsv
    String out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 8
  Int mem_gb    = 32
  # Parquet(zstd)은 배경 0값 압축이 잘 되어 보통 TSV의 1/10~1/20 크기가 되지만,
  # DuckDB 변환 중 임시 작업공간을 넉넉히 확보
  Int disk_size = ceil(size(differential_mods_tsv, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    pip install --quiet duckdb 2>/dev/null || true

    python3 << 'PYEOF'
import duckdb

tsv        = "~{differential_mods_tsv}"
out_prefix = "~{out_prefix}"

duckdb.sql("SET threads TO ~{threads};")
duckdb.sql(f"""
    COPY (
        SELECT * FROM read_csv('{tsv}', delim='\t', header=true)
    ) TO '{out_prefix}.differential_known_mods.parquet' (FORMAT PARQUET, COMPRESSION zstd)
""")

n = duckdb.sql(f"SELECT count(*) FROM '{out_prefix}.differential_known_mods.parquet'").fetchone()[0]
print(f"parquet_rows={n}")
print("Output: {}.differential_known_mods.parquet".format(out_prefix))
PYEOF
  >>>

  output {
    File differential_mods_parquet = "~{out_prefix}.differential_known_mods.parquet"
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

# ── Task 4b-3: known mods 필터링 + 시각화 ─────────────────────────────────────
task filter_and_visualize_known_mods {
  meta { description: "known_mods_to_parquet 결과에서 커버리지/델타 임계값으로 유의 사이트를 걸러내고, genome-wide 개요 + rescue 산점도 + 델타 분포를 렌더링" }

  parameter_meta {
    differential_mods_parquet: { name: "known_mods_to_parquet 결과 Parquet" }
    ref_fai:                    { name: "참조 게놈 .fai (염색체 순서/길이, x축 정렬용)" }
    out_prefix:                  { name: "Output prefix" }
    min_total_coverage:          { name: "genome-wide 개요에 포함할 최소 총 커버리지 (4그룹 mc_cov 합)" }
    min_5mC_delta:                { name: "유의 사이트 판정: |5mC_delta_drugA_vs_ctrl| 임계값 (percent-point)" }
    min_5hmC_ratio_delta:        { name: "유의 사이트 판정: |5hmC_ratio_delta_drugA_vs_ctrl| 임계값 (0-1 비율)" }
    bin_size:                    { name: "genome-wide 개요 플롯의 bin 크기 (bp)" }
    scatter_point_cap:           { name: "산점도 렌더링 시 다운샘플링 상한 (전체 개수는 TSV에 그대로 유지)" }
    runtime_attributes:          { name: "Runtime attribute structure" }
  }

  input {
    File   differential_mods_parquet
    File   ref_fai
    String out_prefix
    Int    min_total_coverage   = 10
    Float  min_5mC_delta        = 10.0
    Float  min_5hmC_ratio_delta = 0.1
    Int    bin_size             = 1000000
    Int    scatter_point_cap    = 200000
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 4
  Int mem_gb    = 8
  Int disk_size = ceil(size(differential_mods_parquet, "GB") * 3 + 10)

  command <<<
    set -euo pipefail

    pip install --quiet duckdb matplotlib 2>/dev/null || true

    python3 << 'PYEOF'
import duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import random

parquet   = "~{differential_mods_parquet}"
fai_path  = "~{ref_fai}"
out_prefix = "~{out_prefix}"
min_cov   = ~{min_total_coverage}
min_5mc   = ~{min_5mC_delta}
min_ratio = ~{min_5hmC_ratio_delta}
bin_size  = ~{bin_size}
point_cap = ~{scatter_point_cap}

con = duckdb.connect()
con.sql("SET threads TO ~{threads};")

# ── 염색체 순서/길이 (x축 정렬 -- Parquet 정렬 순서와 동일한 lexicographic 순서를 가정하지 않고,
#    fai에 등장하는 순서를 그대로 사용) ──────────────────────────────────────
chrom_order, chrom_len = [], {}
with open(fai_path) as f:
    for line in f:
        cols = line.rstrip("\n").split("\t")
        chrom_order.append(cols[0])
        chrom_len[cols[0]] = int(cols[1])

chrom_offset = {}
cum = 0
for c in chrom_order:
    chrom_offset[c] = cum
    cum += chrom_len[c]
genome_size = cum

# ── Tier 1: 커버리지만 걸러 genome-wide 개요용 bin 집계 ──────────────────────
bins = con.sql(f"""
    SELECT chrom,
           CAST(pos / {bin_size} AS INTEGER) AS bin,
           avg(abs(\"5mC_delta_drugA_vs_ctrl\"))        AS mean_5mC_delta,
           avg(abs(\"5hmC_ratio_delta_drugA_vs_ctrl\")) AS mean_5hmC_ratio_delta,
           count(*) AS n_sites
    FROM read_parquet('{parquet}')
    WHERE (control_5mC_cov + drug_a_5mC_cov + drug_b_5mC_cov + drug_ab_5mC_cov) >= {min_cov}
    GROUP BY chrom, bin
""").df()

def plot_manhattan(df, value_col, title, out_png):
    fig, ax = plt.subplots(figsize=(16, 4))
    colors = ["#1f77b4", "#ff7f0e"]
    for i, c in enumerate(chrom_order):
        sub = df[df["chrom"] == c]
        if sub.empty:
            continue
        x = chrom_offset[c] + sub["bin"].astype("int64") * bin_size
        ax.scatter(x, sub[value_col], s=3, color=colors[i % 2], linewidths=0)
    ax.set_xlabel("Genomic position (bp, concatenated chroms)")
    ax.set_ylabel(value_col)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

plot_manhattan(bins, "mean_5mC_delta",
                f"Genome-wide mean |5mC delta| (drug_a vs ctrl), {bin_size}bp bins, cov>={min_cov}",
                f"{out_prefix}.manhattan_5mC_delta.png")
plot_manhattan(bins, "mean_5hmC_ratio_delta",
                f"Genome-wide mean |5hmC ratio delta| (drug_a vs ctrl), {bin_size}bp bins, cov>={min_cov}",
                f"{out_prefix}.manhattan_5hmC_ratio_delta.png")

# ── Tier 2: 커버리지 + 델타 임계값을 모두 만족하는 유의 사이트 ──────────────
sig_query = f"""
    SELECT * FROM read_parquet('{parquet}')
    WHERE (control_5mC_cov + drug_a_5mC_cov + drug_b_5mC_cov + drug_ab_5mC_cov) >= {min_cov}
      AND (abs(\"5mC_delta_drugA_vs_ctrl\") >= {min_5mc}
           OR abs(\"5hmC_ratio_delta_drugA_vs_ctrl\") >= {min_ratio})
"""
con.sql(f"COPY ({sig_query}) TO '{out_prefix}.significant_known_mods.tsv' (FORMAT CSV, DELIMITER '\t', HEADER)")
sig = con.sql(sig_query).df()
n_sig = len(sig)

plot_df = sig if n_sig <= point_cap else sig.sample(n=point_cap, random_state=0)

# rescue 산점도: drug_a 변화가 drug_ab에서 얼마나 되돌아오는지
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(plot_df["5mC_delta_drugA_vs_ctrl"], plot_df["5mC_delta_rescue"],
           s=4, alpha=0.3, linewidths=0)
ax.axhline(0, color="grey", linewidth=0.5)
ax.axvline(0, color="grey", linewidth=0.5)
ax.set_xlabel("5mC delta (drug_a vs control)")
ax.set_ylabel("5mC delta (rescue: drug_ab vs drug_a)")
ax.set_title(f"5mC rescue effect (n={n_sig}, plotted={len(plot_df)})")
fig.tight_layout()
fig.savefig(f"{out_prefix}.rescue_scatter.png", dpi=150)
plt.close(fig)

# 델타 분포 히스토그램
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].hist(sig["5mC_delta_drugA_vs_ctrl"].dropna(), bins=100)
axes[0].set_title("5mC delta (drug_a vs ctrl) distribution")
axes[1].hist(sig["5hmC_ratio_delta_drugA_vs_ctrl"].dropna(), bins=100)
axes[1].set_title("5hmC ratio delta (drug_a vs ctrl) distribution")
fig.tight_layout()
fig.savefig(f"{out_prefix}.delta_histograms.png", dpi=150)
plt.close(fig)

with open(f"{out_prefix}.analysis_summary.txt", "w") as f:
    f.write(f"input_parquet={parquet}\n")
    f.write(f"genome_size_bp={genome_size}\n")
    f.write(f"min_total_coverage={min_cov}\n")
    f.write(f"min_5mC_delta={min_5mc}\n")
    f.write(f"min_5hmC_ratio_delta={min_ratio}\n")
    f.write(f"n_significant_sites={n_sig}\n")
    f.write(f"n_bins_genome_overview={len(bins)}\n")

print(f"n_significant_sites={n_sig}")
print(f"Output: {out_prefix}.significant_known_mods.tsv")
PYEOF
  >>>

  output {
    File significant_known_mods_tsv = "~{out_prefix}.significant_known_mods.tsv"
    File manhattan_5mC_delta_png    = "~{out_prefix}.manhattan_5mC_delta.png"
    File manhattan_5hmC_ratio_png   = "~{out_prefix}.manhattan_5hmC_ratio_delta.png"
    File rescue_scatter_png         = "~{out_prefix}.rescue_scatter.png"
    File delta_histograms_png       = "~{out_prefix}.delta_histograms.png"
    File analysis_summary_txt       = "~{out_prefix}.analysis_summary.txt"
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

# ── Task 4c: UMAP+HDBSCAN unknown modification 클러스터링 ─────────────────────
task cluster_unknown_modifications {
  meta { description: "4개 그룹 unknown 후보 통합 후 IPD 프로파일 + 서열 context로 UMAP+HDBSCAN 클러스터링" }

  parameter_meta {
    unknown_tsvs:       { name: "Array of unknown_candidates.tsv from filter_and_classify (4 comparisons)" }
    group_labels:       { name: "Array of comparison labels (same order)" }
    out_prefix:         { name: "Output prefix" }
    min_cluster_size:   { name: "HDBSCAN min_cluster_size (default 50)" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    Array[File]   unknown_tsvs
    Array[String] group_labels
    String        out_prefix
    Int           min_cluster_size = 50
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 4
  Int mem_gb    = 16
  Int disk_size = ceil(size(unknown_tsvs, "GB") * 5 + 10)

  command <<<
    set -euo pipefail

    # umap-learn / hdbscan がない場合はインストール
    pip install --quiet umap-learn hdbscan seaborn scikit-learn 2>/dev/null || true

    python3 << 'PYEOF'
import sys, os, csv
import numpy as np

# 의존 패키지 동적 임포트
try:
    import umap
    import hdbscan as hdbscan_mod
    import seaborn as sns
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.preprocessing import StandardScaler
    HAS_UMAP = True
except ImportError as e:
    print(f"WARNING: {e}. Clustering will be skipped; output TSV will still be generated.")
    HAS_UMAP = False

tsv_paths   = "~{sep=',' unknown_tsvs}".split(",")
group_labels= "~{sep=',' group_labels}".split(",")
out_prefix  = "~{out_prefix}"
min_cs      = int("~{min_cluster_size}")

# 통합 로드
rows = []
for path, label in zip(tsv_paths, group_labels):
    if not os.path.exists(path):
        continue
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            r["_group"] = label
            rows.append(r)

print(f"total_candidates_loaded={len(rows)}")

if len(rows) == 0:
    print("No candidates. Creating empty output files.")
    open(f"{out_prefix}.clustered_modifications.tsv","w").write("chrom\tpos\tstrand\tcluster_id\tgroup\n")
    open(f"{out_prefix}.cluster_summary.tsv","w").write("cluster_id\tsize\tmean_ipdRatio\ttop_context\tgroup_dist\n")
    sys.exit(0)

# 특징 벡터 구성
def context_encode(ctx):
    """3-mer at center → 4-hot (A,C,G,T) × 3 = 12 dims"""
    mid = max(0, len(ctx)//2 - 1)
    kmer = ctx[mid:mid+3].upper().ljust(3, "N")
    bases = "ACGT"
    vec = []
    for b in kmer:
        vec += [1 if b == x else 0 for x in bases]
    return vec

features = []
valid_rows = []
for r in rows:
    try:
        ipd  = float(r.get("ipdRatio","0") or "0")
        sc   = float(r.get("modificationScore","0") or "0")
        ctx  = r.get("context","NNN")
        ggg  = 1 if r.get("GGG_flag","False") == "True" else 0
        feat = [ipd, sc/100.0, ggg] + context_encode(ctx)
        features.append(feat)
        valid_rows.append(r)
    except (ValueError, KeyError):
        continue

X = np.array(features, dtype=float)
X = StandardScaler().fit_transform(X)

cluster_ids = np.full(len(X), -1, dtype=int)

if HAS_UMAP and len(X) >= min_cs * 2:
    reducer   = umap.UMAP(n_components=2, random_state=42, n_jobs=~{threads})
    embedding = reducer.fit_transform(X)

    clusterer  = hdbscan_mod.HDBSCAN(min_cluster_size=min_cs, prediction_data=True)
    cluster_ids = clusterer.fit_predict(embedding)

    # UMAP scatter plot
    fig, ax = plt.subplots(figsize=(10, 8))
    unique_ids = sorted(set(cluster_ids))
    palette    = sns.color_palette("husl", len(unique_ids))
    for cid, color in zip(unique_ids, palette):
        mask = cluster_ids == cid
        lbl  = f"Cluster {cid}" if cid >= 0 else "Noise"
        ax.scatter(embedding[mask, 0], embedding[mask, 1],
                   c=[color], label=lbl, s=5, alpha=0.6)
    ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
    ax.set_title(f"Unknown Modification Clusters — {out_prefix}")
    ax.legend(markerscale=3, bbox_to_anchor=(1, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}.umap_clusters.pdf", dpi=150)
    plt.close()
    print(f"n_clusters={len(set(cluster_ids)) - (1 if -1 in cluster_ids else 0)}")
    print(f"noise_points={int((cluster_ids == -1).sum())}")
else:
    print("UMAP skipped (too few points or package missing) — all assigned cluster 0")
    cluster_ids = np.zeros(len(valid_rows), dtype=int)

# 출력 TSV
with open(f"{out_prefix}.clustered_modifications.tsv", "w") as out:
    out.write("chrom\tpos\tstrand\tcoverage\tipdRatio\tmodificationScore\tcontext\tGGG_flag\tgroup\tcluster_id\n")
    for r, cid in zip(valid_rows, cluster_ids):
        out.write(f"{r.get('chrom','')}\t{r.get('pos','')}\t{r.get('strand','.')}\t"
                  f"{r.get('coverage','')}\t{r.get('ipdRatio','')}\t"
                  f"{r.get('modificationScore','')}\t{r.get('context','')}\t"
                  f"{r.get('GGG_flag','')}\t{r['_group']}\t{cid}\n")

# 클러스터 요약
from collections import Counter, defaultdict
cluster_rows = defaultdict(list)
for r, cid in zip(valid_rows, cluster_ids):
    cluster_rows[int(cid)].append(r)

with open(f"{out_prefix}.cluster_summary.tsv", "w") as out:
    out.write("cluster_id\tsize\tmean_ipdRatio\tmean_modScore\ttop_context\tGGG_pct\tgroup_distribution\tputative_modification\n")
    for cid in sorted(cluster_rows.keys()):
        rlist = cluster_rows[cid]
        ipds  = [float(r.get("ipdRatio","0") or "0") for r in rlist]
        scs   = [float(r.get("modificationScore","0") or "0") for r in rlist]
        ctxs  = Counter(r.get("context","") for r in rlist)
        ggg_n = sum(1 for r in rlist if r.get("GGG_flag","False")=="True")
        grps  = Counter(r["_group"] for r in rlist)
        top_ctx = ctxs.most_common(1)[0][0] if ctxs else ""
        ggg_pct = ggg_n / len(rlist) * 100 if rlist else 0
        mean_ipd = np.mean(ipds) if ipds else 0
        # 추정 변형 유형 휴리스틱
        if ggg_pct > 50 and mean_ipd > 3:
            putative = "8-oxoG_candidate"
        elif mean_ipd > 5:
            putative = "FapyG_FapyA_candidate"
        elif cid == -1:
            putative = "noise"
        else:
            putative = "unknown"
        out.write(f"{cid}\t{len(rlist)}\t{mean_ipd:.3f}\t{np.mean(scs):.1f}\t{top_ctx}\t"
                  f"{ggg_pct:.1f}\t{dict(grps)}\t{putative}\n")

print(f"Clustering done: {out_prefix}.clustered_modifications.tsv")
PYEOF
  >>>

  output {
    File   clustered_tsv    = "~{out_prefix}.clustered_modifications.tsv"
    File   cluster_summary  = "~{out_prefix}.cluster_summary.tsv"
    File?  umap_plot        = "~{out_prefix}.umap_clusters.pdf"
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

# ── Task 5: 4-group 통합 비교 + rescue analysis ───────────────────────────────
task multigroup_comparison {
  meta { description: "unknown + known 변형 4개 그룹 비교 → rescue 후보 및 5hmC_rescue 추출" }

  parameter_meta {
    unknown_tsvs:          { name: "Array[4] unknown_candidates.tsv (A_vs_ctrl, B_vs_ctrl, AB_vs_ctrl, AB_vs_A)" }
    comparison_labels:     { name: "Array[4] comparison labels matching unknown_tsvs" }
    differential_mods_tsv: { name: "compare_known_mods 결과 TSV" }
    clustered_tsv:         { name: "cluster_unknown_modifications 결과 TSV" }
    out_prefix:            { name: "Output prefix" }
    min_ipd_ratio:         { name: "IPD ratio cutoff for Drug A effect (default 2.0)" }
    runtime_attributes:    { name: "Runtime attribute structure" }
  }

  input {
    Array[File]   unknown_tsvs
    Array[String] comparison_labels
    File          differential_mods_tsv
    File          clustered_tsv
    String        out_prefix
    Float         min_ipd_ratio = 2.0
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 16
  Int disk_size = ceil(size(unknown_tsvs, "GB") * 4 + 10)

  command <<<
    set -euo pipefail

    python3 << 'PYEOF'
import csv, os
from collections import defaultdict

tsvs        = "~{sep=',' unknown_tsvs}".split(",")
labels      = "~{sep=',' comparison_labels}".split(",")
diff_tsv    = "~{differential_mods_tsv}"
clust_tsv   = "~{clustered_tsv}"
out_prefix  = "~{out_prefix}"
min_ipd     = float("~{min_ipd_ratio}")

def load_unknown(path, label):
    sites = {}
    if not os.path.exists(path):
        return sites
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            key = (r["chrom"], r["pos"], r["strand"])
            sites[key] = {"label": label, **r}
    return sites

# 기대 레이블: a_vs_ctrl, b_vs_ctrl, ab_vs_ctrl, ab_vs_a
label_map = {l: d for l, d in zip(labels, [load_unknown(p, l) for p in tsvs])}

# Drug A에서 증가한 위치
a_label  = labels[0]   # a_vs_ctrl
ab_label = labels[2]   # ab_vs_ctrl
aba_label= labels[3]   # ab_vs_a (rescue: Drug A+B vs Drug A)

a_sites   = label_map.get(a_label,  {})
ab_sites  = label_map.get(ab_label, {})
aba_sites = label_map.get(aba_label,{})

# cluster ID 로드
cluster_ids = {}
if os.path.exists(clust_tsv):
    with open(clust_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            key = (r["chrom"], r["pos"], r["strand"])
            cluster_ids[key] = r.get("cluster_id", "-1")

# 5hmC rescue 로드
hmc_rescue_sites = set()
if os.path.exists(diff_tsv):
    with open(diff_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            try:
                hmc_da  = float(r.get("5hmC_delta_drugA_vs_ctrl","0") or "0")
                hmc_res = float(r.get("5hmC_delta_rescue","0") or "0")
                # Drug A에서 5hmC 증가 AND Drug A+B에서 복구(감소)
                if hmc_da > 0.1 and hmc_res < -0.05:
                    hmc_rescue_sites.add((r["chrom"], r["pos"], r["strand"]))
            except (ValueError, KeyError):
                continue

# categorized output
cat_header = ["chrom","pos","strand","ipdRatio_drugA","modScore_drugA","context",
              "GGG_flag","cluster_id","category"]
rescue_header = ["chrom","pos","strand","ipdRatio_drugA","ipdRatio_drugAB","rescue_ratio",
                 "context","GGG_flag","cluster_id"]
hmc_header = ["chrom","pos","strand","5hmC_delta_drugA","5hmC_delta_rescue"]
venn_counts = defaultdict(int)

with (open(f"{out_prefix}.categorized_modifications.tsv","w") as cat_out,
      open(f"{out_prefix}.rescue_candidates.tsv","w") as rescue_out,
      open(f"{out_prefix}.5hmC_rescue.tsv","w") as hmc_out,
      open(f"{out_prefix}.venn_summary.tsv","w") as venn_out):

    cat_out.write("\t".join(cat_header)+"\n")
    rescue_out.write("\t".join(rescue_header)+"\n")
    hmc_out.write("\t".join(hmc_header)+"\n")

    for key, a_row in a_sites.items():
        ipd_a = float(a_row.get("ipdRatio","0") or "0")
        if ipd_a < min_ipd:
            continue
        in_ab   = key in ab_sites
        rescued = key not in aba_sites   # Drug A+B에서 사라진 = rescued
        cid     = cluster_ids.get(key, "-1")
        ctx     = a_row.get("context","")
        ggg     = a_row.get("GGG_flag","")
        sc      = a_row.get("modificationScore","")

        if rescued:
            cat = "rescued_by_B"
            venn_counts["A_only_rescued"] += 1
            ipd_ab = float(aba_sites.get(key, {}).get("ipdRatio","0") or "0")
            rescue_ratio = ipd_ab / ipd_a if ipd_a > 0 else 0
            rescue_out.write("\t".join([key[0], key[1], key[2],
                             f"{ipd_a:.4f}", f"{ipd_ab:.4f}", f"{rescue_ratio:.4f}",
                             ctx, ggg, cid]) + "\n")
        else:
            cat = "A_only_persistent"
            venn_counts["A_persistent"] += 1

        cat_out.write("\t".join([key[0], key[1], key[2],
                     f"{ipd_a:.4f}", sc, ctx, ggg, cid, cat]) + "\n")

    # 5hmC rescue
    if os.path.exists(diff_tsv):
        with open(diff_tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                key = (r["chrom"], r["pos"], r["strand"])
                if key in hmc_rescue_sites:
                    hmc_out.write(f"{r['chrom']}\t{r['pos']}\t{r['strand']}\t"
                                  f"{r.get('5hmC_delta_drugA_vs_ctrl','')}\t"
                                  f"{r.get('5hmC_delta_rescue','')}\n")

    venn_out.write("category\tcount\n")
    for cat, cnt in sorted(venn_counts.items()):
        venn_out.write(f"{cat}\t{cnt}\n")

print(f"Drug A sites (IPD>={min_ipd}): {len([k for k,v in a_sites.items() if float(v.get('ipdRatio','0') or '0') >= min_ipd])}")
print(f"Rescue candidates written to {out_prefix}.rescue_candidates.tsv")
print(f"5hmC rescue sites: {len(hmc_rescue_sites)}")
PYEOF
  >>>

  output {
    File   categorized_tsv     = "~{out_prefix}.categorized_modifications.tsv"
    File   rescue_candidates   = "~{out_prefix}.rescue_candidates.tsv"
    File   fivehmC_rescue_tsv  = "~{out_prefix}.5hmC_rescue.tsv"
    File   venn_summary        = "~{out_prefix}.venn_summary.tsv"
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

# ── Task 6: IPD ratio → BigWig 변환 ──────────────────────────────────────────
task ipd_to_bigwig {
  meta { description: "modifications.csv ipdRatio → sorted BedGraph → BigWig (IGV/UCSC 시각화)" }

  parameter_meta {
    modifications_csv:  { name: "ipdSummary modifications.csv" }
    chrom_sizes:        { name: "Chromosome sizes file (from samtools faidx .fai)" }
    comparison_name:    { name: "Comparison label (출력 파일명에 사용)" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    File   modifications_csv
    File   chrom_sizes
    String comparison_name
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 8
  Int disk_size = ceil(size(modifications_csv, "GB") * 6 + 5)

  command <<<
    set -euo pipefail

    # FAI → chrom.sizes (col1, col2 만 필요)
    awk '{print $1"\t"$2}' "~{chrom_sizes}" > chrom.sizes

    # CSV → BedGraph (chrom, start, end, ipdRatio)
    python3 << 'PYEOF'
import csv, os
csv_path = "~{modifications_csv}"
with open(csv_path) as fin, open("ipdratio.bedgraph","w") as fout:
    reader = csv.DictReader(fin)
    for row in reader:
        try:
            chrom = row.get("refName","") or row.get("seqname","")
            pos   = int(row.get("tpl","0") or row.get("position","0") or "0")
            ipd   = float(row.get("ipdRatio","0") or "0")
            cov   = int(float(row.get("coverage","0") or "0"))
        except (ValueError, KeyError):
            continue
        if cov < 5 or ipd <= 0:
            continue
        fout.write(f"{chrom}\t{pos}\t{pos+1}\t{ipd:.4f}\n")
PYEOF

    # sort → clip → BigWig
    sort -k1,1 -k2,2n ipdratio.bedgraph > ipdratio.sorted.bedgraph
    bedGraphToBigWig ipdratio.sorted.bedgraph chrom.sizes "~{comparison_name}.ipdRatio.bw"
    echo "BigWig created: $(wc -c < ~{comparison_name}.ipdRatio.bw) bytes"
  >>>

  output {
    File   ipdRatio_bw = "~{comparison_name}.ipdRatio.bw"
  }

  runtime {
    docker:                "quay.io/biocontainers/ucsc-bedgraphtobigwig:472--h9b8f530_1"
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

# ── Task 7: MEME motif analysis ───────────────────────────────────────────────
task motif_analysis {
  meta { description: "rescue 후보 주변 서열 추출 후 MEME/STREME motif 탐색" }

  parameter_meta {
    rescue_candidates:  { name: "rescue_candidates.tsv from multigroup_comparison" }
    ref_fasta:          { name: "GRCm39 reference FASTA" }
    drug_a_candidates:  { name: "unknown_candidates.tsv for Drug A (background)" }
    out_prefix:         { name: "Output prefix" }
    window_bp:          { name: "Flanking window around each site (default 250bp)" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    File   rescue_candidates
    File   ref_fasta
    File   drug_a_candidates
    String out_prefix
    Int    window_bp = 250
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 8
  Int mem_gb    = 16
  Int disk_size = ceil(size(ref_fasta, "GB") * 2 + 5)

  command <<<
    set -euo pipefail

    python3 << 'PYEOF'
import csv, os

ref_fasta_path = "~{ref_fasta}"
rescue_path    = "~{rescue_candidates}"
bg_path        = "~{drug_a_candidates}"
out_prefix     = "~{out_prefix}"
win            = int("~{window_bp}")

# 레퍼런스 로드 (전체 로드 — GRCm39 ~2.7GB in memory)
print("Loading reference FASTA...", flush=True)
ref = {}
cur = None
buf = []
with open(ref_fasta_path) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if cur:
                ref[cur] = "".join(buf).upper()
            cur = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
if cur:
    ref[cur] = "".join(buf).upper()
print(f"Reference loaded: {len(ref)} chromosomes", flush=True)

def extract_seq(path, out_fa, max_seqs=5000):
    written = 0
    with open(path) as f, open(out_fa,"w") as fo:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            chrom = r.get("chrom","")
            try:
                pos = int(r.get("pos","0"))
            except ValueError:
                continue
            seq_all = ref.get(chrom,"")
            if not seq_all:
                continue
            s = max(0, pos - win)
            e = min(len(seq_all), pos + win + 1)
            seq = seq_all[s:e]
            if len(seq) < 50:
                continue
            fo.write(f">{chrom}_{pos}\n{seq}\n")
            written += 1
            if written >= max_seqs:
                break
    return written

n_rescue = extract_seq(rescue_path,    f"{out_prefix}.rescue.fa")
n_bg     = extract_seq(bg_path,        f"{out_prefix}.background.fa", max_seqs=10000)
print(f"rescue_seqs={n_rescue}, background_seqs={n_bg}", flush=True)
PYEOF

    # MEME: discovery on rescue sequences
    mkdir -p "~{out_prefix}.meme_out" "~{out_prefix}.streme_out"
    if [ -s "~{out_prefix}.rescue.fa" ]; then
      meme "~{out_prefix}.rescue.fa" \
        -dna -oc "~{out_prefix}.meme_out" \
        -nmotifs 10 -minw 6 -maxw 20 \
        -mod zoops -p ~{threads} \
        -evt 0.05 2>&1 | tail -20 || true

      # STREME: discriminative (rescue vs background)
      if [ -s "~{out_prefix}.background.fa" ]; then
        streme \
          --p "~{out_prefix}.rescue.fa" \
          --n "~{out_prefix}.background.fa" \
          --dna --oc "~{out_prefix}.streme_out" \
          --thresh 0.05 2>&1 | tail -10 || true
      fi
    else
      echo "No rescue sequences — skipping MEME/STREME"
    fi

    tar -czf "~{out_prefix}.meme_results.tar.gz"   "~{out_prefix}.meme_out"
    tar -czf "~{out_prefix}.streme_results.tar.gz" "~{out_prefix}.streme_out"
    echo "Motif analysis complete"
  >>>

  output {
    File   rescue_fa       = "~{out_prefix}.rescue.fa"
    File   background_fa   = "~{out_prefix}.background.fa"
    File   meme_tar        = "~{out_prefix}.meme_results.tar.gz"
    File   streme_tar      = "~{out_prefix}.streme_results.tar.gz"
  }

  runtime {
    docker:                "quay.io/biocontainers/meme:5.5.7--pl5321hca22f42_0"
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

# ── Task 8: 변형 위치 유전체 어노테이션 ──────────────────────────────────────
task annotate_modifications {
  meta { description: "rescue/categorized 후보에 Gencode GTF 및 CpG island 어노테이션 부착" }

  parameter_meta {
    categorized_tsv:    { name: "multigroup_comparison 출력 categorized_modifications.tsv" }
    annotation_gtf_bgz: { name: "Sorted bgzipped Gencode GTF" }
    annotation_gtf_tbi: { name: "Tabix index" }
    cpg_island_tsv:     { name: "CpG island coordinates TSV (chrom,start,end,name)" }
    out_prefix:         { name: "Output prefix" }
    window_bp:          { name: "TSS upstream window for promoter classification (default 2000)" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    File   categorized_tsv
    File   annotation_gtf_bgz
    File   annotation_gtf_tbi
    File   cpg_island_tsv
    String out_prefix
    Int    window_bp = 2000
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 8
  Int disk_size = ceil(size(categorized_tsv, "GB") * 5 + size(annotation_gtf_bgz, "GB") * 3 + 5)

  command <<<
    set -euo pipefail

    # categorized TSV → BED (chrom, start, end, name)
    awk 'NR>1{print $1"\t"$2"\t"$2+1"\t"NR-1"\t"$4"\t"$3}' \
      "~{categorized_tsv}" > candidates.bed

    # CpG island TSV → BED
    awk 'NR>1{print $1"\t"$2"\t"$3"\t"$4}' \
      "~{cpg_island_tsv}" > cpg_islands.bed

    # GTF → gene BED (gene level only)
    zcat "~{annotation_gtf_bgz}" 2>/dev/null || gzip -dc "~{annotation_gtf_bgz}" | \
      awk '$3=="gene"{
        match($0, /gene_name "([^"]+)"/, a)
        match($0, /gene_type "([^"]+)"/, b)
        print $1"\t"$4-1"\t"$5"\t"a[1]"\t"b[1]"\t"$7
      }' > genes.bed 2>/dev/null || \
    tabix "~{annotation_gtf_bgz}" -h | \
      awk '$3=="gene"{
        match($0, /gene_name "([^"]+)"/, a)
        match($0, /gene_type "([^"]+)"/, b)
        print $1"\t"$4-1"\t"$5"\t"a[1]"\t"b[1]"\t"$7
      }' > genes.bed

    # GTF → exon BED
    (zcat "~{annotation_gtf_bgz}" 2>/dev/null || gzip -dc "~{annotation_gtf_bgz}") | \
      awk '$3=="exon"{
        match($0, /gene_name "([^"]+)"/, a)
        print $1"\t"$4-1"\t"$5"\t"a[1]"\t.\t"$7
      }' > exons.bed 2>/dev/null || true

    # bedtools intersect — gene
    bedtools intersect -a candidates.bed -b genes.bed -wao \
      | awk 'BEGIN{OFS="\t"} {
          if($NF==0) {gene_name="intergenic"; gene_type="."; gene_strand="."}
          else       {gene_name=$10; gene_type=$11; gene_strand=$12}
          print $1,$2,$3,$4,$5,$6,gene_name,gene_type,gene_strand
        }' > with_gene.bed

    # bedtools intersect — exon
    bedtools intersect -a with_gene.bed -b exons.bed -wao \
      | awk 'BEGIN{OFS="\t"} {
          if($NF==0) feature="intronic/intergenic"
          else       feature="exonic"
          print $0, feature
        }' > with_exon.bed

    # bedtools intersect — CpG island
    bedtools intersect -a with_exon.bed -b cpg_islands.bed -wao \
      | awk 'BEGIN{OFS="\t"} {
          if($NF==0) cpg_isl="no"
          else       cpg_isl="yes"
          print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,cpg_isl
        }' > annotated.bed

    # promoter flag (upstream ~{window_bp}bp of gene TSS)
    python3 << 'PYEOF'
import csv

win = int("~{window_bp}")
out_prefix = "~{out_prefix}"

# gene TSS map: gene_name -> (chrom, tss, strand)
tss_map = {}
try:
    with open("genes.bed") as f:
        for line in f:
            cols = line.rstrip().split("\t")
            if len(cols) < 6:
                continue
            chrom, start, end, name, gtype, strand = cols[:6]
            start, end = int(start), int(end)
            tss = start if strand == "+" else end
            tss_map[name] = (chrom, tss, strand)
except FileNotFoundError:
    pass

header = ["chrom","pos","end","site_id","ipdRatio","strand",
          "gene_name","gene_type","gene_strand","feature","cpg_island",
          "is_promoter","dist_to_TSS"]

with open("annotated.bed") as fin, \
     open(f"{out_prefix}.annotated_modifications.tsv","w") as fout:
    fout.write("\t".join(header)+"\n")
    for line in fin:
        cols = line.rstrip().split("\t")
        if len(cols) < 11:
            continue
        chrom, pos, end = cols[0], cols[1], cols[2]
        site_id, ipd, strand = cols[3], cols[4], cols[5]
        gene_name, gene_type, gene_strand = cols[6], cols[7], cols[8]
        feature, cpg_isl = cols[9], cols[10]
        # dist to TSS
        is_promoter = "no"
        dist_tss = "."
        if gene_name in tss_map:
            gc, tss, gs = tss_map[gene_name]
            if gc == chrom:
                d = int(pos) - tss
                dist_tss = str(abs(d))
                if abs(d) <= win:
                    is_promoter = "yes"
                    feature = "promoter"
        fout.write("\t".join([chrom,pos,end,site_id,ipd,strand,
                              gene_name,gene_type,gene_strand,
                              feature,cpg_isl,is_promoter,dist_tss])+"\n")

print(f"Annotated: {out_prefix}.annotated_modifications.tsv")
PYEOF
  >>>

  output {
    File   annotated_tsv = "~{out_prefix}.annotated_modifications.tsv"
  }

  runtime {
    docker:                "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2"
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

# ── Task 9: GO/KEGG pathway enrichment ───────────────────────────────────────
task pathway_enrichment {
  meta { description: "rescue/annotated 변형 유전자 GO/KEGG enrichment (gseapy)" }

  parameter_meta {
    annotated_tsv:      { name: "annotate_modifications 출력 TSV" }
    out_prefix:         { name: "Output prefix" }
    species:            { name: "Species for gseapy (default: Mouse)" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    File   annotated_tsv
    String out_prefix
    String species = "Mouse"
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 8
  Int disk_size = ceil(size(annotated_tsv, "GB") * 3 + 5)

  command <<<
    set -euo pipefail

    python3 << 'PYEOF'
import csv, os, sys

try:
    import gseapy as gp
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_GSEAPY = True
except ImportError:
    print("gseapy not found — installing")
    os.system("pip install --quiet gseapy")
    try:
        import gseapy as gp
        HAS_GSEAPY = True
    except ImportError:
        HAS_GSEAPY = False

out_prefix = "~{out_prefix}"
species    = "~{species}"

# 유전자 목록 추출 (rescue 후보 위치의 유전자)
genes = set()
with open("~{annotated_tsv}") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        g = r.get("gene_name","")
        cat = r.get("feature","")
        if g and g != "intergenic" and g != ".":
            genes.add(g)

gene_list = sorted(genes)
print(f"Input genes: {len(gene_list)}")

with open(f"{out_prefix}.gene_list.txt","w") as f:
    f.write("\n".join(gene_list)+"\n")

if not HAS_GSEAPY or len(gene_list) < 5:
    print("Skipping enrichment (gseapy unavailable or too few genes)")
    open(f"{out_prefix}.go_enrichment.tsv","w").write("Term\tGenes\tP-value\n")
    open(f"{out_prefix}.kegg_enrichment.tsv","w").write("Term\tGenes\tP-value\n")
    sys.exit(0)

# GO enrichment (Biological Process)
gene_sets_go   = "GO_Biological_Process_2021"
gene_sets_kegg = "KEGG_2019_Mouse" if species == "Mouse" else "KEGG_2021_Human"

for gs_name, out_name in [(gene_sets_go, "go"), (gene_sets_kegg, "kegg")]:
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gs_name,
            organism=species,
            outdir=None,
            cutoff=0.05
        )
        if enr is not None and hasattr(enr, "results") and enr.results is not None:
            df = enr.results
            df = df[df["Adjusted P-value"] < 0.05].head(30)
            df.to_csv(f"{out_prefix}.{out_name}_enrichment.tsv", sep="\t", index=False)
            print(f"{out_name} significant terms: {len(df)}")

            # dot plot
            if len(df) > 0:
                try:
                    fig, ax = plt.subplots(figsize=(8, max(4, len(df)*0.3 + 2)))
                    df_plot = df.head(20).sort_values("Adjusted P-value")
                    import numpy as np
                    ax.scatter(
                        -np.log10(df_plot["Adjusted P-value"].astype(float)),
                        range(len(df_plot)),
                        s=50, c="steelblue", alpha=0.7
                    )
                    ax.set_yticks(range(len(df_plot)))
                    ax.set_yticklabels(df_plot["Term"].tolist(), fontsize=8)
                    ax.set_xlabel("-log10(adj. P-value)")
                    ax.set_title(f"{gs_name} Enrichment — {out_prefix}")
                    plt.tight_layout()
                    plt.savefig(f"{out_prefix}.{out_name}_enrichment.pdf", dpi=150)
                    plt.close()
                except Exception as e:
                    print(f"Plot failed: {e}")
        else:
            open(f"{out_prefix}.{out_name}_enrichment.tsv","w").write("Term\tGenes\tP-value\n")
    except Exception as e:
        print(f"{gs_name} enrichment failed: {e}")
        open(f"{out_prefix}.{out_name}_enrichment.tsv","w").write("Term\tGenes\tP-value\n")
PYEOF
  >>>

  output {
    File   gene_list_txt    = "~{out_prefix}.gene_list.txt"
    File   go_enrichment    = "~{out_prefix}.go_enrichment.tsv"
    File   kegg_enrichment  = "~{out_prefix}.kegg_enrichment.tsv"
    File?  go_plot          = "~{out_prefix}.go_enrichment.pdf"
    File?  kegg_plot        = "~{out_prefix}.kegg_enrichment.pdf"
  }

  runtime {
    docker:                "quay.io/biocontainers/gseapy:1.3.0--py311heb3b1e3_0"
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

# ── Task 10: 통합 HTML 보고서 생성 ────────────────────────────────────────────
task modification_report {
  meta { description: "모든 분석 결과 통합 HTML 보고서 + 최종 rescue 우선순위 목록" }

  parameter_meta {
    experiment_id:          { name: "Experiment ID" }
    rescue_candidates:      { name: "rescue_candidates.tsv" }
    fivehmC_rescue_tsv:     { name: "5hmC_rescue.tsv" }
    venn_summary:           { name: "venn_summary.tsv" }
    cluster_summary:        { name: "cluster_summary.tsv" }
    annotated_tsv:          { name: "annotated_modifications.tsv" }
    go_enrichment:          { name: "go_enrichment.tsv" }
    kegg_enrichment:        { name: "kegg_enrichment.tsv" }
    differential_mods_tsv:  { name: "differential_known_mods.tsv" }
    group_names:            { name: "Array[4] group names [control, drug_a, drug_b, drug_ab]" }
    run_hifimeth_qc:        { name: "Whether hifimeth QC was run (default false)" }
    out_prefix:             { name: "Output prefix" }
    runtime_attributes:     { name: "Runtime attribute structure" }
  }

  input {
    String        experiment_id
    File          rescue_candidates
    File          fivehmC_rescue_tsv
    File          venn_summary
    File          cluster_summary
    File          annotated_tsv
    File          go_enrichment
    File          kegg_enrichment
    File          differential_mods_tsv
    Array[String] group_names
    Boolean       run_hifimeth_qc = false
    String        out_prefix
    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 8
  Int disk_size = 10

  command <<<
    set -euo pipefail

    python3 << 'PYEOF'
import csv, os, sys, json
from datetime import date

exp_id        = "~{experiment_id}"
out_prefix    = "~{out_prefix}"
group_names   = "~{sep=',' group_names}".split(",")
hifimeth_qc   = "~{run_hifimeth_qc}".lower() == "true"
today         = date.today().isoformat()

def load_tsv_rows(path, max_rows=None):
    rows = []
    if not os.path.exists(path):
        return [], []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = reader.fieldnames or []
        for i, r in enumerate(reader):
            if max_rows and i >= max_rows:
                break
            rows.append(r)
    return headers, rows

def rows_to_html_table(headers, rows, max_rows=20, cls=""):
    if not rows:
        return "<p><em>(no data)</em></p>"
    shown = rows[:max_rows]
    html = f'<table class="table {cls}"><thead><tr>'
    for h in headers:
        html += f"<th>{h}</th>"
    html += "</tr></thead><tbody>"
    for r in shown:
        html += "<tr>" + "".join(f"<td>{r.get(h,'')}</td>" for h in headers) + "</tr>"
    if len(rows) > max_rows:
        html += f'<tr><td colspan="{len(headers)}"><em>... {len(rows)-max_rows} more rows</em></td></tr>'
    html += "</tbody></table>"
    return html

# Load key stats
hdr_diff, diff_rows = load_tsv_rows("~{differential_mods_tsv}", max_rows=200)
_, rescue_rows    = load_tsv_rows("~{rescue_candidates}")
_, hmc_rows       = load_tsv_rows("~{fivehmC_rescue_tsv}")
_, venn_rows      = load_tsv_rows("~{venn_summary}")
_, cluster_rows   = load_tsv_rows("~{cluster_summary}")
hdr_ann, ann_rows = load_tsv_rows("~{annotated_tsv}", max_rows=200)
hdr_go,  go_rows  = load_tsv_rows("~{go_enrichment}")
hdr_kegg,kegg_rows= load_tsv_rows("~{kegg_enrichment}")

n_rescue   = len(rescue_rows)
n_hmc      = len(hmc_rows)
n_clusters = sum(1 for r in cluster_rows if r.get("cluster_id","") != "-1")
gene_set   = set(r.get("gene_name","") for r in ann_rows
                 if r.get("gene_name","") not in ("","intergenic","."))
n_genes    = len(gene_set)

# rescue priority (top 20 by ipdRatio)
try:
    rescue_sorted = sorted(rescue_rows,
                           key=lambda r: -float(r.get("ipdRatio_drugA","0") or "0"))
except Exception:
    rescue_sorted = rescue_rows
priority_list = rescue_sorted[:50]

# venn text
venn_text = " | ".join(f"{r.get('category','')}: {r.get('count','')}" for r in venn_rows)

html = f"""<!DOCTYPE html>
<html lang="ko"><head>
<meta charset="UTF-8">
<title>Base Modification Discovery Report — {exp_id}</title>
<style>
body{{font-family:Arial,sans-serif;margin:2em;background:#f8f9fa;}}
h1{{color:#1a237e;}} h2{{color:#283593;border-bottom:2px solid #3f51b5;padding-bottom:4px;}}
h3{{color:#3f51b5;}}
.card{{background:#fff;border-radius:8px;box-shadow:0 2px 6px rgba(0,0,0,0.1);
      padding:1em 1.5em;margin:1em 0;}}
.stat-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:1em;margin:1em 0;}}
.stat-box{{background:#e8eaf6;border-radius:8px;padding:1em;text-align:center;}}
.stat-num{{font-size:2em;font-weight:bold;color:#1a237e;}}
.stat-label{{font-size:0.85em;color:#555;}}
table{{border-collapse:collapse;width:100%;font-size:0.85em;}}
th{{background:#3f51b5;color:#fff;padding:6px 8px;text-align:left;}}
td{{padding:5px 8px;border-bottom:1px solid #ddd;}}
tr:nth-child(even){{background:#f5f5f5;}}
.tag{{display:inline-block;padding:2px 8px;border-radius:4px;font-size:0.8em;font-weight:bold;}}
.tag-rescue{{background:#e8f5e9;color:#2e7d32;}}
.tag-hmC{{background:#e3f2fd;color:#1565c0;}}
.note{{color:#888;font-size:0.85em;font-style:italic;}}
</style>
</head><body>

<h1>Base Modification Discovery Report</h1>
<div class="card">
<p><strong>Experiment:</strong> {exp_id} &nbsp;|&nbsp;
   <strong>Generated:</strong> {today} &nbsp;|&nbsp;
   <strong>Groups:</strong> {", ".join(group_names)}</p>
</div>

<h2>Summary Statistics</h2>
<div class="stat-grid">
  <div class="stat-box"><div class="stat-num">{n_rescue}</div>
    <div class="stat-label">Unknown Rescue Candidates<br>(Drug A → Drug A+B reversal)</div></div>
  <div class="stat-box"><div class="stat-num">{n_hmc}</div>
    <div class="stat-label">5hmC Rescue Sites<br>(TET/ROS-linked)</div></div>
  <div class="stat-box"><div class="stat-num">{n_clusters}</div>
    <div class="stat-label">Unknown Mod Clusters<br>(UMAP+HDBSCAN)</div></div>
  <div class="stat-box"><div class="stat-num">{n_genes}</div>
    <div class="stat-label">Genes with Modifications</div></div>
  <div class="stat-box"><div class="stat-num">{len(go_rows)}</div>
    <div class="stat-label">GO Terms (adj.P&lt;0.05)</div></div>
</div>
<p class="note">Venn: {venn_text or "(no data)"}</p>

<h2>Top Rescue Candidates (Unknown Modifications)</h2>
<div class="card">
<p class="tag tag-rescue">Drug A ↑ → Drug A+B reversal</p>
{rows_to_html_table(["chrom","pos","strand","ipdRatio_drugA","ipdRatio_drugAB","rescue_ratio","context","GGG_flag","cluster_id"],
                    priority_list[:20])}
</div>

<h2>5hmC Rescue (TET Activity Reversal)</h2>
<div class="card">
<p class="tag tag-hmC">5hmC increased in Drug A, rescued by Drug B</p>
{rows_to_html_table(["chrom","pos","strand","5hmC_delta_drugA","5hmC_delta_rescue"],
                    hmc_rows, max_rows=20)}
</div>

<h2>Known Modification Differentials (5mC / 5hmC / 6mA, top 50)</h2>
<div class="card">
{rows_to_html_table(hdr_diff, diff_rows, max_rows=50)}
</div>

<h2>Unknown Modification Clusters</h2>
<div class="card">
{rows_to_html_table(list(cluster_rows[0].keys()) if cluster_rows else [],
                    cluster_rows)}
<p class="note">Putative types: 8-oxoG (GGG context, high IPD), FapyG/FapyA (purine, asymmetric, very high IPD), unknown</p>
</div>

<h2>Annotated Modification Sites (top 50)</h2>
<div class="card">
{rows_to_html_table(hdr_ann, ann_rows, max_rows=50)}
</div>

<h2>GO Biological Process Enrichment</h2>
<div class="card">
{rows_to_html_table(hdr_go, go_rows, max_rows=20)}
<p class="note">Prioritized pathways: BER (base excision repair), Nrf2/KEAP1, NF-κB, TET/DNA demethylation</p>
</div>

<h2>KEGG Pathway Enrichment</h2>
<div class="card">
{rows_to_html_table(hdr_kegg, kegg_rows, max_rows=20)}
</div>

<h2>Experimental Validation Recommendations</h2>
<div class="card">
<ul>
<li><strong>8-oxoG candidates (GGG-context clusters):</strong> Dot blot with anti-8-oxoG antibody; OGG1-based lesion enrichment + sequencing</li>
<li><strong>5hmC rescue sites:</strong> hMeDIP-seq targeted validation; TAB-seq for single-base 5hmC mapping</li>
<li><strong>FapyG/FapyA candidates (high IPD, purine):</strong> LC-MS/MS quantification; Formamidopyrimidine DNA glycosylase (Fpg) digestion assay</li>
<li><strong>All rescue candidates (top priority):</strong> CUT&amp;RUN with mod-specific antibodies at priority loci</li>
{"<li><strong>5mC QC (hifimeth):</strong> hifimeth concordance table available</li>" if hifimeth_qc else ""}
</ul>
</div>

</body></html>"""

with open(f"{out_prefix}.base_mod_report.html","w") as f:
    f.write(html)
print(f"HTML report: {out_prefix}.base_mod_report.html")

# rescue priority TSV (final deliverable)
priority_header = ["chrom","pos","strand","ipdRatio_drugA","ipdRatio_drugAB","rescue_ratio",
                   "context","GGG_flag","cluster_id"]
with open(f"{out_prefix}.rescue_priority_list.tsv","w") as f:
    f.write("\t".join(priority_header)+"\n")
    for r in priority_list:
        f.write("\t".join(r.get(h,"") for h in priority_header)+"\n")
print(f"Priority list: {out_prefix}.rescue_priority_list.tsv ({len(priority_list)} sites)")
PYEOF
  >>>

  output {
    File   report_html          = "~{out_prefix}.base_mod_report.html"
    File   rescue_priority_list = "~{out_prefix}.rescue_priority_list.tsv"
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
