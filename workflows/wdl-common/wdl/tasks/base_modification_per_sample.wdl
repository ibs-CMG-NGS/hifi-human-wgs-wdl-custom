version 1.0

import "../structs.wdl"

# per-sample 5hmC/5mC 경로용 task 모음 (base_modification_per_sample.wdl).
# 기존 base_modification.wdl task는 call cache 해시 보존을 위해 수정하지 않는다.

# ── Task 1: 미정렬 raw HiFi BAM에 jasmine 직접 적용 ────────────────────────────
task jasmine_unaligned {
  meta { description: "미정렬 raw HiFi BAM(fi/fp/ri/rp) → jasmine 5mC/5hmC/6mA 호출. kinetics는 출력에서 제거되어 원본의 ~18% 크기" }

  parameter_meta {
    raw_hifi_bam:       { name: "Raw HiFi BAM (kinetics 포함, 미정렬)" }
    sample_id:          { name: "Sample ID" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    File   raw_hifi_bam
    String sample_id
    Int    threads = 32
    RuntimeAttributes runtime_attributes
  }

  Int mem_gb    = 48
  Int disk_size = ceil(size(raw_hifi_bam, "GB") * 0.4 + 20)

  command <<<
    set -euo pipefail
    jasmine --version

    # 기본값: --5hmC, --6mA 활성, kinetics 트랙 제거(--keep-kinetics 미지정)
    jasmine \
      --log-level INFO \
      --num-threads ~{threads} \
      "~{raw_hifi_bam}" \
      "~{sample_id}.jasmine.bam"
  >>>

  output {
    File jasmine_bam = "~{sample_id}.jasmine.bam"
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

# ── Task 2: MM 태그 없는 리드 제거 + pbmm2 정렬 ────────────────────────────────
task pbmm2_align_modbam {
  meta { description: "jasmine BAM에서 MM 태그 없는 리드(~0.4%, modkit 0.6.3 데드락 유발) 제거 후 GRCm39 정렬. 필터 사본은 task 안에서 삭제" }

  parameter_meta {
    jasmine_bam:        { name: "jasmine 출력 BAM (미정렬, MM/ML 포함)" }
    ref_fasta:          { name: "Reference FASTA (GRCm39)" }
    ref_fai:            { name: "FAI index" }
    sample_id:          { name: "Sample ID" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    File   jasmine_bam
    File   ref_fasta
    File   ref_fai
    String sample_id
    Int    threads = 32
    RuntimeAttributes runtime_attributes
  }

  Int mem_gb    = 80
  Int disk_size = ceil(size(jasmine_bam, "GB") * 3 + size(ref_fasta, "GB") * 2 + 20)

  command <<<
    set -euo pipefail
    pbmm2 --version
    samtools --version | head -1

    # pbmm2 sort spill이 루트 /tmp(여유 ~100G)로 가지 않도록 work dir 안에 둔다
    mkdir -p tmp
    export TMPDIR="$PWD/tmp"

    total=$(samtools view -c -@ 8 "~{jasmine_bam}")
    samtools view -@ 8 -b -d MM -o mm_only.bam "~{jasmine_bam}"
    kept=$(samtools view -c -@ 8 mm_only.bam)
    echo "reads_total=${total} reads_with_MM=${kept}"

    pbmm2 align \
      --preset HIFI \
      --sort \
      --sort-memory 4G \
      --sort-threads 8 \
      --sample "~{sample_id}" \
      --num-threads ~{threads} \
      --log-level INFO \
      "~{ref_fasta}" \
      mm_only.bam \
      "~{sample_id}.GRCm39.modbam.bam"

    rm -f mm_only.bam
    rm -rf tmp
  >>>

  output {
    File modbam = "~{sample_id}.GRCm39.modbam.bam"
    File modbai = "~{sample_id}.GRCm39.modbam.bam.bai"
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

# ── Task 3: CpG 5mC/5hmC per-sample pileup (DSS 입력 형식) ─────────────────────
task modbam_cpg_pileup {
  meta { description: "modkit pileup --cpg --combine-strands, 전 샘플 공통 고정 threshold → 5mC/5hmC 4컬럼(chr, pos, N_valid, N_mod) bgzip" }

  parameter_meta {
    modbam:             { name: "정렬된 MM/ML BAM (MM 없는 리드 제거됨)" }
    modbai:             { name: "BAI index" }
    ref_fasta:          { name: "Reference FASTA" }
    ref_fai:            { name: "FAI index" }
    sample_id:          { name: "Sample ID" }
    min_coverage:       { name: "최소 N_valid_cov" }
    filter_threshold:   { name: "modkit 고정 pass threshold (샘플별 자동 추정 대신 공통값)" }
    runtime_attributes: { name: "Runtime attribute structure" }
  }

  input {
    File   modbam
    File   modbai
    File   ref_fasta
    File   ref_fai
    String sample_id
    Int    min_coverage     = 5
    Float  filter_threshold = 0.7
    Int    threads          = 8
    RuntimeAttributes runtime_attributes
  }

  Int mem_gb    = 24
  Int disk_size = ceil(size(modbam, "GB") * 0.2 + size(ref_fasta, "GB") + 10)

  command <<<
    set -euo pipefail
    modkit --version

    # modkit 0.6.3은 jasmine 26 BAM(MM: C+m?;A+a.;T-a.;C+h?;G-h?)에서 --cpg / --motif /
    # --modified-bases 중 하나라도 주면 모든 read가 "failed processing"으로 빠진다.
    # 그래서 필터 없이 pileup한 뒤 5mC/5hmC 행만 골라 참조 CpG에서 두 가닥을 합친다
    # (--cpg --combine-strands와 같은 결과: + 가닥 C(p)와 - 가닥 G(p+1)를 p로 합산).
    cat > cpg_combine.py <<'PY'
    import subprocess, sys
    ref, mc, m5, h5 = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
    fai = {l.split("\t")[0]: [int(x) for x in l.split("\t")[1:5]] for l in open(ref + ".fai")}
    def load(c):
        n, off, lb, lw = fai[c]
        with open(ref, "rb") as f:
            f.seek(off)
            raw = f.read(n + (n + lb - 1) // lb * (lw - lb))
        return raw.replace(b"\n", b"").replace(b"\r", b"").upper()[:n]
    outs = {m: subprocess.Popen(f"bgzip -c > {p}", shell=True, stdin=subprocess.PIPE, text=True)
            for m, p in (("m", m5), ("h", h5))}
    acc, cur, seq = {}, None, b""
    def flush(limit):
        for k in sorted(k for k in acc if k[1] < limit):
            nv, nm = acc.pop(k)
            if nv >= mc:
                outs[k[0]].stdin.write(f"{cur}\t{k[1]}\t{nv}\t{nm}\n")
    for line in sys.stdin:
        f = line.split("\t")
        if f[0] != cur:
            if cur is not None:
                flush(float("inf"))
            cur, seq = f[0], load(f[0])
        p = int(f[1])
        if f[5] == "+" and seq[p:p + 2] == b"CG":
            k = (f[3], p)
        elif f[5] == "-" and p >= 1 and seq[p - 1:p + 1] == b"CG":
            k = (f[3], p - 1)
        else:
            continue
        a = acc.setdefault(k, [0, 0])
        a[0] += int(f[9]); a[1] += int(f[11])
        if len(acc) > 200000:
            flush(p - 10000)
    flush(float("inf"))
    for o in outs.values():
        o.stdin.close(); o.wait()
    PY

    # bedMethyl: 1 chrom, 2 start, 4 mod_code(m/h/a), 6 strand, 10 N_valid_cov, 12 N_mod
    modkit pileup \
      "~{modbam}" \
      - \
      --reference "~{ref_fasta}" \
      --filter-threshold ~{filter_threshold} \
      --threads ~{threads} \
      --log-filepath "~{sample_id}.modkit.log" \
    | awk -F'\t' '$4=="m" || $4=="h"' \
    | python3 cpg_combine.py "~{ref_fasta}" ~{min_coverage} "~{sample_id}.5mC.tsv.gz" "~{sample_id}.5hmC.tsv.gz"

    echo "5mC sites:  $(bgzip -dc ~{sample_id}.5mC.tsv.gz | wc -l)"
    echo "5hmC sites: $(bgzip -dc ~{sample_id}.5hmC.tsv.gz | wc -l)"
  >>>

  output {
    File fivemC_tsv  = "~{sample_id}.5mC.tsv.gz"
    File fivehmC_tsv = "~{sample_id}.5hmC.tsv.gz"
    File modkit_log  = "~{sample_id}.modkit.log"
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
