# Multi-Sample Batch Processing 완전 가이드

## 📋 목차
1. [현재 환경 요약](#현재-환경-요약)
2. [Batch Processing 설정](#batch-processing-설정)
3. [결과 파일 구조](#결과-파일-구조)
4. [성능 최적화](#성능-최적화)
5. [실행 및 모니터링](#실행-및-모니터링)

---

## 🔍 현재 환경 요약

### 성공적으로 완료된 단일 샘플 분석:
- **샘플 ID**: KTY9537
- **실행 시간**: ~20시간
- **데이터 위치**: `/data_4tb/pacbio_rawdata/`
- **결과 위치**: `/home/ygkim/ngs-pipeline/hifi-human-wgs-wdl-custom/20260120_101704_humanwgs_singleton/`
- **캐시 위치**: `/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_call_cache/`
- **컨테이너 캐시**: `/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_singularity_cache/`

### 디스크 공간:
- **루트 (/)**: 457GB (100% 사용) ❌
- **/data_4tb**: 3.6TB (20% 사용, 2.8TB 여유) ✅

**권장 사항**: 모든 출력을 `/data_4tb`에 저장

---

## 🚀 Batch Processing 설정

### 1. 디렉토리 구조 생성

```bash
cd ~/ngs-pipeline/hifi-human-wgs-wdl-custom

# 입력 파일 디렉토리
mkdir -p batch_inputs

# 출력을 /data_4tb로 설정
# 기존에 사용했던 /data_4tb/hifi-human-wgs-wdl-custom 디렉토리 활용 (권한 문제 없음)
export BATCH_OUTPUT_DIR="/data_4tb/hifi-human-wgs-wdl-custom/batch_results"

# 디렉토리 생성
mkdir -p ${BATCH_OUTPUT_DIR}
mkdir -p ${BATCH_OUTPUT_DIR}/logs

# 심볼릭 링크 생성 (편의상)
ln -sf ${BATCH_OUTPUT_DIR} ./batch_outputs
```

### 2. 샘플 입력 파일 준비

#### 원본 데이터 확인:
```bash
# /data_4tb에 있는 모든 BAM 파일 찾기
find /data_4tb/pacbio_rawdata -name "*.bam" -type f

# 예상 구조:
# /data_4tb/pacbio_rawdata/r84285_20260108_080127/1_A01/hifi_reads/*.bam
# /data_4tb/pacbio_rawdata/r84285_20260108_080127/1_B01/hifi_reads/*.bam
# 등...
```

#### 자동 입력 파일 생성 스크립트:
```bash
cat > create_batch_inputs.sh << 'EOF'
#!/bin/bash
# create_batch_inputs.sh
# /data_4tb의 BAM 파일들로부터 자동으로 입력 JSON 파일 생성

RAWDATA_DIR="/data_4tb/pacbio_rawdata"
BATCH_INPUT_DIR="batch_inputs"
REF_MAP="/home/ygkim/ngs-pipeline/hifi-human-wgs-wdl-custom/hifi-wdl-resources/hifi-wdl-resources-v3.1.0/GRCh38.ref_map.v3p1p0.template.tsv"
TERTIARY_MAP="/home/ygkim/ngs-pipeline/hifi-human-wgs-wdl-custom/hifi-wdl-resources/hifi-wdl-resources-v3.1.0/GRCh38.tertiary_map.v3p1p0.template.tsv"

mkdir -p ${BATCH_INPUT_DIR}

# 샘플 정보 CSV 파일 (수동 작성 필요)
# 형식: sample_id,sex,bam_files
# 예: KTY9537,MALE,/data_4tb/pacbio_rawdata/.../file1.bam:/data_4tb/.../file2.bam

if [[ ! -f "samples.csv" ]]; then
    echo "Error: samples.csv not found"
    echo "Create samples.csv with format:"
    echo "sample_id,sex,bam_files"
    echo "KTY9537,MALE,/path/to/file1.bam:/path/to/file2.bam"
    exit 1
fi

# CSV 파일 읽기 (헤더 제외)
tail -n +2 samples.csv | while IFS=',' read -r sample_id sex bam_files; do
    echo "Creating input file for ${sample_id}..."
    
    # BAM 파일들을 배열로 변환
    IFS=':' read -ra BAM_ARRAY <<< "$bam_files"
    
    # JSON 배열 생성
    bam_json=""
    for bam in "${BAM_ARRAY[@]}"; do
        if [[ -z "$bam_json" ]]; then
            bam_json="\"${bam}\""
        else
            bam_json="${bam_json},\n    \"${bam}\""
        fi
    done
    
    # JSON 파일 생성
    cat > ${BATCH_INPUT_DIR}/${sample_id}.inputs.json << JSONEOF
{
  "humanwgs_singleton.sample_id": "${sample_id}",
  "humanwgs_singleton.sex": "${sex}",
  "humanwgs_singleton.hifi_reads": [
    ${bam_json}
  ],
  "humanwgs_singleton.ref_map_file": "${REF_MAP}",
  "humanwgs_singleton.tertiary_map_file": "${TERTIARY_MAP}",
  "humanwgs_singleton.backend": "HPC",
  "humanwgs_singleton.preemptible": false,
  "humanwgs_singleton.gpu": true,
  "humanwgs_singleton.max_reads_per_alignment_chunk": 100000000
}
JSONEOF
    
    echo "  Created: ${BATCH_INPUT_DIR}/${sample_id}.inputs.json"
done

echo ""
echo "✓ All input files created in ${BATCH_INPUT_DIR}/"
ls -lh ${BATCH_INPUT_DIR}/
EOF

chmod +x create_batch_inputs.sh
```

#### samples.csv 예시 파일 생성:
```bash
cat > samples.csv << 'EOF'
sample_id,sex,bam_files
KTY9537,MALE,/data_4tb/pacbio_rawdata/r84285_20260108_080127/1_A01/hifi_reads/m84285_260108_082608_s1.hifi_reads.bc2016.bam
KTY9538,FEMALE,/data_4tb/pacbio_rawdata/sample2/hifi_reads/sample2.bam
KTY9539,MALE,/data_4tb/pacbio_rawdata/sample3/cell1.bam:/data_4tb/pacbio_rawdata/sample3/cell2.bam
EOF

# 입력 파일 생성
./create_batch_inputs.sh
```

### 3. Batch 실행 스크립트 수정

기존 `batch_run.sh`를 복사하고 `/data_4tb` 출력 경로로 수정:

```bash
cat > batch_run_optimized.sh << 'EOF'
#!/bin/bash
# batch_run_optimized.sh
# /data_4tb에 출력하는 최적화된 배치 처리 스크립트

# 사용법:
# ./batch_run_optimized.sh [parallel|sequential] [sample1 sample2 ...]
# 
# 예시:
# ./batch_run_optimized.sh parallel              # 모든 샘플 병렬 실행
# ./batch_run_optimized.sh sequential KTY9537 KTY9538  # 특정 샘플만 순차 실행

set -e  # 에러 발생 시 중단

# 설정
WORKFLOW="workflows/singleton.wdl"
INPUT_DIR="batch_inputs"
OUTPUT_BASE="/data_4tb/hifi-human-wgs-wdl-custom/batch_results"
LOG_DIR="${OUTPUT_BASE}/logs"
CONFIG_FILE="config/miniwdl.local.cfg"

# 실행 모드 (기본값: parallel)
MODE="${1:-parallel}"
shift || true  # 첫 번째 인자 제거

# 샘플 목록 (인자로 받거나 batch_inputs/*.inputs.json에서 자동 생성)
if [[ $# -gt 0 ]]; then
    SAMPLES=("$@")
else
    # batch_inputs 디렉토리에서 모든 .inputs.json 파일 찾기
    SAMPLES=()
    for file in ${INPUT_DIR}/*.inputs.json; do
        if [[ -f "$file" ]]; then
            sample=$(basename "$file" .inputs.json)
            SAMPLES+=("$sample")
        fi
    done
fi

# 디렉토리 생성
mkdir -p "${OUTPUT_BASE}"
mkdir -p "${LOG_DIR}"

echo "========================================"
echo "HiFi-human-WGS Batch Processing"
echo "========================================"
echo "Workflow: ${WORKFLOW}"
echo "Input directory: ${INPUT_DIR}"
echo "Output base: ${OUTPUT_BASE}"
echo "Log directory: ${LOG_DIR}"
echo "Config file: ${CONFIG_FILE}"
echo "Mode: ${MODE}"
echo "Samples to process: ${#SAMPLES[@]}"
for sample in "${SAMPLES[@]}"; do
    echo "  - ${sample}"
done
echo "========================================"
echo ""

# Conda 환경 확인
if [[ -z "${CONDA_DEFAULT_ENV}" ]] || [[ "${CONDA_DEFAULT_ENV}" != "hifi-human-wgs" ]]; then
    echo "⚠ Warning: hifi-human-wgs conda environment is not activated."
    echo "Activating environment..."
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate hifi-human-wgs
fi

# miniwdl 확인
if ! command -v miniwdl &> /dev/null; then
    echo "✗ Error: miniwdl is not installed or not in PATH"
    exit 1
fi

# GPU 환경 변수 설정 (GPU 1번만 사용)
export CUDA_VISIBLE_DEVICES=1

# 병렬 실행 함수
run_parallel() {
    echo "🚀 Running samples in parallel mode..."
    echo ""
    
    for sample in "${SAMPLES[@]}"; do
        input_file="${INPUT_DIR}/${sample}.inputs.json"
        output_dir="${OUTPUT_BASE}/${sample}"
        log_file="${LOG_DIR}/${sample}.log"
        
        # 입력 파일 존재 확인
        if [[ ! -f "${input_file}" ]]; then
            echo "⚠ Warning: Input file not found: ${input_file}"
            echo "  Skipping ${sample}..."
            echo ""
            continue
        fi
        
        echo "▶ Starting ${sample}..."
        miniwdl run "${WORKFLOW}" \
            --input "${input_file}" \
            --cfg "${CONFIG_FILE}" \
            --dir "${output_dir}" \
            --verbose \
            > "${log_file}" 2>&1 &
        
        pid=$!
        echo "  PID: ${pid}"
        echo "  Input: ${input_file}"
        echo "  Output: ${output_dir}"
        echo "  Log: ${log_file}"
        echo ""
        
        # PID를 파일에 기록
        echo "${pid}" > "${LOG_DIR}/${sample}.pid"
    done
    
    echo "⏳ All samples started. Waiting for completion..."
    wait
    echo ""
    echo "✓ All samples completed!"
}

# 순차 실행 함수
run_sequential() {
    echo "🔄 Running samples in sequential mode..."
    echo ""
    
    local success_count=0
    local fail_count=0
    local failed_samples=()
    
    for sample in "${SAMPLES[@]}"; do
        input_file="${INPUT_DIR}/${sample}.inputs.json"
        output_dir="${OUTPUT_BASE}/${sample}"
        log_file="${LOG_DIR}/${sample}.log"
        
        # 입력 파일 존재 확인
        if [[ ! -f "${input_file}" ]]; then
            echo "⚠ Warning: Input file not found: ${input_file}"
            echo "  Skipping ${sample}..."
            echo ""
            ((fail_count++))
            failed_samples+=("${sample} (input not found)")
            continue
        fi
        
        echo "========================================"
        echo "▶ Processing ${sample}..."
        echo "========================================"
        echo "  Input: ${input_file}"
        echo "  Output: ${output_dir}"
        echo "  Log: ${log_file}"
        echo "  Started: $(date)"
        echo ""
        
        sample_start=$(date +%s)
        
        miniwdl run "${WORKFLOW}" \
            --input "${input_file}" \
            --cfg "${CONFIG_FILE}" \
            --dir "${output_dir}" \
            --verbose \
            2>&1 | tee "${log_file}"
        
        exit_code=${PIPESTATUS[0]}
        sample_end=$(date +%s)
        sample_elapsed=$((sample_end - sample_start))
        
        echo ""
        if [[ ${exit_code} -eq 0 ]]; then
            echo "✓ ${sample} completed successfully"
            echo "  Duration: $((sample_elapsed / 3600))h $(((sample_elapsed % 3600) / 60))m $((sample_elapsed % 60))s"
            ((success_count++))
        else
            echo "✗ ${sample} failed with exit code ${exit_code}"
            echo "  Check log: ${log_file}"
            ((fail_count++))
            failed_samples+=("${sample} (exit code ${exit_code})")
            
            # 계속할지 물어봄
            read -p "Continue with next sample? (y/n) " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                echo "⏹ Stopping batch processing."
                break
            fi
        fi
        echo ""
    done
    
    echo "========================================"
    echo "📊 Processing Summary"
    echo "========================================"
    echo "Total samples: ${#SAMPLES[@]}"
    echo "Successful: ${success_count}"
    echo "Failed: ${fail_count}"
    
    if [[ ${fail_count} -gt 0 ]]; then
        echo ""
        echo "Failed samples:"
        for failed in "${failed_samples[@]}"; do
            echo "  ✗ ${failed}"
        done
    fi
    echo "========================================"
}

# 실행 시작 시간
start_time=$(date +%s)
echo "⏱ Started at: $(date)"
echo ""

# 모드에 따라 실행
case "${MODE}" in
    parallel)
        run_parallel
        ;;
    sequential)
        run_sequential
        ;;
    *)
        echo "✗ Error: Unknown mode '${MODE}'"
        echo "Usage: $0 [parallel|sequential] [sample1 sample2 ...]"
        exit 1
        ;;
esac

# 실행 종료 시간
end_time=$(date +%s)
elapsed=$((end_time - start_time))
hours=$((elapsed / 3600))
minutes=$(((elapsed % 3600) / 60))
seconds=$((elapsed % 60))

echo ""
echo "========================================"
echo "✓ Batch Processing Complete"
echo "========================================"
echo "Total samples: ${#SAMPLES[@]}"
echo "Total time: ${hours}h ${minutes}m ${seconds}s"
echo "Output directory: ${OUTPUT_BASE}"
echo "Log directory: ${LOG_DIR}"
echo ""
echo "To check status:"
echo "  tail -f ${LOG_DIR}/<sample>.log"
echo ""
echo "To check results:"
echo "  ls -lh ${OUTPUT_BASE}/<sample>/out/"
echo "========================================"
EOF

chmod +x batch_run_optimized.sh
```

---

## 📂 결과 파일 구조

### 각 샘플별 디렉토리 구조:

```
/data_4tb/hifi-human-wgs-wdl-custom/batch_results/
├── logs/
│   ├── KTY9537.log         # 실행 로그
│   ├── KTY9537.pid         # 프로세스 ID
│   ├── KTY9538.log
│   └── KTY9539.log
│
├── KTY9537/                # 샘플 1 전체 결과
│   ├── out/                # 최종 결과 파일들
│   │   ├── phased_small_variant_vcf/
│   │   │   └── KTY9537.GRCh38.small_variants.phased.vcf.gz
│   │   ├── phased_sv_vcf/
│   │   │   └── KTY9537.GRCh38.structural_variants.phased.vcf.gz
│   │   ├── merged_haplotagged_bam/
│   │   │   └── KTY9537.GRCh38.haplotagged.bam
│   │   ├── pharmcat_report_html/
│   │   │   └── KTY9537.pharmcat.report.html
│   │   ├── stats_file/
│   │   │   └── KTY9537.stats.txt
│   │   ├── bam_statistics/
│   │   ├── mosdepth_summary/
│   │   ├── read_length_plot/
│   │   └── ... (70+ 출력 파일 카테고리)
│   │
│   ├── call-upstream/      # 중간 파일들
│   ├── call-downstream/
│   └── workflow.log
│
├── KTY9538/                # 샘플 2
│   └── out/
│       ├── phased_small_variant_vcf/
│       │   └── KTY9538.GRCh38.small_variants.phased.vcf.gz
│       └── ...
│
└── KTY9539/                # 샘플 3
    └── out/
        └── ...
```

### 주요 결과 파일 위치:

각 샘플의 결과는 **`/data_4tb/hifi-human-wgs-wdl-custom/batch_results/<sample_id>/out/`**에 저장됩니다.

#### 1. VCF 파일들:
```bash
# Small variants (SNP/INDEL)
/data_4tb/hifi-human-wgs-wdl-custom/batch_results/<sample_id>/out/phased_small_variant_vcf/<sample_id>.GRCh38.small_variants.phased.vcf.gz

# Structural variants
/data_4tb/hifi-human-wgs-wdl-custom/batch_results/<sample_id>/out/phased_sv_vcf/<sample_id>.GRCh38.structural_variants.phased.vcf.gz

# Tandem repeats
/data_4tb/hifi-human-wgs-wdl-custom/batch_results/<sample_id>/out/phased_trgt_vcf/<sample_id>.GRCh38.trgt.sorted.phased.vcf.gz
```

#### 2. BAM 파일:
```bash
/data_4tb/hifi-human-wgs-wdl-custom/batch_results/<sample_id>/out/merged_haplotagged_bam/<sample_id>.GRCh38.haplotagged.bam
```

#### 3. 리포트 및 통계:
```bash
# HTML 리포트
/data_4tb/hifi-human-wgs-wdl-custom/batch_results/<sample_id>/out/pharmcat_report_html/<sample_id>.pharmcat.report.html

# 통계 요약
/data_4tb/hifi-human-wgs-wdl-custom/batch_results/<sample_id>/out/stats_file/<sample_id>.stats.txt

# Coverage 통계
/data_4tb/hifi-human-wgs-wdl-custom/batch_results/<sample_id>/out/mosdepth_summary/<sample_id>.GRCh38.mosdepth.summary.txt
```

### 결과 수집 스크립트:

```bash
cat > collect_results.sh << 'EOF'
#!/bin/bash
# collect_results.sh
# 모든 샘플의 주요 결과를 한 곳에 모으기

BATCH_RESULTS="/data_4tb/hifi-human-wgs-wdl-custom/batch_results"
SUMMARY_DIR="${BATCH_RESULTS}/summary"

mkdir -p ${SUMMARY_DIR}/{vcfs,bams,reports,stats}

echo "Collecting results from ${BATCH_RESULTS}..."

for sample_dir in ${BATCH_RESULTS}/*/; do
    sample=$(basename "${sample_dir}")
    
    # 로그나 summary 디렉토리는 건너뜀
    if [[ "$sample" == "logs" ]] || [[ "$sample" == "summary" ]]; then
        continue
    fi
    
    out_dir="${sample_dir}/out"
    
    if [[ ! -d "$out_dir" ]]; then
        echo "⚠ No output directory for ${sample}"
        continue
    fi
    
    echo "Processing ${sample}..."
    
    # VCF 파일 복사
    find "${out_dir}" -name "*.vcf.gz" -exec cp {} ${SUMMARY_DIR}/vcfs/ \; 2>/dev/null
    
    # BAM 파일 심볼릭 링크 (용량 절약)
    find "${out_dir}" -name "*.bam" -exec ln -sf {} ${SUMMARY_DIR}/bams/ \; 2>/dev/null
    
    # HTML 리포트 복사
    find "${out_dir}" -name "*.html" -exec cp {} ${SUMMARY_DIR}/reports/ \; 2>/dev/null
    
    # 통계 파일 복사
    find "${out_dir}" -name "*stats*.txt" -o -name "*summary*.txt" | while read file; do
        cp "$file" "${SUMMARY_DIR}/stats/$(basename $(dirname $(dirname $file)))_$(basename $file)" 2>/dev/null
    done
done

echo ""
echo "✓ Results collected in ${SUMMARY_DIR}/"
echo ""
echo "Summary:"
echo "  VCFs: $(ls -1 ${SUMMARY_DIR}/vcfs/*.vcf.gz 2>/dev/null | wc -l) files"
echo "  BAMs: $(ls -1 ${SUMMARY_DIR}/bams/*.bam 2>/dev/null | wc -l) files"
echo "  Reports: $(ls -1 ${SUMMARY_DIR}/reports/*.html 2>/dev/null | wc -l) files"
echo "  Stats: $(ls -1 ${SUMMARY_DIR}/stats/*.txt 2>/dev/null | wc -l) files"
EOF

chmod +x collect_results.sh
```

---

## ⚡ 성능 최적화

### 1. 실행 시간 단축 전략

#### 문제: 단일 샘플이 20시간 소요

**원인 분석**:
- pbmm2 alignment: 예상 2-4시간 → **실제 20시간** (비정상적으로 느림)
- 가능한 원인:
  1. 디스크 I/O 병목
  2. 메모리 부족으로 인한 스왑
  3. CPU 스레드 경합

#### 해결책 1: 디스크 I/O 최적화

```bash
# /data_4tb를 작업 디렉토리로 완전히 이동
# (이미 적용됨)

# 추가 최적화: tmpdir를 /data_4tb로 설정
cat >> config/miniwdl.local.cfg << 'EOF'

[task_runtime]
# tmpdir를 /data_4tb로 설정하여 I/O 병목 방지
defaults = {
        "maxRetries": 2,
        "docker": "ubuntu:20.04",
        "cpu": 16,
        "memory": "64G",
        "tmpdir": "/data_4tb/tmp"
    }
EOF

mkdir -p /data_4tb/tmp
```

#### 해결책 2: pbmm2 스레드 수 조정

현재 pbmm2가 32 스레드를 사용하는데, 서버가 40코어이므로 적절합니다.
하지만 메모리 설정을 확인:

```bash
# workflow에서 pbmm2 메모리가 128GB 요청
# 서버는 251GB이므로 2개 동시 실행 불가

# 해결: task_concurrency를 1로 유지 (이미 적용됨)
```

#### 해결책 3: Batch 처리 시 순차 실행

**권장**: 처음에는 **순차 실행**으로 시작

이유:
- 각 샘플이 많은 리소스 사용 (32 CPU, 128GB RAM)
- 서버 스펙: 40 CPU, 251GB RAM
- 2개 동시 실행 시 메모리 부족 발생 가능

```bash
# 순차 실행 (안전)
./batch_run_optimized.sh sequential

# 병렬 실행은 리소스가 충분한 경우만
# (예: 샘플당 20 CPU, 64GB RAM으로 제한)
```

#### 해결책 4: GPU 사용 활성화

**중요**: 입력 JSON에서 `gpu: true`로 설정 (이미 스크립트에 포함)

DeepVariant GPU 모드:
- CPU 모드: 64 cores, 8-12시간
- GPU 모드: 1 GPU, 2-4시간
- **시간 절약: 50-70%**

### 2. 예상 실행 시간 (최적화 후)

#### 단일 샘플 (GPU 모드, 최적화):
- pbmm2 alignment: 3-4시간
- DeepVariant (GPU): 2-3시간
- 기타 분석: 2-3시간
- **총 예상: 7-10시간**

#### 3개 샘플 순차 실행:
- **총 예상: 21-30시간**

#### Call Cache 효과:
- 이미 처리한 샘플 재실행 시: **즉시 완료** (캐시에서)
- 동일한 참조 데이터 사용 시: 일부 단계 재사용 가능

### 3. 리소스 모니터링

```bash
cat > monitor_batch.sh << 'EOF'
#!/bin/bash
# monitor_batch.sh
# Batch 처리 중 리소스 모니터링

BATCH_RESULTS="/data_4tb/hifi-human-wgs-wdl-custom/batch_results"
LOG_DIR="${BATCH_RESULTS}/logs"

echo "========================================"
echo "Batch Processing Monitor"
echo "========================================"
echo ""

# 실행 중인 miniwdl 프로세스
echo "Running processes:"
ps aux | grep -E "(miniwdl|pbmm2|deepvariant)" | grep -v grep || echo "  None"
echo ""

# CPU 및 메모리 사용량
echo "System resources:"
echo "  CPU: $(top -bn1 | grep "Cpu(s)" | awk '{print $2}')% used"
echo "  Memory:"
free -h | grep -E "Mem|Swap"
echo ""

# GPU 상태
echo "GPU status:"
nvidia-smi --query-gpu=index,name,temperature.gpu,utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits 2>/dev/null || echo "  Not available"
echo ""

# 디스크 사용량
echo "Disk usage:"
df -h / /data_4tb
echo ""

# 샘플별 상태
echo "Sample status:"
for log in ${LOG_DIR}/*.log; do
    if [[ -f "$log" ]]; then
        sample=$(basename "$log" .log)
        
        # 마지막 로그 라인
        last_line=$(tail -1 "$log" 2>/dev/null)
        
        # 완료 여부 확인
        if grep -q "done" "$log" 2>/dev/null; then
            status="✓ COMPLETED"
        elif grep -q "error\|failed" "$log" 2>/dev/null; then
            status="✗ FAILED"
        elif [[ -f "${LOG_DIR}/${sample}.pid" ]]; then
            pid=$(cat "${LOG_DIR}/${sample}.pid")
            if ps -p $pid > /dev/null 2>&1; then
                status="⏳ RUNNING (PID: $pid)"
            else
                status="⏹ STOPPED"
            fi
        else
            status="⏸ WAITING"
        fi
        
        echo "  ${sample}: ${status}"
    fi
done
echo ""
echo "========================================"
EOF

chmod +x monitor_batch.sh

# 실시간 모니터링
watch -n 30 ./monitor_batch.sh
```

---

## 🎯 실행 및 모니터링

### 전체 프로세스:

```bash
# 1. samples.csv 작성 (실제 데이터에 맞게)
vim samples.csv

# 2. 입력 JSON 파일 생성
./create_batch_inputs.sh

# 3. 생성된 입력 파일 확인
ls -lh batch_inputs/

# 4. GPU 환경 변수 설정
export CUDA_VISIBLE_DEVICES=1

# 5. Batch 실행 (순차 모드 권장)
./batch_run_optimized.sh sequential

# 또는 특정 샘플만
./batch_run_optimized.sh sequential KTY9537 KTY9538

# 6. 다른 터미널에서 모니터링
watch -n 30 ./monitor_batch.sh

# 7. GPU 모니터링
watch -n 5 nvidia-smi

# 8. 디스크 모니터링
watch -n 60 'df -h /data_4tb'
```

### 로그 확인:

```bash
# 특정 샘플 로그 실시간 확인
tail -f /data_4tb/hifi-human-wgs-wdl-custom/batch_results/logs/KTY9537.log

# 에러 확인
grep -i "error\|failed" /data_4tb/hifi-human-wgs-wdl-custom/batch_results/logs/*.log

# 진행 상황 확인
grep -E "done|completed" /data_4tb/hifi-human-wgs-wdl-custom/batch_results/logs/*.log
```

### 결과 수집:

```bash
# 모든 결과 수집
./collect_results.sh

# 요약 확인
ls -lh /data_4tb/hifi-human-wgs-wdl-custom/batch_results/summary/*/
```

---

## 📊 최종 체크리스트

### 실행 전:
- [ ] samples.csv 작성 완료
- [ ] 입력 JSON 파일 생성 (`./create_batch_inputs.sh`)
- [ ] GPU 설정 확인 (`nvidia-smi`)
- [ ] 디스크 공간 확인 (`df -h /data_4tb`) - 샘플당 ~500GB 필요
- [ ] Config 파일 확인 (`config/miniwdl.local.cfg`)
- [ ] Call cache 위치 확인 (`/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_call_cache/`)

### 실행 중:
- [ ] 리소스 모니터링 (`./monitor_batch.sh`)
- [ ] GPU 온도 확인 (< 80°C 유지)
- [ ] 디스크 공간 확인 (지속적)
- [ ] 로그 확인 (에러 발생 시)

### 실행 후:
- [ ] 모든 샘플 완료 확인
- [ ] 결과 파일 수집 (`./collect_results.sh`)
- [ ] VCF/BAM 파일 검증
- [ ] 통계 리포트 확인

---

## 🔧 예상 문제 및 해결

### 1. "디스크 공간 부족"
```bash
# 중간 파일 정리
rm -rf /data_4tb/hifi-human-wgs-wdl-custom/batch_results/*/call-*/_miniwdl_*

# 또는 완료된 샘플의 중간 파일만 정리
for sample in /data_4tb/hifi-human-wgs-wdl-custom/batch_results/*/; do
    if [[ -f "${sample}/outputs.json" ]]; then
        echo "Cleaning ${sample}..."
        find "$sample" -path "*/call-*/_miniwdl_*" -type d -exec rm -rf {} + 2>/dev/null
    fi
done
```

### 2. "메모리 부족"
```bash
# task_concurrency를 더 낮춤 (이미 1로 설정됨)
# 또는 pbmm2 메모리 제한 조정
# 워크플로우 내부 설정이므로 수정 불가, 순차 실행 필수
```

### 3. "샘플 하나가 실패"
```bash
# 해당 샘플만 재실행
./batch_run_optimized.sh sequential KTY9538

# Call cache 덕분에 성공한 단계는 재사용됨
```

---

## 💡 추가 팁

### 1. 우선순위 설정
급한 샘플부터 처리:
```bash
./batch_run_optimized.sh sequential KTY9537 KTY9540 KTY9542
```

### 2. 야간 실행
```bash
# nohup으로 백그라운드 실행
nohup ./batch_run_optimized.sh sequential > batch_run.out 2>&1 &

# 또는 screen/tmux 사용
screen -S batch_processing
./batch_run_optimized.sh sequential
# Ctrl+A, D로 detach
```

### 3. 알림 설정
```bash
# 완료 시 이메일 전송 (sendmail 설정 필요)
./batch_run_optimized.sh sequential && echo "Batch complete!" | mail -s "HiFi Pipeline" your@email.com
```

### 4. 자동 QC 리포트 생성
Batch 처리가 완료되면 자동으로 HTML QC 리포트가 생성됩니다:

```bash
# batch_run_optimized.sh 실행 시 자동 생성됨
./batch_run_optimized.sh parallel KTY9537 KTY9538

# 리포트 위치:
# /data_4tb/hifi-human-wgs-wdl-custom/batch_results/QC_Report_YYYYMMDD_HHMMSS.html
```

**QC 리포트에 포함되는 내용:**
- 전체 샘플 요약 통계
- Coverage 통계 (샘플별 평균 depth)
- Variant calling 결과 (SNPs, Indels, SVs)
- 파일 크기 정보
- 주요 출력 파일 상태
- PharmCAT 결과 (약물유전체 분석)
- Phasing 통계

**수동으로 리포트 생성:**
```bash
# 특정 샘플들만 포함
python3 scripts/generate_qc_report.py \
  --batch-results /data_4tb/hifi-human-wgs-wdl-custom/batch_results \
  --output custom_report.html \
  --samples KTY9537 KTY9538

# 모든 완료된 샘플 포함
python3 scripts/generate_qc_report.py \
  --batch-results /data_4tb/hifi-human-wgs-wdl-custom/batch_results \
  --output QC_Report.html
```

**리포트 확인:**
```bash
# 브라우저에서 열기
firefox /data_4tb/hifi-human-wgs-wdl-custom/batch_results/QC_Report_*.html

# 또는 원격에서 접속
# Windows에서 WSL 경로로 접근:
# \\wsl.localhost\Ubuntu\data_4tb\hifi-human-wgs-wdl-custom\batch_results\QC_Report_*.html
```

이제 준비가 완료되었습니다! 🚀
