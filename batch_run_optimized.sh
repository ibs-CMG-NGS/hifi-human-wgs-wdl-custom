#!/bin/bash
# batch_run_optimized.sh
# /data_4tb에 출력하는 최적화된 배치 처리 스크립트

# 사용법:
# ./batch_run_optimized.sh [parallel|sequential] [sample1 sample2 ...]
# 
# 예시:
# ./batch_run_optimized.sh parallel              # 모든 샘플 병렬 실행
# ./batch_run_optimized.sh sequential KTY9537 KTY9538  # 특정 샘플만 순차 실행

# set -e 제거 - 수동으로 에러 처리
set -o pipefail  # 파이프에서 에러 감지

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
        
        # 파이프 처리 개선: set +e로 에러 무시, 수동으로 exit code 캡처
        set +e
        miniwdl run "${WORKFLOW}" \
            --input "${input_file}" \
            --cfg "${CONFIG_FILE}" \
            --dir "${output_dir}" \
            --verbose \
            2>&1 | tee "${log_file}"
        
        exit_code=${PIPESTATUS[0]}
        set -e
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
