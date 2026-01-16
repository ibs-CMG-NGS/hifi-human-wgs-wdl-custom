#!/bin/bash
# batch_run.sh
# 여러 샘플을 순차적으로 또는 병렬로 처리하는 배치 스크립트

# 사용법:
# ./batch_run.sh [parallel|sequential]
# 
# 예시:
# ./batch_run.sh parallel    # 샘플들을 백그라운드에서 병렬 실행
# ./batch_run.sh sequential  # 샘플들을 순차적으로 실행

# 설정
WORKFLOW="workflows/singleton.wdl"
INPUT_DIR="inputs"
OUTPUT_DIR="outputs"
LOG_DIR="logs"

# 처리할 샘플 목록
SAMPLES=("sample1" "sample2" "sample3" "sample4")

# 실행 모드 (기본값: parallel)
MODE="${1:-parallel}"

# 디렉토리 생성
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

echo "========================================"
echo "HiFi-human-WGS Batch Processing Script"
echo "========================================"
echo "Workflow: ${WORKFLOW}"
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Log directory: ${LOG_DIR}"
echo "Mode: ${MODE}"
echo "Number of samples: ${#SAMPLES[@]}"
echo "========================================"
echo ""

# Conda 환경 활성화 확인
if [[ -z "${CONDA_DEFAULT_ENV}" ]] || [[ "${CONDA_DEFAULT_ENV}" != "hifi-human-wgs" ]]; then
    echo "Warning: hifi-human-wgs conda environment is not activated."
    echo "Please run: conda activate hifi-human-wgs"
    exit 1
fi

# miniwdl 설치 확인
if ! command -v miniwdl &> /dev/null; then
    echo "Error: miniwdl is not installed or not in PATH"
    exit 1
fi

# 병렬 실행 함수
run_parallel() {
    echo "Running samples in parallel mode..."
    echo ""
    
    for sample in "${SAMPLES[@]}"; do
        input_file="${INPUT_DIR}/${sample}.inputs.json"
        
        # 입력 파일 존재 확인
        if [[ ! -f "${input_file}" ]]; then
            echo "Warning: Input file not found: ${input_file}"
            echo "Skipping ${sample}..."
            echo ""
            continue
        fi
        
        echo "Starting ${sample}..."
        miniwdl run "${WORKFLOW}" \
            --input "${input_file}" \
            --dir "${OUTPUT_DIR}/${sample}" \
            > "${LOG_DIR}/${sample}.log" 2>&1 &
        
        # 프로세스 ID 저장
        echo "  PID: $!"
        echo "  Input: ${input_file}"
        echo "  Output: ${OUTPUT_DIR}/${sample}"
        echo "  Log: ${LOG_DIR}/${sample}.log"
        echo ""
    done
    
    echo "All samples started. Waiting for completion..."
    wait
    echo ""
    echo "All samples completed!"
}

# 순차 실행 함수
run_sequential() {
    echo "Running samples in sequential mode..."
    echo ""
    
    for sample in "${SAMPLES[@]}"; do
        input_file="${INPUT_DIR}/${sample}.inputs.json"
        
        # 입력 파일 존재 확인
        if [[ ! -f "${input_file}" ]]; then
            echo "Warning: Input file not found: ${input_file}"
            echo "Skipping ${sample}..."
            echo ""
            continue
        fi
        
        echo "========================================"
        echo "Processing ${sample}..."
        echo "========================================"
        echo "  Input: ${input_file}"
        echo "  Output: ${OUTPUT_DIR}/${sample}"
        echo "  Log: ${LOG_DIR}/${sample}.log"
        echo ""
        
        miniwdl run "${WORKFLOW}" \
            --input "${input_file}" \
            --dir "${OUTPUT_DIR}/${sample}" \
            2>&1 | tee "${LOG_DIR}/${sample}.log"
        
        exit_code=${PIPESTATUS[0]}
        
        if [[ ${exit_code} -eq 0 ]]; then
            echo ""
            echo "✓ ${sample} completed successfully"
            echo ""
        else
            echo ""
            echo "✗ ${sample} failed with exit code ${exit_code}"
            echo "Check log file: ${LOG_DIR}/${sample}.log"
            echo ""
            
            # 순차 실행 중 실패 시 계속할지 물어봄
            read -p "Continue with next sample? (y/n) " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                echo "Stopping batch processing."
                exit ${exit_code}
            fi
        fi
    done
    
    echo "========================================"
    echo "All samples processed!"
    echo "========================================"
}

# 실행 시작 시간 기록
start_time=$(date +%s)

# 모드에 따라 실행
case "${MODE}" in
    parallel)
        run_parallel
        ;;
    sequential)
        run_sequential
        ;;
    *)
        echo "Error: Unknown mode '${MODE}'"
        echo "Usage: $0 [parallel|sequential]"
        exit 1
        ;;
esac

# 실행 종료 시간 계산
end_time=$(date +%s)
elapsed=$((end_time - start_time))
hours=$((elapsed / 3600))
minutes=$(((elapsed % 3600) / 60))
seconds=$((elapsed % 60))

echo ""
echo "========================================"
echo "Batch processing summary"
echo "========================================"
echo "Total samples: ${#SAMPLES[@]}"
echo "Elapsed time: ${hours}h ${minutes}m ${seconds}s"
echo "Output directory: ${OUTPUT_DIR}"
echo "Log directory: ${LOG_DIR}"
echo ""
echo "To check the status of each sample, see:"
echo "  ${LOG_DIR}/*.log"
echo "========================================"
