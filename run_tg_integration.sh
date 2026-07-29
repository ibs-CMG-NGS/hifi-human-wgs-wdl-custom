#!/bin/bash
# Transgene Integration Analysis 실행 스크립트
# 사용법: bash run_tg_integration.sh <sample_id> <inputs_json> <output_dir>
#
# 예시:
#   bash run_tg_integration.sh TG-6283 TG6283.tg_integration.inputs.json \
#     /mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6283_integration_wdl

set -o pipefail

source /home/ygkim/program/anaconda3/etc/profile.d/conda.sh
conda activate hifi-human-wgs

SAMPLE_ID="${1:-TG-6283}"
INPUT="${2:-TG6283.tg_integration.inputs.json}"
OUTPUT_DIR="${3:-/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/${SAMPLE_ID}_integration_wdl}"
LOG_FILE="${OUTPUT_DIR}.tg_integration.log"

WORK_DIR="/data_4tb/hifi-human-wgs-wdl-custom"

echo "========================================"
echo "Transgene Integration Analysis"
echo "Sample:  ${SAMPLE_ID}"
echo "Input:   ${INPUT}"
echo "Output:  ${OUTPUT_DIR}"
echo "Log:     ${LOG_FILE}"
echo "시작: $(date)"
echo "========================================"

cd "${WORK_DIR}" || exit 1

miniwdl run workflows/transgene_integration.wdl \
    --input "${INPUT}" \
    --dir "${OUTPUT_DIR}" \
    --verbose \
    < /dev/null 2>&1 | tee "${LOG_FILE}"

exit_code=${PIPESTATUS[0]}

echo ""
echo "========================================"
if [[ ${exit_code} -eq 0 ]]; then
    echo "완료: $(date)"
    echo "리포트: ${OUTPUT_DIR}/_LAST/out/report_txt/"
else
    echo "실패 (exit code: ${exit_code}): $(date)"
    echo "로그: ${LOG_FILE}"
fi
echo "========================================"

exit ${exit_code}
