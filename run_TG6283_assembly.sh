#!/bin/bash
# TG6283 de novo assembly 분석 스크립트
# reference-based 분석 완료 후 assembly tasks만 추가 실행
# call cache로 pbmm2/deepvariant 등 reference tasks는 재실행하지 않음

set -o pipefail

source /home/ygkim/program/anaconda3/etc/profile.d/conda.sh
conda activate hifi-human-wgs

WORK_DIR="/data_4tb/hifi-human-wgs-wdl-custom"
INPUT="TG6283.assembly.inputs.json"
OUTPUT_DIR="/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6283_B01"
LOG_FILE="/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6283_B01.assembly.log"

echo "========================================"
echo "TG6283 de novo assembly 분석"
echo "시작: $(date)"
echo "로그: ${LOG_FILE}"
echo "========================================"

cd "${WORK_DIR}" || exit 1

miniwdl run workflows/singleton.wdl \
    --input "${INPUT}" \
    --dir "${OUTPUT_DIR}" \
    --verbose \
    2>&1 | tee "${LOG_FILE}"

exit_code=${PIPESTATUS[0]}

echo ""
echo "========================================"
if [[ ${exit_code} -eq 0 ]]; then
    echo "완료: $(date)"
else
    echo "실패 (exit code: ${exit_code}): $(date)"
    echo "로그 확인: ${LOG_FILE}"
fi
echo "========================================"

exit ${exit_code}
