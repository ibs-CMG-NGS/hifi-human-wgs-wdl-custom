#!/bin/bash
# TG6283 분석 재개 스크립트
# 중단 지점: call-upstream > call-deepvariant > deepvariant_call_variants_gpu
# 캐시 히트 예상: deepvariant_make_examples (8개 샤드) → call_variants_gpu부터 재실행

set -o pipefail

WORK_DIR="/data_4tb/hifi-human-wgs-wdl-custom"
INPUT="TG6283.inputs.json"
OUTPUT_DIR="/media/ygkim/e847f448-7179-490b-89b1-19eb760a361e/tg-integration-analysis/TG6283"
LOG_FILE="/media/ygkim/e847f448-7179-490b-89b1-19eb760a361e/tg-integration-analysis/TG6283.rerun.log"

echo "========================================"
echo "TG6283 재분석 재개"
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
