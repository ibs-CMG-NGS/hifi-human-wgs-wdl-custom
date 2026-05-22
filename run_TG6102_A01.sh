#!/bin/bash
# TG-6102 A01 humanwgs_singleton 파이프라인 실행 스크립트
# reference-based mapping은 call cache에서 가져오고,
# run_assembly: true로 de novo assembly만 추가 실행
#
# 실행 방법:
#   nohup bash run_TG6102_A01.sh > /mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6102_A01.nohup.log 2>&1 &

set -o pipefail

source /home/ygkim/program/anaconda3/etc/profile.d/conda.sh
conda activate hifi-human-wgs

WORK_DIR="/data_4tb/hifi-human-wgs-wdl-custom"
INPUT="TG6102.inputs.json"
OUTPUT_DIR="/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6102_A01"
LOG_FILE="/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6102_A01.assembly.log"

echo "========================================"
echo "TG-6102 A01 humanwgs_singleton 분석"
echo "  (ref 태스크: call cache 재사용, assembly만 신규 실행)"
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
    echo ""
    echo "Assembly 결과 확인:"
    ls -lh "${OUTPUT_DIR}/_LAST/out/assembly_hap1_fasta/" 2>/dev/null || echo "  hap1_fasta 없음"
    ls -lh "${OUTPUT_DIR}/_LAST/out/assembly_hap2_fasta/" 2>/dev/null || echo "  hap2_fasta 없음"
    echo ""
    echo "다음 단계: transgene integration 분석"
    echo "  bash run_tg_integration.sh TG-6102 TG6102.tg_integration.inputs.json \\"
    echo "    /mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6102_A01_integration_analysis"
else
    echo "실패 (exit code: ${exit_code}): $(date)"
    echo "로그 확인: ${LOG_FILE}"
fi
echo "========================================"

exit ${exit_code}
