#!/bin/bash
# TG-6425 C01 — humanwgs_singleton (reference-based + de novo assembly)
# 실행: bash run_TG6425_C01.sh

set -o pipefail

source /home/ygkim/program/anaconda3/etc/profile.d/conda.sh
conda activate hifi-human-wgs

WORK_DIR="/data_4tb/hifi-human-wgs-wdl-custom"
OUTPUT_DIR="/mnt/Ext_1tb_A/tg-integration-denovo-assembly/TG6425_C01"
LOG_FILE="${OUTPUT_DIR}.nohup.log"

echo "========================================"
echo "TG-6425 C01 humanwgs_singleton"
echo "Output: ${OUTPUT_DIR}"
echo "시작: $(date)"
echo "========================================"

cd "${WORK_DIR}" || exit 1

miniwdl run workflows/singleton.wdl \
    --input TG6425.inputs.json \
    --dir "${OUTPUT_DIR}" \
    --verbose \
    2>&1 | tee "${LOG_FILE}"

exit_code=${PIPESTATUS[0]}

echo ""
echo "========================================"
if [[ ${exit_code} -eq 0 ]]; then
    echo "완료: $(date)"
    echo "결과: ${OUTPUT_DIR}/_LAST/out/"
else
    echo "실패 (exit code: ${exit_code}): $(date)"
fi
echo "========================================"
exit ${exit_code}
