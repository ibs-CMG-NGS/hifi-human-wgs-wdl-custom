#!/bin/bash
# Base Modification Discovery — kinetics 기반 5mC/5hmC/6mA + unknown modification 탐색
# 사용법: bash run_base_mod_discovery.sh [inputs_json] [output_dir]
#
# 예시 (파일럿: 그룹당 1개 BAM):
#   bash run_base_mod_discovery.sh base_mod_discovery.pilot.inputs.json \
#     /data_4tb/hifi-human-wgs-wdl-custom/output/base_mod_discovery_pilot

set -o pipefail

source /home/ygkim/program/anaconda3/etc/profile.d/conda.sh
conda activate hifi-human-wgs

WORK_DIR="/data_4tb/hifi-human-wgs-wdl-custom"
INPUT_JSON="${1:-base_mod_discovery.inputs.json}"
OUTPUT_DIR="${2:-/data_4tb/hifi-human-wgs-wdl-custom/output/base_mod_discovery}"
LOG_FILE="${OUTPUT_DIR}.nohup.log"

echo "========================================"
echo "Base Modification Discovery Pipeline"
echo "Output: ${OUTPUT_DIR}"
echo "시작: $(date)"
echo "========================================"

cd "${WORK_DIR}" || exit 1

if [[ ! -f "${INPUT_JSON}" ]]; then
    echo "ERROR: ${INPUT_JSON} not found."
    echo "base_mod_discovery.inputs.template.json 을 복사하여 /mnt/ BAM 경로를 채운 뒤 다시 실행하세요."
    exit 1
fi

miniwdl run workflows/base_modification_discovery.wdl \
    --input "${INPUT_JSON}" \
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
