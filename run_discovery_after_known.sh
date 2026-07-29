#!/bin/bash
# A(known-mod) 완료를 기다린 뒤 B(전체 discovery: known + unknown) 자동 실행
set -o pipefail
source /home/ygkim/program/anaconda3/etc/profile.d/conda.sh
conda activate hifi-human-wgs
cd /data_4tb/hifi-human-wgs-wdl-custom

A_PID="${1:-314238}"

echo "[$(date)] A(known) PID ${A_PID} 완료 대기 중..."
while kill -0 "${A_PID}" 2>/dev/null; do
  sleep 60
done
echo "[$(date)] A 종료 감지. B(discovery) 시작."

miniwdl run workflows/base_modification_discovery.wdl \
  --input base_mod_discovery.pilot.inputs.json \
  --dir /mnt/JJ_dis_8tb/base_mod_discovery_pilot_output/20260706_discovery \
  --verbose 2>&1 | tee output/base_mod_discovery_full.nohup.log

echo "[$(date)] B(discovery) 종료 (exit ${PIPESTATUS[0]})"
