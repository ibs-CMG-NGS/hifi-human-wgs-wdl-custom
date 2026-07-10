#!/bin/bash
# v2 재시작: compare_known_mods를 sort-merge 스트리밍 + 6mA 배경(0%) 필터링 방식으로
# 재작성한 뒤 재실행. 업스트림(align/merge/filter/pileup/jasmine/cpg_qc)은 call cache로 재사용.
set -o pipefail
source /home/ygkim/program/anaconda3/etc/profile.d/conda.sh
conda activate hifi-human-wgs
cd /data_4tb/hifi-human-wgs-wdl-custom

miniwdl run workflows/base_modification_known.wdl \
  --input base_mod_known.pilot.inputs.json \
  --dir pipeline_runs/20260708_known_v2 \
  --verbose 2>&1 | tee output/base_mod_known_v3.nohup.log

echo "[$(date)] base_modification_known (v3) 종료 (exit ${PIPESTATUS[0]})"
