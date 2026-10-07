#!/usr/bin/env bash
# run_6mA_pileup.sh가 끝낸(.done) 샘플을 차례로 aggregate_6mA.py로 집계 (1코어).
# 12샘플이 모두 집계되면 종료. 사용: setsid nohup bash run_6mA_aggregate.sh > <DIR>/run_6mA_aggregate.log 2>&1 < /dev/null &
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DIR="${DIR:-/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/6mA}"
PY=/home/ygkim/program/anaconda3/bin/python
mapfile -t SIDS < <(tail -n +2 "${HERE}/samples.tsv" | cut -f1)
while :; do
  left=0
  for sid in "${SIDS[@]}"; do
    agg="${DIR}/${sid}.6mA.agg.tsv.gz"
    [ -f "${agg}" ] && continue
    left=$((left + 1))
    [ -f "${DIR}/${sid}.6mA.tsv.gz.done" ] || continue
    echo "[$(date '+%F %T')] aggregate ${sid}"
    if "${PY}" "${HERE}/aggregate_6mA.py" "${DIR}/${sid}.6mA.tsv.gz" "${agg}.tmp.gz" > "${DIR}/${sid}.agg.log" 2>&1; then
      mv "${agg}.tmp.gz" "${agg}"
      echo "[$(date '+%F %T')] done ${sid}"
    else
      echo "[$(date '+%F %T')] FAIL ${sid} — see ${DIR}/${sid}.agg.log"; exit 1
    fi
  done
  [ "${left}" -eq 0 ] && break
  sleep 600
done
echo "[$(date '+%F %T')] all aggregated"
