#!/usr/bin/env bash
# 원본 HiFi BAM에서 샘플마다 앞쪽 N개 read를 뽑아 마우스+마이코플라스마 패널(7종)에 정렬하고,
# read별 정렬 대상과 6mA 호출 밀도를 기록한다 (6mA 신호의 오염 기원 점검).
# 패널: A. laidlawii, M. hominis, M. hyorhinis, M. fermentans, M. arginini, M. salivarium, M. orale
# 사용: bash mycoplasma_check.sh [N_READS] > <OUT>/mycoplasma_check.log 2>&1
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
N="${1:-100000}"
THREADS="${THREADS:-4}"
OUT="${OUT:-/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/mycoplasma}"
RES=/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39
MMI="${RES}/mouse_GRCm39_plus_myco.map-hifi.mmi"
CACHE=/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_singularity_cache
SAMTOOLS=(apptainer exec --bind /mnt,/data_4tb "${CACHE}/docker___quay.io_biocontainers_samtools_1.21--h50ea8bc_0.sif" samtools)
MINIMAP2=(apptainer exec --bind /mnt,/data_4tb "${CACHE}/docker___quay.io_biocontainers_minimap2_2.30--h577a1d6_0.sif" minimap2)
PY=/home/ygkim/program/anaconda3/bin/python

mkdir -p "${OUT}"
tail -n +2 "${HERE}/samples.tsv" | while IFS=$'\t' read -r sid run group bam raw; do
  res="${OUT}/${sid}.reads.tsv.gz"
  [ -f "${res}" ] && { echo "[$(date '+%F %T')] skip ${sid}"; continue; }
  echo "[$(date '+%F %T')] ${sid}: first ${N} reads of ${raw}"
  # head가 파이프를 먼저 닫으므로 samtools의 SIGPIPE는 무시
  { "${SAMTOOLS[@]}" view -h "${raw}" || true; } | awk -v n="${N}" '/^@/{print; next} {print; if (++c >= n) exit}' \
    | "${SAMTOOLS[@]}" fastq -T MM,ML - 2>/dev/null \
    | "${MINIMAP2[@]}" -a -y -x map-hifi --secondary=no -t "${THREADS}" "${MMI}" - 2>/dev/null \
    | "${PY}" "${HERE}/mycoplasma_classify.py" "${res}"
done
"${PY}" "${HERE}/mycoplasma_classify.py" --summary "${OUT}"
