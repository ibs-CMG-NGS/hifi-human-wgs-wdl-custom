#!/usr/bin/env bash
# 12개 H2O2 샘플의 WGS haplotagged BAM에서 6mA per-sample pileup 생성.
# WGS BAM에는 장비 내 jasmine 2.3.0이 호출한 6mA(A+a./T-a.)가 이미 있어 재호출 불필요.
# - 6mA는 확률 ≥0.94 호출만 기록된 implicit 형식이라 modkit 자동 threshold 추정이 실패함
#   → 전 샘플 공통 고정 threshold(0.7) 사용 (샘플 간 비교 공정성)
# - 출력은 DSS 입력용 5컬럼(chr, pos, strand, N_valid_cov, N_mod)만 bgzip (샘플당 ~6G)
# 사용: nohup bash run_6mA_pileup.sh > run_6mA_pileup.log 2>&1 &
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
SAMPLES="${HERE}/samples.tsv"
OUTDIR="${OUTDIR:-/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/6mA}"
REF="/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/mouse_GRCm39.fasta"
SIF="/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_singularity_cache/docker___ontresearch_modkit_latest.sif"
THREADS="${THREADS:-8}"
MIN_COV="${MIN_COV:-5}"
FILTER_THRESHOLD="0.7"

mkdir -p "${OUTDIR}"
EXEC=(apptainer exec --bind /mnt,/data_4tb "${SIF}")

tail -n +2 "${SAMPLES}" | while IFS=$'\t' read -r sid run group bam raw; do
  out="${OUTDIR}/${sid}.6mA.tsv.gz"
  if [ -f "${out}.done" ]; then
    echo "[$(date '+%F %T')] skip ${sid} (done)"
    continue
  fi
  [ -f "${bam}" ] || { echo "[$(date '+%F %T')] ERROR missing BAM ${bam}"; exit 1; }
  echo "[$(date '+%F %T')] start ${sid} (run ${run}, group ${group})"

  # columns: 1 chrom, 2 start, 4 mod_code, 6 strand, 10 N_valid_cov, 12 N_mod
  "${EXEC[@]}" modkit pileup "${bam}" - \
      --reference "${REF}" \
      --threads "${THREADS}" \
      --filter-threshold "${FILTER_THRESHOLD}" \
      --log-filepath "${OUTDIR}/${sid}.modkit.log" \
    | awk -F'\t' -v OFS='\t' -v mc="${MIN_COV}" '$4=="a" && $10>=mc {print $1,$2,$6,$10,$12}' \
    | "${EXEC[@]}" bgzip -c > "${out}.tmp"

  mv "${out}.tmp" "${out}"
  n=$("${EXEC[@]}" bgzip -dc "${out}" | wc -l)
  echo "${n}" > "${out}.done"
  echo "[$(date '+%F %T')] done ${sid}: ${n} sites, $(du -h "${out}" | cut -f1)"
done

echo "[$(date '+%F %T')] all samples finished"
