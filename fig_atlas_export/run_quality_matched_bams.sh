#!/usr/bin/env bash
# 12샘플의 chr1–3 품질 맞춘 BAM 생성 (낮은 우선순위, 동시 3개)
# 사용: run_quality_matched_bams.sh OUT_DIR
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "${HERE}/.." && pwd)"
OUT="$1"
KEEP=/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/fig_atlas_work/qmatch/keep_prob_bins.tsv
SIF="${REPO}/miniwdl_singularity_cache/docker___quay.io_biocontainers_samtools_1.21--h50ea8bc_0.sif"
PY=/home/ygkim/program/anaconda3/bin/python
mkdir -p "$OUT"
one() {
  sid="$1"; bam="$2"; o="$OUT/${sid}.qm.chr1-3.bam"
  [ -f "$o.bai" ] && { echo "skip $sid"; return; }
  nice -n 10 apptainer exec --bind /mnt,/data_4tb "$SIF" samtools view -h "$bam" chr1 chr2 chr3 \
    | nice -n 10 "$PY" "$HERE/quality_matched_subsample.py" "$KEEP" "$sid" \
    | nice -n 10 apptainer exec --bind /mnt,/data_4tb "$SIF" samtools view -b -@2 -o "$o.tmp" -
  mv "$o.tmp" "$o"
  apptainer exec --bind /mnt,/data_4tb "$SIF" samtools index "$o"
  echo "[$(date '+%F %T')] done $sid"
}
export -f one; export OUT KEEP SIF PY HERE
tail -n +2 "${REPO}/base_mod_per_sample/samples.tsv" | cut -f1,4 | xargs -P 3 -n 2 bash -c 'one "$0" "$1"'
