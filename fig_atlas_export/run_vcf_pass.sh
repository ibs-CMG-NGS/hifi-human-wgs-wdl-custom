#!/usr/bin/env bash
# 12샘플 vcf_sample_pass.py를 4개씩 병렬 실행
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "${HERE}/.." && pwd)"
OUTD="$1"
QC=/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/figure_data/sequencing_qc.tsv
PY=/home/ygkim/program/anaconda3/bin/python
tail -n +2 "${REPO}/base_mod_per_sample/samples.tsv" | while IFS=$'\t' read -r sid run group bam raw; do
  o="${bam%%/out/*}/out"
  # sequencing_qc.tsv 8번째 열 = depth_autosomal
  dp=$(awk -F'\t' -v s="$sid" 'NR==1{for(i=1;i<=NF;i++) if($i=="depth_autosomal") c=i; next} $1==s{print $c}' "$QC")
  echo "$PY ${HERE}/vcf_sample_pass.py $sid $o/phased_small_variant_vcf/${sid}.GRCm39.small_variants.phased.vcf.gz $o/phased_sv_vcf/${sid}.GRCm39.structural_variants.phased.vcf.gz $dp ${OUTD}/${sid}"
done | xargs -P 4 -I{} bash -c "{}"
