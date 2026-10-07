#!/usr/bin/env bash
# Fig 2e 입력 생성: 12샘플 DeepVariant VCF(PASS+RefCall)의 상염색체 biallelic SNV를 병합해
# mutation_spectrum.py로 샘플 고유 SNV 스펙트럼 계산.
# 사용: bash mutation_spectrum.sh OUT_PREFIX
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$1"
BCF=/home/ygkim/program/anaconda3/envs/wgs-tertiary-pipeline/bin/bcftools
PY=/home/ygkim/program/anaconda3/bin/python
TMP="${OUT}.tmp"; mkdir -p "${TMP}"
REG=$(printf 'chr%s,' $(seq 1 19)); REG=${REG%,}
list=()
while IFS=$'\t' read -r sid run group bam raw; do
  vcf="${bam%%/out/*}/out/phased_small_variant_vcf/${sid}.GRCm39.small_variants.phased.vcf.gz"
  o="${TMP}/${sid}.snv.bcf"
  if [ ! -f "${o}.csi" ]; then
    # 원본 VCF 옆에 인덱스가 없으므로 -t(스트리밍 필터) 사용
    "${BCF}" view -t "${REG}" -v snps -m2 -M2 -Ob -o "${o}" "${vcf}" || { echo "FAIL ${sid}"; exit 1; }
    "${BCF}" index "${o}" || { echo "FAIL index ${sid}"; exit 1; }
  fi
  list+=("${o}")
  echo "[$(date '+%T')] ${sid} extracted"
done < <(tail -n +2 "${HERE}/samples.tsv")
"${BCF}" merge -m none --force-samples "${list[@]}" -Ou \
  | "${BCF}" query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%AD\t%DP]\n' \
  | "${PY}" "${HERE}/mutation_spectrum.py" "${HERE}/samples.tsv" "${OUT}"
