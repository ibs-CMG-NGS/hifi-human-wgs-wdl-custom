#!/usr/bin/env bash
# 샘플 하나의 상염색체 PASS 이형접합 indel 수와 결실/삽입 구성
# 사용: het_indel_count.sh SAMPLE_ID VCF
BCF=/home/ygkim/program/anaconda3/envs/wgs-tertiary-pipeline/bin/bcftools
AUT=$(printf 'chr%s,' $(seq 1 19)); AUT=${AUT%,}
"$BCF" query -t "$AUT" -i 'FILTER="PASS"' -f '%REF\t%ALT[\t%GT]\n' "$2" | awk -F'\t' -v s="$1" '
  { if (length($1)==1 && length($2)==1) next
    split($3, a, /[|\/]/); if (a[1]==a[2]) next
    n++; if (length($1) > length($2)) del++; else ins++
    if (length($1) - length($2) == 1) del1++ }
  END { print s "\t" n "\t" del "\t" ins "\t" del1 }'
