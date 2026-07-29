#!/bin/bash
# Integration site 주변 SV 기능적 어노테이션 분석
# 사용법: bash sv_functional_annotation.sh
set -euo pipefail

ANNO=/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39
GTF=$ANNO/gencode.vM36.annotation.sorted.gtf.bgz
CPG=$ANNO/cpgIslandExt.sorted.mm39.tsv
REPEAT=$ANNO/simpleRepeat.txt.gz
OUT=/data_4tb/hifi-human-wgs-wdl-custom/output/sv_integration_annotation
mkdir -p "$OUT"

declare -A SV_VCFS
SV_VCFS["TG-6102"]="/mnt/Ext_1tb_B/tg-integration-analysis/TG6102/20260227_103422_humanwgs_singleton/out/phased_sv_vcf/TG-6102.GRCm39.structural_variants.phased.vcf.gz"
SV_VCFS["TG-6283"]="/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6283_B01/20260519_092741_humanwgs_singleton/out/phased_sv_vcf/TG-6283.GRCm39.structural_variants.phased.vcf.gz"
SV_VCFS["TG-6903"]="/data_4tb/hifi-human-wgs-wdl-custom/output/TG6903_D01/20260529_113141_humanwgs_singleton/out/phased_sv_vcf/TG-6903.GRCm39.structural_variants.phased.vcf.gz"

declare -A SITES
SITES["TG-6102"]="chr6:90534206"
SITES["TG-6283"]="chr2:5559779"
SITES["TG-6903"]="chr4:68426842"

WINDOW=500000

echo "========================================================================"
echo "TG Integration Site — Complex SV Functional Annotation"
echo "Window: ±${WINDOW} bp | Reference: GRCm39 | Annotation: Gencode vM36"
echo "Generated: $(date '+%Y-%m-%d')"
echo "========================================================================"
echo ""

for SAMPLE in TG-6102 TG-6283 TG-6903; do
    VCF=${SV_VCFS[$SAMPLE]}
    IFS=: read CHR POS <<< "${SITES[$SAMPLE]}"
    START=$((POS - WINDOW))
    END=$((POS + WINDOW))

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "[ $SAMPLE ]  Integration site: $CHR:$POS  (window: $CHR:$START-$END)"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    # 1. PASS SV 추출
    SVBED=$OUT/${SAMPLE}.integration_svs.bed
    zcat "$VCF" 2>/dev/null | awk -v c="$CHR" -v s="$START" -v e="$END" '
      $0~/^#/ {next}
      $1==c && $2>=s && $2<=e && $7~/PASS/ {
        svtype=""; svlen="0"; endpos=$2
        n=split($8,a,";")
        for(i=1;i<=n;i++){
          if(a[i]~/^SVTYPE=/) svtype=substr(a[i],8)
          if(a[i]~/^SVLEN=/)  { svlen=substr(a[i],7); if(svlen<0) svlen=-svlen }
          if(a[i]~/^END=/)    endpos=substr(a[i],5)
        }
        bedend = (endpos>$2) ? endpos : $2+svlen
        if(bedend<=$2) bedend=$2+1
        printf "%s\t%d\t%d\t%s\t%s\t%s\n", $1,$2,bedend,svtype,svlen,$NF
      }' | sort -k1,1 -k2,2n > "$SVBED"

    SV_COUNT=$(wc -l < "$SVBED")
    echo ""
    printf "  %-10s %d 개\n" "PASS SV:" "$SV_COUNT"
    echo ""

    if [ "$SV_COUNT" -eq 0 ]; then
        echo "  (분석할 SV 없음)"
        echo ""
        continue
    fi

    # 2. SV 유형별 요약
    echo "  [ SV 유형별 분포 ]"
    awk '{print $4}' "$SVBED" | sort | uniq -c | sort -rn | \
      awk '{printf "    %-10s %d\n", $2, $1}'
    echo ""

    # 3. 유전자 overlap (Gencode GTF)
    echo "  [ 유전자 overlap (Gencode vM36) ]"
    printf "  %-6s %-10s %-10s %-25s %-20s %-15s %s\n" \
      "SV_POS" "SVTYPE" "SVLEN" "GENE" "GENE_TYPE" "FEATURE" "DISTANCE"

    while IFS=$'\t' read -r schr sstart send svtype svlen gt; do
        # GTF에서 해당 구간 overlap 유전자 조회
        tabix "$GTF" ${schr}:${sstart}-${send} 2>/dev/null | \
          awk -v ss="$sstart" -v se="$send" -v svt="$svtype" -v svl="$svlen" -v sp="$sstart" '$3~/gene|exon|UTR|CDS/ {
            n=""; t=""
            for(i=1;i<=NF;i++){
              if($i=="gene_name"){n=$(i+1); gsub(/[";]/,"",n)}
              if($i=="gene_type"){t=$(i+1); gsub(/[";]/,"",t)}
            }
            # overlap 여부
            if($5>=ss && $4<=se) {
              dist=0
            } else {
              dist=($4>se) ? $4-se : ss-$5
            }
            printf "  %-6s %-10s %-10s %-25s %-20s %-15s %d bp\n", sp, svt, svl, n, t, $3, dist
          }' | sort -k7 -n | awk '!seen[$4]++' | head -5
    done < "$SVBED"
    echo ""

    # 4. CpG island overlap
    echo "  [ CpG island overlap ]"
    CpG_HITS=0
    while IFS=$'\t' read -r schr sstart send svtype svlen gt; do
        hits=$(awk -v c="$schr" -v s="$sstart" -v e="$send" \
          'NR>1 && $1==c && $3>=s && $2<=e {print $1,$2,$3,$5}' "$CPG" | wc -l)
        if [ "$hits" -gt 0 ]; then
            CpG_HITS=$((CpG_HITS + hits))
            awk -v c="$schr" -v s="$sstart" -v e="$send" -v svt="$svtype" \
              'NR>1 && $1==c && $3>=s && $2<=e {
                printf "    %s:%s-%s (obs/exp=%.2f)\n", $1,$2,$3,$5
              }' "$CPG" | head -3
        fi
    done < "$SVBED"
    [ "$CpG_HITS" -eq 0 ] && echo "    (없음)"
    echo ""

    # 5. Simple repeat overlap
    echo "  [ Simple repeat overlap ]"
    REP_HITS=0
    while IFS=$'\t' read -r schr sstart send svtype svlen gt; do
        hits=$(zcat "$REPEAT" 2>/dev/null | awk -v c="$schr" -v s="$sstart" -v e="$send" \
          'NR>1 && $2==c && $4>=s && $3<=e {print $5}' | wc -l)
        REP_HITS=$((REP_HITS + hits))
    done < "$SVBED"
    printf "    %d 개 SV가 simple repeat 구간과 중첩\n" "$REP_HITS"
    echo ""

    # 6. Integration site 바로 인접 SV (±10kb) 강조
    echo "  [ Integration site 직접 인접 SV (±10kb) ]"
    awk -v p="$POS" '(p-$2)^2<=100000000 || ($3-p)^2<=100000000 {
      dist=($2<=p && $3>=p) ? 0 : (p<$2 ? $2-p : p-$3)
      if(dist<=10000) printf "    %s:%s-%s  %s  len=%s bp  dist=%d bp\n",$1,$2,$3,$4,$5,dist
    }' "$SVBED"
    echo ""

done

echo "========================================================================"
echo "출력 파일: $OUT/"
echo "========================================================================"
