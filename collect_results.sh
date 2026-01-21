#!/bin/bash
# collect_results.sh
# 모든 샘플의 주요 결과를 한 곳에 모으기

BATCH_RESULTS="/data_4tb/hifi-human-wgs-wdl-custom/batch_results"
SUMMARY_DIR="${BATCH_RESULTS}/summary"

mkdir -p ${SUMMARY_DIR}/{vcfs,bams,reports,stats}

echo "Collecting results from ${BATCH_RESULTS}..."

for sample_dir in ${BATCH_RESULTS}/*/; do
    sample=$(basename "${sample_dir}")
    
    # 로그나 summary 디렉토리는 건너뜀
    if [[ "$sample" == "logs" ]] || [[ "$sample" == "summary" ]]; then
        continue
    fi
    
    out_dir="${sample_dir}/out"
    
    if [[ ! -d "$out_dir" ]]; then
        echo "⚠ No output directory for ${sample}"
        continue
    fi
    
    echo "Processing ${sample}..."
    
    # VCF 파일 복사
    find "${out_dir}" -name "*.vcf.gz" -exec cp {} ${SUMMARY_DIR}/vcfs/ \; 2>/dev/null
    
    # BAM 파일 심볼릭 링크 (용량 절약)
    find "${out_dir}" -name "*.bam" -exec ln -sf {} ${SUMMARY_DIR}/bams/ \; 2>/dev/null
    
    # HTML 리포트 복사
    find "${out_dir}" -name "*.html" -exec cp {} ${SUMMARY_DIR}/reports/ \; 2>/dev/null
    
    # 통계 파일 복사
    find "${out_dir}" -name "*stats*.txt" -o -name "*summary*.txt" | while read file; do
        cp "$file" "${SUMMARY_DIR}/stats/$(basename $(dirname $(dirname $file)))_$(basename $file)" 2>/dev/null
    done
done

echo ""
echo "✓ Results collected in ${SUMMARY_DIR}/"
echo ""
echo "Summary:"
echo "  VCFs: $(ls -1 ${SUMMARY_DIR}/vcfs/*.vcf.gz 2>/dev/null | wc -l) files"
echo "  BAMs: $(ls -1 ${SUMMARY_DIR}/bams/*.bam 2>/dev/null | wc -l) files"
echo "  Reports: $(ls -1 ${SUMMARY_DIR}/reports/*.html 2>/dev/null | wc -l) files"
echo "  Stats: $(ls -1 ${SUMMARY_DIR}/stats/*.txt 2>/dev/null | wc -l) files"
