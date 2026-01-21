#!/bin/bash
# create_batch_inputs.sh
# /data_4tb의 BAM 파일들로부터 자동으로 입력 JSON 파일 생성

RAWDATA_DIR="/data_4tb/pacbio_rawdata"
BATCH_INPUT_DIR="batch_inputs"
REF_MAP="/home/ygkim/ngs-pipeline/hifi-human-wgs-wdl-custom/hifi-wdl-resources/hifi-wdl-resources-v3.1.0/GRCh38.ref_map.v3p1p0.template.tsv"
TERTIARY_MAP="/home/ygkim/ngs-pipeline/hifi-human-wgs-wdl-custom/hifi-wdl-resources/hifi-wdl-resources-v3.1.0/GRCh38.tertiary_map.v3p1p0.template.tsv"

mkdir -p ${BATCH_INPUT_DIR}

# 샘플 정보 CSV 파일 (수동 작성 필요)
# 형식: sample_id,sex,bam_files
# 예: KTY9537,MALE,/data_4tb/pacbio_rawdata/.../file1.bam:/data_4tb/.../file2.bam

if [[ ! -f "samples.csv" ]]; then
    echo "Error: samples.csv not found"
    echo "Create samples.csv with format:"
    echo "sample_id,sex,bam_files"
    echo "KTY9537,MALE,/path/to/file1.bam:/path/to/file2.bam"
    exit 1
fi

# CSV 파일 읽기 (헤더 제외)
tail -n +2 samples.csv | while IFS=',' read -r sample_id sex bam_files; do
    echo "Creating input file for ${sample_id}..."
    
    # BAM 파일들을 배열로 변환
    IFS=':' read -ra BAM_ARRAY <<< "$bam_files"
    
    # JSON 배열 생성
    bam_json=""
    for bam in "${BAM_ARRAY[@]}"; do
        if [[ -z "$bam_json" ]]; then
            bam_json="\"${bam}\""
        else
            bam_json="${bam_json},\n    \"${bam}\""
        fi
    done
    
    # JSON 파일 생성
    cat > ${BATCH_INPUT_DIR}/${sample_id}.inputs.json << JSONEOF
{
  "humanwgs_singleton.sample_id": "${sample_id}",
  "humanwgs_singleton.sex": "${sex}",
  "humanwgs_singleton.hifi_reads": [
    ${bam_json}
  ],
  "humanwgs_singleton.ref_map_file": "${REF_MAP}",
  "humanwgs_singleton.tertiary_map_file": "${TERTIARY_MAP}",
  "humanwgs_singleton.backend": "HPC",
  "humanwgs_singleton.preemptible": false,
  "humanwgs_singleton.gpu": true,
  "humanwgs_singleton.max_reads_per_alignment_chunk": 100000000
}
JSONEOF
    
    echo "  Created: ${BATCH_INPUT_DIR}/${sample_id}.inputs.json"
done

echo ""
echo "✓ All input files created in ${BATCH_INPUT_DIR}/"
ls -lh ${BATCH_INPUT_DIR}/
