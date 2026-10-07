#!/usr/bin/env bash
# 샘플 BAM의 한 구간에서 primary read별 rq(read 품질), np(pass 수), 길이를 뽑는다.
# 사용: read_quality_profile.sh SAMPLE_ID BAM REGION OUT_TSV
S=/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_singularity_cache/docker___quay.io_biocontainers_samtools_1.21--h50ea8bc_0.sif
apptainer exec --bind /mnt,/data_4tb "$S" samtools view -F 0x900 "$2" "$3" | awk -v s="$1" '
  { rq = ""; np = ""
    for (i = 12; i <= NF; i++) { if ($i ~ /^rq:f:/) rq = substr($i, 6); if ($i ~ /^np:i:/) np = substr($i, 6) }
    print s "\t" rq "\t" np "\t" length($10) }' > "$4"
