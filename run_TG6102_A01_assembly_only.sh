#!/bin/bash
# TG-6102 A01 de novo assembly only (hifiasm)
# reference-based mapping은 /mnt/Ext_1tb_B/tg-integration-analysis/TG6102/ 결과 재사용
#
# 실행 방법:
#   nohup bash run_TG6102_A01_assembly_only.sh \
#     > /mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6102_A01.assembly.log 2>&1 &

set -euo pipefail

SAMPLE="TG-6102"
RAW_BAM="/mnt/Ext_1tb_B/r84285_20260219_052427/1_A01/hifi_reads/m84285_260219_053344_s1.hifi_reads.bc2024.bam"
OUTDIR="/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6102_A01/assembly"
THREADS=32

SAMTOOLS_SIF="/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_singularity_cache/docker___quay.io_biocontainers_samtools_1.21--h50ea8bc_0.sif"
HIFIASM_SIF="/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_singularity_cache/docker___quay.io_biocontainers_hifiasm_0.25.0--h5ca1c30_0.sif"

APPTAINER="apptainer exec --bind /mnt,/data_4tb"

mkdir -p "$OUTDIR"
LOG="$OUTDIR/assembly.log"
exec > >(tee -a "$LOG") 2>&1

echo "========================================"
echo "TG-6102 A01 de novo assembly (hifiasm)"
echo "Sample:  $SAMPLE"
echo "시작: $(date)"
echo "========================================"

# ── Step 1: BAM → FASTQ ──────────────────────────────────────────────────────
FASTQ="$OUTDIR/${SAMPLE}.hifi_reads.fastq"
if [ ! -f "$FASTQ" ]; then
    echo ""
    echo "[Step 1] BAM → FASTQ 변환"
    $APPTAINER "$SAMTOOLS_SIF" \
        samtools fastq -@ 8 "$RAW_BAM" > "$FASTQ"
    echo "  완료: $(wc -l < "$FASTQ" | awk '{print $1/4}') reads"
else
    echo "[Step 1] FASTQ 이미 존재함, 건너뜀"
fi

# ── Step 2: hifiasm assembly ─────────────────────────────────────────────────
HAP1_GFA="$OUTDIR/${SAMPLE}.hifiasm.bp.hap1.p_ctg.gfa"
if [ ! -f "$HAP1_GFA" ]; then
    echo ""
    echo "[Step 2] hifiasm assembly"
    cd "$OUTDIR"
    $APPTAINER "$HIFIASM_SIF" \
        hifiasm \
        -o "${SAMPLE}.hifiasm" \
        -t "$THREADS" \
        "$FASTQ"
    echo "  완료: $(date)"
else
    echo "[Step 2] GFA 이미 존재함, 건너뜀"
fi

# ── Step 3: GFA → FASTA ──────────────────────────────────────────────────────
echo ""
echo "[Step 3] GFA → FASTA 변환"
for HAP in hap1 hap2; do
    GFA="$OUTDIR/${SAMPLE}.hifiasm.bp.${HAP}.p_ctg.gfa"
    FA="$OUTDIR/${SAMPLE}.hifiasm.bp.${HAP}.p_ctg.fa"
    if [ ! -f "$FA" ]; then
        awk '/^S/{print ">"$2; print $3}' "$GFA" > "$FA"
        CTGS=$(grep -c "^>" "$FA")
        SIZE=$(awk '/^>/{next}{total+=length($0)}END{printf "%.1f Mb", total/1e6}' "$FA")
        echo "  ${HAP}: ${CTGS} contigs, ${SIZE}"
    else
        echo "  ${HAP}: 이미 존재함"
    fi
done

echo ""
echo "========================================"
echo "완료: $(date)"
echo "결과:"
ls -lh "$OUTDIR/${SAMPLE}.hifiasm.bp.hap"*.p_ctg.fa
echo ""
echo "다음 단계: transgene integration 분석"
echo "  cd /data_4tb/hifi-human-wgs-wdl-custom"
echo "  bash run_tg_integration.sh TG-6102 TG6102.tg_integration.inputs.json \\"
echo "    /mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6102_A01_integration_analysis"
echo "========================================"
