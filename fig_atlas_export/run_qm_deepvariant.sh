#!/usr/bin/env bash
# 품질 맞춘 chr1–3 BAM에 DeepVariant 1.9.0(원 파이프라인과 같은 make_examples 인자)과
# 조건별 GLnexus(원 파이프라인과 같은 config)를 실행한다.
# make_examples는 CPU(nice, 샘플당 MAKE_JOBS개), call_variants는 GPU(샘플을 GPU 0/1에 번갈아 배정).
# 사용: run_qm_deepvariant.sh QM_DIR   (QM_DIR/bams/*.qm.chr1-3.bam → QM_DIR/dv/, QM_DIR/joint/)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "${HERE}/.." && pwd)"
QM="$1"
REF="${REPO}/hifi-wdl-resources/GRCm39/mouse_GRCm39.fasta"
DV_SIF="${REPO}/miniwdl_singularity_cache/docker___google_deepvariant_1.9.0-gpu.sif"
GLX_SIF=$(ls "${REPO}"/miniwdl_singularity_cache/docker___quay.io_pacbio_glnexus@*.sif | head -1)
MAKE_JOBS="${MAKE_JOBS:-6}"
BIND=(--bind /mnt,/data_4tb)
mkdir -p "$QM/dv" "$QM/joint"
awk -F'\t' -v OFS='\t' '$1=="chr1"||$1=="chr2"||$1=="chr3"{print $1, 0, $2}' "$REF.fai" > "$QM/regions.bed"

dv_one() {
  sid="$1"; gpu="$2"; d="$QM/dv/$sid"
  [ -f "$d/$sid.vcf.gz.tbi" ] && { echo "skip $sid"; return 0; }
  rm -rf "$d"; mkdir -p "$d/ex" "$d/gv"
  echo "[$(date '+%F %T')] make_examples $sid"
  seq 0 $((MAKE_JOBS - 1)) | nice -n 10 apptainer exec "${BIND[@]}" "$DV_SIF" parallel --jobs "$MAKE_JOBS" --halt 2 \
    /opt/deepvariant/bin/make_examples \
      --checkpoint /opt/models/pacbio --norealign_reads --trim_reads_for_pileup \
      --vsc_min_fraction_indels 0.12 --pileup_image_width 147 --track_ref_reads --phase_reads \
      --partition_size=25000 --max_reads_per_partition=600 --alt_aligned_pileup=diff_channels \
      --sort_by_haplotypes --parse_sam_aux_fields --enable_methylation_calling \
      --methylation_calling_threshold 0.5 --min_mapping_quality=1 --mode calling \
      --ref "$REF" --regions "$QM/regions.bed" --reads "$QM/bams/$sid.qm.chr1-3.bam" \
      --examples "$d/ex/make_examples.tfrecord@${MAKE_JOBS}.gz" \
      --gvcf "$d/gv/gvcf.tfrecord@${MAKE_JOBS}.gz" --task {} > "$d/make_examples.log" 2>&1
  echo "[$(date '+%F %T')] call_variants $sid (GPU $gpu)"
  CUDA_VISIBLE_DEVICES="$gpu" apptainer exec --nv "${BIND[@]}" "$DV_SIF" /opt/deepvariant/bin/call_variants \
    --outfile "$d/call_variants_output.tfrecord.gz" --examples "$d/ex/make_examples.tfrecord@${MAKE_JOBS}.gz" \
    --checkpoint /opt/models/pacbio > "$d/call_variants.log" 2>&1
  nice -n 10 apptainer exec "${BIND[@]}" "$DV_SIF" /opt/deepvariant/bin/postprocess_variants --cpus 2 \
    --vcf_stats_report=false --ref "$REF" --infile "$d/call_variants_output.tfrecord.gz" \
    --outfile "$d/$sid.vcf.gz" --nonvariant_site_tfrecord_path "$d/gv/gvcf.tfrecord@${MAKE_JOBS}.gz" \
    --gvcf_outfile "$d/$sid.g.vcf.gz" > "$d/postprocess.log" 2>&1
  rm -rf "$d/ex" "$d/gv" "$d"/call_variants_output*.tfrecord.gz
  echo "[$(date '+%F %T')] done $sid"
}

# 샘플 2개씩 동시에 (GPU 0, 1)
mapfile -t SIDS < <(tail -n +2 "${REPO}/base_mod_per_sample/samples.tsv" | cut -f1)
for ((i = 0; i < ${#SIDS[@]}; i += 2)); do
  dv_one "${SIDS[i]}" 0 &
  [ $((i + 1)) -lt ${#SIDS[@]} ] && dv_one "${SIDS[i + 1]}" 1 &
  wait
done

# 조건별 GLnexus (원 파이프라인 config: DeepVariant_unfiltered 기반, required_dp 1)
cat > "$QM/joint/config.yml" <<'YML'
unifier_config:
  min_AQ1: 0
  min_AQ2: 0
  min_GQ: 0
  monoallelic_sites_for_lost_alleles: true
  max_alleles_per_site: 32
genotyper_config:
  required_dp: 1
  revise_genotypes: false
  allow_partial_data: true
  more_PL: true
  trim_uncalled_alleles: true
  liftover_fields:
    - orig_names: [MIN_DP, DP]
      name: DP
      description: '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">'
      type: int
      combi_method: min
      number: basic
      count: 1
      ignore_non_variants: true
    - orig_names: [AD]
      name: AD
      description: '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">'
      type: int
      number: alleles
      combi_method: min
      default_type: zero
      count: 0
    - orig_names: [GQ]
      name: GQ
      description: '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'
      type: int
      number: basic
      combi_method: min
      count: 1
      ignore_non_variants: true
    - orig_names: [PL]
      name: PL
      description: '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype Likelihoods">'
      type: int
      number: genotype
      combi_method: missing
      count: 0
      ignore_non_variants: true
YML
# 이 BAM들의 gVCF에는 methylation 호출이 만든 NoCall 행(ALT ".,<*>", CpG 위치 수백만 개)이 있어
# GLnexus가 "allele is not a DNA sequence"로 거부한다. 원 파이프라인 gVCF에서는 같은 위치가
# 참조 블록이므로, 이 행을 0/0 참조 블록 형식(ALT <*>, GT:GQ:MIN_DP:PL)으로 바꾼 사본을 쓴다.
for sid in "${SIDS[@]}"; do
  d="$QM/dv/$sid"; f="$d/$sid.glx.g.vcf.gz"
  [ -f "$f.tbi" ] && continue
  zcat "$d/$sid.g.vcf.gz" | awk -F'\t' -v OFS='\t' '
    /^#/ {print; next}
    $5 ~ /^\.,/ {
      n = split($9, k, ":"); split($10, v, ":")
      for (i = 1; i <= n; i++) x[k[i]] = v[i]
      split(x["PL"], pl, ",")
      print $1, $2, $3, $4, "<*>", 0, ".", "END=" $2, "GT:GQ:MIN_DP:PL", "0/0:0:" x["DP"] ":" pl[1] "," pl[2] "," pl[3]
      delete x; next
    }
    {print}' > "$d/$sid.glx.g.vcf"
  apptainer exec "${BIND[@]}" "$GLX_SIF" bcftools view -Oz -o "$f" "$d/$sid.glx.g.vcf"
  apptainer exec "${BIND[@]}" "$GLX_SIF" bcftools index -t "$f"
  rm -f "$d/$sid.glx.g.vcf"
done

for g in A B C D; do
  o="$QM/joint/$g.vcf.gz"
  [ -f "$o.tbi" ] && continue
  gv=$(tail -n +2 "${REPO}/base_mod_per_sample/samples.tsv" | awk -F'\t' -v g=$g -v q="$QM" '$3==g{printf "%s/dv/%s/%s.glx.g.vcf.gz ", q, $1, $1}')
  rm -rf "$QM/joint/$g.DB"
  ( cd "$QM/joint" && nice -n 10 apptainer exec "${BIND[@]}" "$GLX_SIF" glnexus_cli --threads 4 --mem-gbytes 16 \
      --dir "$g.DB" --config ./config.yml --bed "$QM/regions.bed" $gv > "$g.bcf" 2> "$g.glnexus.log" )
  apptainer exec "${BIND[@]}" "$GLX_SIF" bcftools view -Oz -o "$o" "$QM/joint/$g.bcf"
  apptainer exec "${BIND[@]}" "$GLX_SIF" bcftools index -t "$o"
  rm -rf "$QM/joint/$g.DB" "$QM/joint/$g.bcf"
  echo "[$(date '+%F %T')] joint $g done"
done
echo "[$(date '+%F %T')] all done"
