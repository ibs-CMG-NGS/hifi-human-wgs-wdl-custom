#!/usr/bin/env bash
# base_modification_per_sample.wdl을 샘플 하나씩 순차 실행 (5mC/5hmC per-sample).
# - run 2/3 (8개): raw BAM → jasmine → MM 필터 + 정렬 → CpG pileup
# - run 1 (4개): 파일럿의 jasmine+정렬 BAM(mm_tagged)으로 pileup만
# 샘플이 성공하면 미정렬 jasmine 출력(원본×0.18, 샘플당 ~50-70G)을 삭제해
# OIM 최대 사용량을 ~0.6T로 유지한다. 정렬 BAM은 남긴다.
# 사용: setsid nohup bash run_5hmC_per_sample.sh > <OUTROOT>/run_5hmC_per_sample.log 2>&1 < /dev/null &
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "${HERE}/.." && pwd)"
SAMPLES="${HERE}/samples.tsv"
OUTROOT="${OUTROOT:-/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/5hmC}"
WDL="${REPO}/workflows/base_modification_per_sample.wdl"
REF="${REPO}/hifi-wdl-resources/GRCm39/mouse_GRCm39.fasta"
PILOT=/mnt/OIM_HDD2/base_mod_discovery_pilot_output

# run 1 파일럿 BAM (그룹당 BAM 1개였으므로 샘플 단위와 동일)
declare -A EXISTING=(
  [H2O2-A01-ctrl1]="${PILOT}/20260706_fix/20260706_153742_base_modification_discovery/call-filter_ctrl/work/control.mm_tagged.bam"
  [H2O2-B01]="${PILOT}/20260706_fix/20260706_153742_base_modification_discovery/call-filter_drug_a/work/drug_a.mm_tagged.bam"
  [H2O2-C01]="${PILOT}/20260706_known/20260706_214242_base_modification_known/call-filter_drug_b/work/drug_b.mm_tagged.bam"
  [H2O2-D01]="${PILOT}/20260706_fix/20260706_153742_base_modification_discovery/call-filter_drug_ab/work/drug_ab.mm_tagged.bam"
)

source /home/ygkim/program/anaconda3/etc/profile.d/conda.sh
conda activate hifi-human-wgs
cd "${REPO}"
mkdir -p "${OUTROOT}"

# 순서: run 2, run 3 (jasmine 필요, 장시간) → run 1 (pileup만)
ORDER=$(tail -n +2 "${SAMPLES}" | awk -F'\t' '{k=($2==1)?9:$2; print k"\t"$0}' | sort -s -k1,1n | cut -f2-)

while IFS=$'\t' read -r sid run group bam raw; do
  sdir="${OUTROOT}/${sid}"
  if [ -f "${sdir}/DONE" ]; then
    echo "[$(date '+%F %T')] skip ${sid} (done)"
    continue
  fi
  mkdir -p "${sdir}"

  inputs="${sdir}/inputs.json"
  if [ -n "${EXISTING[${sid}]:-}" ]; then
    src="${EXISTING[${sid}]}"
    bam_json="\"base_modification_per_sample.existing_modbam\": \"${src}\", \"base_modification_per_sample.existing_modbai\": \"${src}.bai\","
  else
    bam_json="\"base_modification_per_sample.raw_hifi_bam\": \"${raw}\","
  fi
  cat > "${inputs}" <<EOF
{
  "base_modification_per_sample.sample_id": "${sid}",
  ${bam_json}
  "base_modification_per_sample.ref_fasta": "${REF}",
  "base_modification_per_sample.ref_fai": "${REF}.fai",
  "base_modification_per_sample.min_coverage": 5,
  "base_modification_per_sample.filter_threshold": 0.7,
  "base_modification_per_sample.default_runtime_attributes": {
    "backend": "HPC",
    "container_namespace": null,
    "container_registry": "quay.io/pacbio",
    "cpuPlatform": "",
    "gpuType": "",
    "max_retries": 0,
    "preemptible_tries": 0,
    "zones": ""
  }
}
EOF

  echo "[$(date '+%F %T')] start ${sid} (run ${run}, group ${group}) $( [ -n "${EXISTING[${sid}]:-}" ] && echo 'pileup-only' || echo 'jasmine+align+pileup')"
  miniwdl run "${WDL}" --input "${inputs}" --dir "${sdir}" > "${sdir}/miniwdl.log" 2>&1
  rc=$?
  if [ ${rc} -ne 0 ]; then
    echo "[$(date '+%F %T')] FAIL ${sid} (exit ${rc}) — see ${sdir}/miniwdl.log; stopping"
    exit ${rc}
  fi

  # 정렬 성공 후 미정렬 jasmine 출력 삭제 (out/ 링크 포함)
  find "${sdir}" -path '*/call-jasmine_unaligned/work/*.jasmine.bam' -type f -print -delete
  find "${sdir}" -path '*/out/jasmine_bam/*' \( -type l -o -type f \) -delete 2>/dev/null || true

  touch "${sdir}/DONE"
  echo "[$(date '+%F %T')] done ${sid}; OIM free: $(df -h /mnt/OIM_HDD2 | awk 'NR==2{print $4}')"
done <<< "${ORDER}"

echo "[$(date '+%F %T')] all samples finished"
