#!/bin/bash
# =============================================================================
# DeepVariant 메틸화 인식(methylation-aware) 재실행 — 단일 웰 테스트 (A01)
# -----------------------------------------------------------------------------
# 목적: deepvariant.wdl make_examples에 추가한 --enable_methylation_calling 플래그가
#       small variant VCF의 MF/MD/MT FORMAT 필드를 실제로 채우는지 검증.
#
# 방식: 전체 파이프라인(수 시간) 대신 DeepVariant 서브워크플로우만 단독 실행.
#       기존 실행에서 이미 만들어진 merged aligned BAM(정렬 완료본)을 입력으로 재사용
#       → 정렬/SV/phasing/메틸화 pileup을 전부 건너뛰어 GPU로 ~1시간 수준.
#
# 주의: 기존 원본 실행 결과는 절대 건드리지 않음(새 출력 디렉토리에만 기록).
# =============================================================================
set -euo pipefail

# ── 설정 ────────────────────────────────────────────────────────────────────
REPO="/data_4tb/hifi-human-wgs-wdl-custom"
WDL="${REPO}/workflows/wdl-common/wdl/workflows/deepvariant/deepvariant.wdl"
INPUTS="${REPO}/rerun_methylation/deepvariant_A01_methyl.inputs.json"
OUTDIR="/mnt/hdd1_1tb/h2o2-wgs/methyl_rerun/A01"   # 여유 디스크(432G) 위치
CONDA_ENV="hifi-human-wgs"

echo "=========================================="
echo " DeepVariant methylation-aware 재실행 (A01 테스트)"
echo "  WDL     : ${WDL}"
echo "  inputs  : ${INPUTS}"
echo "  출력    : ${OUTDIR}"
echo "=========================================="

# ── conda 환경 ──────────────────────────────────────────────────────────────
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
command -v miniwdl >/dev/null || { echo "ERROR: miniwdl 없음 (conda env: ${CONDA_ENV})"; exit 1; }

# ── 0) 사전점검: 플래그가 실제로 들어갔는지 + WDL 문법 검사 ──────────────────
echo "[0/2] WDL 검증"
if ! grep -q -- "--enable_methylation_calling" "${WDL}"; then
  echo "ERROR: ${WDL} 에 --enable_methylation_calling 플래그가 없습니다. 먼저 WDL 수정 필요."; exit 1
fi
miniwdl check "${WDL}"    # 문법/타입 오류 있으면 여기서 중단

# ── 1) 실행 ─────────────────────────────────────────────────────────────────
echo "[1/2] miniwdl run 시작 (GPU DeepVariant)"
mkdir -p "${OUTDIR}"
miniwdl run "${WDL}" \
  --input "${INPUTS}" \
  --dir "${OUTDIR}" \
  --verbose 2>&1 | tee "${OUTDIR}/run.log"

# ── 2) 검증: MF/MD/MT가 이제 채워졌는지 확인 ────────────────────────────────
echo "[2/2] 결과 검증 — MF/MD/MT 필드 채워짐 여부"
OUTVCF=$(find "${OUTDIR}" -iname "*.small_variants.vcf.gz" ! -iname "*.g.vcf.gz" 2>/dev/null | head -1)
[ -z "${OUTVCF}" ] && OUTVCF=$(find "${OUTDIR}" -iname "*.small_variants*.vcf.gz" 2>/dev/null | head -1)
echo "  출력 VCF: ${OUTVCF}"

echo "  --- FORMAT 필드에 MF/MD/MT 포함 여부 (이형접합 레코드 3개) ---"
bcftools view -H "${OUTVCF}" 2>/dev/null | awk -F'\t' '$9 ~ /MF|MD|MT/ {print $9"  =>  "$10; c++} c>=3{exit}'

echo "  --- MF 값이 실제로 채워진 레코드 수 (예: MF에 . 아닌 값) ---"
bcftools query -f '[%MF]\n' "${OUTVCF}" 2>/dev/null | grep -vcE '^\.?(,\.)*$' || true

echo "=========================================="
echo " 완료. 위에서 FORMAT에 MF:MD:MT가 보이고 값이 채워졌으면 성공."
echo " 성공 시 → 나머지 11개 웰로 확장하거나, phasing(HiPhase)까지 이어서 재실행."
echo "=========================================="
