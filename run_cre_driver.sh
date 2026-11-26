#!/bin/bash
# 2026 Cre driver 프로젝트 — pv-cre / camk-cre 의 2차 분석 (humanwgs_singleton).
#
# 준비됐지만 **아직 돌리지 않았다** (2026-11-26). 돌리기 전에 확인할 것:
#   - /data_4tb 여유 용량. 같은 run 의 1_C01 / 1_D01 이 전송되면 디스크를
#     더 채운다. **C01/D01 은 이 프로젝트가 아니다** — 또 다른 TG 라인이고
#     Cre 드라이버 TG 도 아니다(2026-11-26 확인). 이 스크립트는 1_A01·1_B01
#     만 쓰므로 전송 중이어도 되고, 전송을 기다릴 이유도 없다
#   - run_assembly 를 true 로 뒀다 (아래 주석)
#
# 출력은 **/data24** 로 간다. 샘플당 final 산출물이 67–95 GB 이고 miniwdl
# 작업 디렉터리는 그보다 크다. /data_3tb 는 409 G, /data_4tb 는 시퀀서가
# 쓰는 중이고, /data24 는 14 T 남아 있다.
#
# 컨테이너 캐시도 /data24 다. /data_4tb 에 두었다가 2026-11-12 에 99% 가 차서
# 캐시가 지워졌다 — 그걸 반복하지 않는다.
set -o pipefail
source /home/ygkim/program/anaconda3/etc/profile.d/conda.sh
conda activate hifi-human-wgs

WORK_DIR="/data_4tb/hifi-human-wgs-wdl-custom"
OUT_ROOT="/data24/processed/wgs/TG/secondary/2026-kyg-cre-driver"

# miniwdl 이 미리 받아 둔 이미지를 쓰게 한다. 이게 없으면 첫 실행이 32개
# 이미지를 전부 다시 내려받는다.
export MINIWDL__SINGULARITY__IMAGE_CACHE=/data24/software/miniwdl_singularity_cache

cd "${WORK_DIR}" || exit 1

# main 에 transgene_integration 수정(1bb167f)이 **들어와 있어야 한다.**
# 2026-11-26 에 fast-forward 머지했다. 확인:
git -C "${WORK_DIR}" merge-base --is-ancestor 1bb167f HEAD \
  || { echo "ERROR: HEAD 에 1bb167f(TG 접합 수정)가 없다. 이대로 돌리면 strand 버그를 다시 맞는다"; exit 1; }

for S in pv-cre camk-cre; do
    OUTPUT_DIR="${OUT_ROOT}/${S}_wgs"
    LOG_FILE="${OUTPUT_DIR}.nohup.log"
    mkdir -p "${OUTPUT_DIR}" || exit 1
    echo "=== ${S} humanwgs_singleton start $(date) ==="
    miniwdl run workflows/singleton.wdl \
        --input "${S}.inputs.json" \
        --dir "${OUTPUT_DIR}" \
        --verbose 2>&1 | tee "${LOG_FILE}"
    echo "=== ${S} exit ${PIPESTATUS[0]} at $(date) ==="
done

# 끝난 뒤 3차로 넘기는 경로 (outputs.json 을 프로젝트 config 에 적는다):
#   config/projects/2026-kyg-cre-driver.yaml 의 secondary_outputs_json
#   -> ${OUT_ROOT}/<sample>_wgs/<timestamp>_humanwgs_singleton/outputs.json
# 그다음:
#   PYTHONPATH=src python -m core.cli ingest    --project config/projects/2026-kyg-cre-driver.yaml
#   PYTHONPATH=src python -m core.cli preflight --project config/projects/2026-kyg-cre-driver.yaml
