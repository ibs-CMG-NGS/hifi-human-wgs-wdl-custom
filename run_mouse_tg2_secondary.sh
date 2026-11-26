#!/bin/bash
# 2026-11 마우스 TG 통합 프로젝트 — 2차 분석 (humanwgs_singleton).
# 3차 프로젝트: config/projects/2026-kyg-mouse-tg-2.yaml (wgs-tertiary-transgene-analysis)
#
# 준비됐지만 **아직 돌리지 않았다** (2026-11-26).
#
# **네 샘플을 다 모은 뒤 한 번에 돌린다** (2026-11-26 사용자 결정). 나눠
# 돌려도 3차 결과는 같지만 — case 는 서로 독립이고 샘플 간 유일한 결합은
# 대조군에서 나온 control_junctions 하나다 — sample_set_hash 가 달라져 같은
# project_id 에 번들이 둘 생기고, 그러면 어느 쪽이 정본인지를 따로 말해야
# 한다. 한 번에 돌려서 그 문제를 없앤다.
#
# 돌리기 전에:
#   - 1_C01 / 1_D01 전송 완료. 2026-11-26 14:27 기준 C01 9.5 GB 진행,
#     D01 은 시작 전. /data_4tb 여유 용량을 같이 본다
#   - 그 두 웰의 입력 json 생성:
#       python3 scripts/well_inputs.py /data_4tb/pacbio_rawdata/r84285_20261124_025446 \
#           1_C01 1_D01 --write
#     **샘플 이름을 손으로 적지 말 것** — 웰 metadata 의 WellSample Name 이
#     이름이고, 손으로 적다가 틀린 적이 두 번 있다 (P8). 생성기는 A01/B01 의
#     기존 json 을 바이트 동일하게 재현하는 것으로 검증했다
#   - 생성한 json 을 `git add -f` 로 추적에 넣는다 (.gitignore 가
#     *.inputs.json 을 무시한다)
#   - 아래 SAMPLES 에 두 이름을 추가한다
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
OUT_ROOT="/data24/processed/wgs/TG/secondary/2026-kyg-mouse-tg-2"

# 인자로 샘플을 주면 그것만 돈다. 안 주면 아래 목록을 쓴다.
# 네 샘플 모두 well_inputs.py 로 metadata 에서 읽은 이름이다 (2026-11-26).
SAMPLES=("$@")
if [ ${#SAMPLES[@]} -eq 0 ]; then
    SAMPLES=(pv-cre camk-cre aldh1l1-egfp hba-floxed)
fi

# miniwdl 이 미리 받아 둔 이미지를 쓰게 한다. 이게 없으면 첫 실행이 33개
# 이미지를 전부 다시 내려받는다.
export MINIWDL__SINGULARITY__IMAGE_CACHE=/data24/software/miniwdl_singularity_cache

cd "${WORK_DIR}" || exit 1

# main 에 transgene_integration 수정(1bb167f)이 **들어와 있어야 한다.**
# 2026-11-26 에 fast-forward 머지했다. 확인:
git -C "${WORK_DIR}" merge-base --is-ancestor 1bb167f HEAD \
  || { echo "ERROR: HEAD 에 1bb167f(TG 접합 수정)가 없다. 이대로 돌리면 strand 버그를 다시 맞는다"; exit 1; }

# 입력이 하나라도 없으면 **시작 전에** 멈춘다. 몇 시간 돌린 뒤 마지막
# 샘플에서 죽는 것이 제일 나쁘다.
missing=0
for S in "${SAMPLES[@]}"; do
    [ -r "${S}.inputs.json" ] || { echo "ERROR: ${S}.inputs.json 이 없다"; missing=1; }
done
[ "${missing}" -eq 0 ] || { echo "well_inputs.py 로 만들고 다시 돌린다"; exit 1; }

echo "=== 샘플 ${#SAMPLES[@]}개: ${SAMPLES[*]}"
for S in "${SAMPLES[@]}"; do
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
#   config/projects/2026-kyg-mouse-tg-2.yaml 의 secondary_outputs_json
#   -> ${OUT_ROOT}/<sample>_wgs/<timestamp>_humanwgs_singleton/outputs.json
# 그다음:
#   PYTHONPATH=src python -m core.cli ingest    --project config/projects/2026-kyg-mouse-tg-2.yaml
#   PYTHONPATH=src python -m core.cli preflight --project config/projects/2026-kyg-mouse-tg-2.yaml
