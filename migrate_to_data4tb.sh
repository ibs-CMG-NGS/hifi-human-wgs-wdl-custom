#!/bin/bash
# 파이프라인 작업 디렉토리를 /data_4tb로 이동하고 심볼릭 링크 생성

set -e  # 에러 발생 시 중단

CURRENT_DIR=~/ngs-pipeline/hifi-human-wgs-wdl-custom
TARGET_DIR=/data_4tb/hifi-human-wgs-wdl-custom

echo "=========================================="
echo "HiFi-WGS 파이프라인 작업공간 이동"
echo "=========================================="
echo "원본: $CURRENT_DIR"
echo "대상: $TARGET_DIR"
echo ""

# 1. 대상 디렉토리 생성
echo "[1/6] 대상 디렉토리 생성..."
mkdir -p $TARGET_DIR
echo "✓ 완료"

# 2. 실행 디렉토리들 이동 (가장 큰 용량)
echo ""
echo "[2/6] 실행 디렉토리 이동 중..."
if ls $CURRENT_DIR/2026*_humanwgs_singleton &> /dev/null; then
    echo "  - 모든 실행 디렉토리 이동 중..."
    for dir in $CURRENT_DIR/2026*_humanwgs_singleton; do
        if [ -d "$dir" ]; then
            echo "    이동: $(basename $dir) ($(du -sh $dir | cut -f1))"
            mv "$dir" $TARGET_DIR/
        fi
    done
    echo "✓ 완료"
else
    echo "  (실행 디렉토리 없음 - 건너뜀)"
fi

# 3. Call cache 이동 (완료된 작업 보존 - 매우 중요!)
echo ""
echo "[3/6] Call cache 이동 중..."
if [ -L "$CURRENT_DIR/miniwdl_call_cache" ]; then
    echo "  (이미 심볼릭 링크 존재 - 건너뜀)"
elif [ -d "$CURRENT_DIR/miniwdl_call_cache" ]; then
    mv $CURRENT_DIR/miniwdl_call_cache $TARGET_DIR/
    ln -sf $TARGET_DIR/miniwdl_call_cache $CURRENT_DIR/miniwdl_call_cache
    echo "✓ 완료 (심볼릭 링크 생성)"
else
    mkdir -p $TARGET_DIR/miniwdl_call_cache
    ln -sf $TARGET_DIR/miniwdl_call_cache $CURRENT_DIR/miniwdl_call_cache
    echo "✓ 새로 생성 (심볼릭 링크 생성)"
fi

# 4. Singularity 캐시 이동 (17GB)
echo ""
echo "[4/6] Singularity 캐시 이동 중..."
if [ -L "$CURRENT_DIR/miniwdl_singularity_cache" ]; then
    echo "  (이미 심볼릭 링크 존재 - 건너뜀)"
elif [ -d "$CURRENT_DIR/miniwdl_singularity_cache" ]; then
    mv $CURRENT_DIR/miniwdl_singularity_cache $TARGET_DIR/
    ln -sf $TARGET_DIR/miniwdl_singularity_cache $CURRENT_DIR/miniwdl_singularity_cache
    echo "✓ 완료 (심볼릭 링크 생성)"
else
    mkdir -p $TARGET_DIR/miniwdl_singularity_cache
    ln -sf $TARGET_DIR/miniwdl_singularity_cache $CURRENT_DIR/miniwdl_singularity_cache
    echo "✓ 새로 생성 (심볼릭 링크 생성)"
fi

# 5. 심볼릭 링크 확인
echo ""
echo "[5/6] 심볼릭 링크 확인..."
if [ -L "$CURRENT_DIR/miniwdl_call_cache" ] && [ -L "$CURRENT_DIR/miniwdl_singularity_cache" ]; then
    echo "✓ 심볼릭 링크 정상"
else
    echo "⚠ 심볼릭 링크 확인 필요"
fi

# 6. 검증
echo ""
echo "[6/6] 이동 완료 확인..."
echo ""
df -h $TARGET_DIR
echo ""
du -sh $TARGET_DIR/*
echo ""
echo "=========================================="
echo "✅ 이동 완료!"
echo "=========================================="
echo ""
echo "다음 단계:"
echo "1. 설정 파일 업데이트: ./update_config_paths.sh"
echo "2. 파이프라인 재실행"
echo ""
