#!/bin/bash
# miniwdl 설정 파일의 경로를 /data_4tb로 업데이트

set -e

CONFIG_FILE="config/miniwdl.local.cfg"
MINIWDL_CFG=".miniwdl.cfg"
TARGET_DIR="/data_4tb/hifi-human-wgs-wdl-custom"

echo "=========================================="
echo "miniwdl 설정 파일 경로 업데이트"
echo "=========================================="

# config/miniwdl.local.cfg 업데이트
echo "[1/2] $CONFIG_FILE 업데이트..."
sed -i 's|dir = "$PWD/miniwdl_call_cache"|dir = "'"$TARGET_DIR"'/miniwdl_call_cache"|g' $CONFIG_FILE
sed -i 's|image_cache = "$PWD/miniwdl_singularity_cache"|image_cache = "'"$TARGET_DIR"'/miniwdl_singularity_cache"|g' $CONFIG_FILE
echo "✓ 완료"

# .miniwdl.cfg 업데이트
echo "[2/2] $MINIWDL_CFG 업데이트..."
sed -i 's|dir = "$PWD/miniwdl_call_cache"|dir = "'"$TARGET_DIR"'/miniwdl_call_cache"|g' $MINIWDL_CFG
sed -i 's|image_cache = "$PWD/miniwdl_singularity_cache"|image_cache = "'"$TARGET_DIR"'/miniwdl_singularity_cache"|g' $MINIWDL_CFG
echo "✓ 완료"

echo ""
echo "=========================================="
echo "✅ 설정 파일 업데이트 완료!"
echo "=========================================="
echo ""
echo "확인:"
grep -E "dir =|image_cache =" $CONFIG_FILE
echo ""
echo "이제 파이프라인을 실행할 수 있습니다:"
echo "  cd ~/ngs-pipeline/hifi-human-wgs-wdl-custom"
echo "  export CUDA_VISIBLE_DEVICES=1"
echo "  miniwdl run workflows/singleton.wdl --input sample1.inputs.json --cfg config/miniwdl.local.cfg --verbose"
echo ""
