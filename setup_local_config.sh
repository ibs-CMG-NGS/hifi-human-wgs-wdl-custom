#!/bin/bash
# setup_local_config.sh - 40코어 서버에 맞는 설정 파일 생성

echo "🔧 Setting up local configuration for 40-core server..."
echo ""

# 1. 로컬 miniwdl 설정 생성
if [ ! -f "config/miniwdl.local.cfg" ]; then
    cp config/miniwdl.cfg config/miniwdl.local.cfg
    echo "✅ Created config/miniwdl.local.cfg"
else
    echo "ℹ️  config/miniwdl.local.cfg already exists"
fi

# 2. 로컬 입력 템플릿 복사
if [ ! -f "sample.local.inputs.json" ]; then
    if [ -f "sample.local.inputs.json.example" ]; then
        cp sample.local.inputs.json.example sample.local.inputs.json
        echo "✅ Created sample.local.inputs.json"
    else
        echo "⚠️  sample.local.inputs.json.example not found"
    fi
else
    echo "ℹ️  sample.local.inputs.json already exists"
fi

echo ""
echo "📊 Server Specifications:"
echo "  - CPU Cores: $(nproc)"
echo "  - Total Memory: $(free -h | awk '/^Mem:/{print $2}')"
echo ""

# 3. Singularity 확인
if command -v singularity &> /dev/null; then
    echo "✅ Singularity found: $(singularity --version)"
else
    echo "⚠️  Singularity not found. Please install Singularity."
fi

# 4. GPU 확인
if command -v nvidia-smi &> /dev/null; then
    echo "✅ GPU detected:"
    nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
    echo ""
    echo "🎉 GPU found! Setting 'gpu: true' is HIGHLY RECOMMENDED"
    echo "   DeepVariant with GPU is 3-5x faster than CPU mode!"
else
    echo "ℹ️  No GPU detected (nvidia-smi not available)"
fi

echo ""
echo "📋 Next steps:"
echo "1. Edit sample.local.inputs.json with your sample information:"
echo "   - Update sample_id"
echo "   - Update hifi_reads path"
echo "   - Verify reference file paths"
echo ""
echo "2. Verify/adjust these parameters for 40-core server:"
echo "   - total_deepvariant_tasks: 32 (instead of 64)"
echo "   - deepvariant_tasks_per_shard: 8"
echo "   - gpu: true (HIGHLY RECOMMENDED - you have 2x RTX 2080 Ti!)"
echo ""
echo "3. Check GPU setup (recommended):"
echo "   chmod +x check_gpu_setup.sh"
echo "   ./check_gpu_setup.sh"
echo ""
echo "4. Run the workflow:"
echo "   miniwdl run --cfg config/miniwdl.local.cfg \\"
echo "     workflows/singleton.wdl \\"
echo "     -i sample.local.inputs.json"
echo ""
echo "5. Monitor resources during execution:"
echo "   # Terminal 1: GPU monitoring"
echo "   watch -n 5 nvidia-smi"
echo "   "
echo "   # Terminal 2: CPU/Memory monitoring"
echo "   watch -n 5 'free -h; echo; ps aux --sort=-%mem | head -10'"
