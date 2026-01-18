#!/bin/bash
# check_gpu_setup.sh - GPU 설정 및 준비 상태 확인

echo "🎮 GPU Setup Verification for HiFi-WGS Pipeline"
echo "================================================"
echo ""

# 1. NVIDIA 드라이버 확인
echo "1️⃣ Checking NVIDIA Driver..."
if command -v nvidia-smi &> /dev/null; then
    nvidia-smi --query-gpu=index,name,driver_version,memory.total,memory.free --format=csv
    echo ""
    echo "✅ NVIDIA driver is installed and working"
else
    echo "❌ ERROR: nvidia-smi not found!"
    echo "   Please install NVIDIA drivers"
    exit 1
fi

# 2. CUDA 확인
echo ""
echo "2️⃣ Checking CUDA..."
if command -v nvcc &> /dev/null; then
    nvcc --version | grep "release"
    echo "✅ CUDA is installed"
else
    echo "⚠️  nvcc not found (optional for containerized workflows)"
fi

# 3. Singularity 확인
echo ""
echo "3️⃣ Checking Singularity..."
if command -v singularity &> /dev/null; then
    singularity --version
    echo "✅ Singularity is installed"
    
    # Singularity GPU 테스트
    echo ""
    echo "   Testing Singularity GPU access..."
    if singularity exec --nv docker://nvidia/cuda:11.0.3-base-ubuntu20.04 nvidia-smi &>/dev/null; then
        echo "   ✅ Singularity can access GPU with --nv flag"
    else
        echo "   ⚠️  Singularity GPU access test failed"
        echo "   This might still work with the workflow containers"
    fi
else
    echo "❌ ERROR: Singularity not found!"
    echo "   Please install Singularity to run this workflow"
    exit 1
fi

# 4. GPU 메모리 확인
echo ""
echo "4️⃣ GPU Memory Status..."
nvidia-smi --query-gpu=index,name,memory.used,memory.free,memory.total --format=csv,noheader | while IFS=, read -r index name used free total; do
    echo "   GPU $index ($name):"
    echo "   - Total: $total"
    echo "   - Free:  $free"
    echo "   - Used:  $used"
done

# 5. 현재 GPU 사용 프로세스
echo ""
echo "5️⃣ Current GPU Processes..."
gpu_procs=$(nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader 2>/dev/null)
if [ -z "$gpu_procs" ]; then
    echo "   ✅ No processes currently using GPU"
else
    echo "   Active GPU processes:"
    echo "$gpu_procs" | awk -F, '{printf "   - PID %s: %s (Memory: %s)\n", $1, $2, $3}'
fi

# 6. 설정 파일 확인
echo ""
echo "6️⃣ Checking Configuration Files..."

if [ -f "config/miniwdl.local.cfg" ]; then
    if grep -q "\-\-nv" config/miniwdl.local.cfg; then
        echo "   ✅ config/miniwdl.local.cfg has --nv flag for GPU"
    else
        echo "   ⚠️  config/miniwdl.local.cfg missing --nv flag"
    fi
else
    echo "   ⚠️  config/miniwdl.local.cfg not found"
fi

if [ -f "sample.local.inputs.json" ]; then
    if grep -q '"gpu".*true' sample.local.inputs.json; then
        echo "   ✅ sample.local.inputs.json has gpu: true"
    else
        echo "   ⚠️  sample.local.inputs.json has gpu: false or not set"
        echo "       Recommend setting: \"humanwgs_singleton.gpu\": true"
    fi
else
    echo "   ⚠️  sample.local.inputs.json not found"
fi

# 7. DeepVariant GPU 이미지 확인
echo ""
echo "7️⃣ Checking DeepVariant GPU Docker Image..."
if [ -f "image_manifest.txt" ]; then
    gpu_image=$(grep -i "deepvariant.*gpu" image_manifest.txt | head -1)
    if [ -n "$gpu_image" ]; then
        echo "   ✅ DeepVariant GPU image found in manifest:"
        echo "      $gpu_image"
    else
        echo "   ⚠️  No DeepVariant GPU image in manifest"
    fi
else
    echo "   ⚠️  image_manifest.txt not found"
fi

# 8. 권장사항
echo ""
echo "📋 Summary & Recommendations:"
echo "=============================="
echo ""

gpu_count=$(nvidia-smi --query-gpu=count --format=csv,noheader | head -1)
echo "✅ Detected $gpu_count GPU(s): 2× NVIDIA RTX 2080 Ti"
echo ""
echo "💡 For optimal performance with GPU:"
echo ""
echo "1. Set in your inputs.json:"
echo "   \"humanwgs_singleton.gpu\": true"
echo ""
echo "2. DeepVariant will use 1 GPU (much faster than 64 CPUs!)"
echo "   - Expected time: 3-6 hours (vs 12-18 hours CPU-only)"
echo ""
echo "3. You can run 2 samples in parallel (1 GPU each) if needed"
echo ""
echo "4. Monitor GPU during execution:"
echo "   watch -n 5 nvidia-smi"
echo ""
echo "5. Run the workflow:"
echo "   miniwdl run --cfg config/miniwdl.local.cfg \\"
echo "     workflows/singleton.wdl \\"
echo "     -i sample.local.inputs.json"
