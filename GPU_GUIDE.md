# GPU 활용 가이드 - HiFi WGS Pipeline

## 🎮 서버 GPU 정보
- **GPU 개수**: 2개
- **모델**: NVIDIA GeForce RTX 2080 Ti
- **VRAM**: 11GB per GPU
- **CUDA Compute Capability**: 7.5

## ⚡ GPU 사용의 이점

### CPU vs GPU 모드 비교

| 항목 | CPU 모드 | GPU 모드 |
|------|----------|----------|
| **필요 CPU 코어** | 64 cores | 8-16 cores |
| **필요 메모리** | 256 GB | 64 GB |
| **DeepVariant 시간** | 8-12시간 | **2-4시간** ⚡ |
| **전체 파이프라인** | 12-18시간 | **4-8시간** ⚡ |
| **서버 적합성** | ⚠️ 부족 (40 cores) | ✅ **완벽** |

### 핵심 포인트
- ✅ **CPU 부족 문제 해결**: 64 코어 대신 GPU 1개로 해결
- ✅ **메모리 절약**: 256GB → 64GB로 요구사항 대폭 감소
- ✅ **실행 시간 단축**: 50-70% 시간 절약
- ✅ **동시 실행 가능**: GPU 2개 → 샘플 2개 동시 처리 가능

## 🔧 GPU 활성화 방법

### 1. GPU 상태 확인
```bash
# GPU 확인
nvidia-smi

# 상세 정보
nvidia-smi --query-gpu=index,name,driver_version,memory.total,memory.free --format=csv

# 실시간 모니터링
watch -n 2 nvidia-smi
```

### 2. 설정 파일에서 GPU 활성화

**`sample.local.inputs.json` 또는 `sample1.inputs.json`:**
```json
{
  "humanwgs_singleton.gpu": true,
  "humanwgs_singleton.total_deepvariant_tasks": 32,
  "humanwgs_singleton.deepvariant_tasks_per_shard": 8
}
```

**`config/miniwdl.local.cfg` 확인 (이미 설정됨):**
```ini
[singularity]
run_options = [
        "--containall",
        "--nv"    # <- GPU 접근 활성화
    ]
```

### 3. GPU 준비 상태 확인
```bash
chmod +x check_gpu_setup.sh
./check_gpu_setup.sh
```

## 🚀 실행 예제

### 기본 실행 (GPU 사용)
```bash
# 1. 설정 파일 생성
cp sample.local.inputs.json.example my_sample.inputs.json

# 2. 경로 및 gpu: true 설정
nano my_sample.inputs.json

# 3. GPU 체크
./check_gpu_setup.sh

# 4. 실행
miniwdl run --cfg config/miniwdl.local.cfg \
  workflows/singleton.wdl \
  -i my_sample.inputs.json

# 5. 모니터링 (별도 터미널)
watch -n 2 nvidia-smi
```

### 2개 샘플 동시 실행 (GPU 2개 활용)
```bash
# Terminal 1: Sample 1 (GPU 0 사용)
CUDA_VISIBLE_DEVICES=0 miniwdl run --cfg config/miniwdl.local.cfg \
  workflows/singleton.wdl \
  -i sample1.inputs.json

# Terminal 2: Sample 2 (GPU 1 사용)  
CUDA_VISIBLE_DEVICES=1 miniwdl run --cfg config/miniwdl.local.cfg \
  workflows/singleton.wdl \
  -i sample2.inputs.json
```

## 📊 GPU 모니터링

### 실시간 모니터링
```bash
# 기본 모니터링
watch -n 2 nvidia-smi

# GPU 사용률과 메모리만 표시
watch -n 2 'nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv'

# 프로세스 상세 정보
watch -n 2 'nvidia-smi pmon -c 1'

# GPU 온도 모니터링
watch -n 2 'nvidia-smi --query-gpu=temperature.gpu,power.draw,power.limit --format=csv'
```

### 로깅
```bash
# GPU 사용률 로그 저장
nvidia-smi --query-gpu=timestamp,utilization.gpu,memory.used,memory.total \
  --format=csv -l 10 > gpu_usage.log &

# 백그라운드 로깅 중지
pkill -f "nvidia-smi.*loop"
```

## 🔍 문제 해결

### GPU가 인식되지 않는 경우
```bash
# NVIDIA 드라이버 확인
nvidia-smi

# Singularity GPU 테스트
singularity exec --nv \
  docker://nvidia/cuda:11.0.3-base-ubuntu20.04 \
  nvidia-smi

# 드라이버 재로드 (필요 시)
sudo rmmod nvidia_uvm
sudo modprobe nvidia_uvm
```

### GPU 메모리 부족 오류
```bash
# 실행 중인 GPU 프로세스 확인
nvidia-smi

# 좀비 프로세스 종료 (주의!)
nvidia-smi | grep python | awk '{print $5}' | xargs -r kill -9
```

### Singularity --nv 플래그 오류
```bash
# Singularity 버전 확인 (3.5+ 필요)
singularity --version

# 설정 확인
grep "run_options" config/miniwdl.local.cfg

# 수동 테스트
singularity exec --nv docker://nvidia/cuda:11.0.3-base-ubuntu20.04 nvidia-smi
```

## 💡 GPU 최적화 팁

### 1. GPU 메모리 최적화
```json
// inputs.json
{
  "humanwgs_singleton.deepvariant_tasks_per_shard": 8,  // GPU 메모리에 맞춰 조정
  "humanwgs_singleton.total_deepvariant_tasks": 32      // 전체 작업 수
}
```

### 2. GPU 전력 제한 (과열 방지)
```bash
# RTX 2080 Ti의 전력 제한 설정 (선택사항)
sudo nvidia-smi -pl 250  # 250W로 제한 (기본 280W)
```

### 3. GPU 클럭 고정 (일관된 성능)
```bash
# Persistence mode 활성화
sudo nvidia-smi -pm 1

# 성능 확인
nvidia-smi -q -d PERFORMANCE
```

## 📈 성능 벤치마크 (예상)

### DeepVariant 단계만 (30x WGS)
- **CPU (64 cores)**: ~6-8시간
- **CPU (40 cores)**: ~10-12시간
- **GPU (1x RTX 2080 Ti)**: ~2-3시간 ⚡

### 전체 파이프라인 (30x WGS)
- **CPU only (40 cores)**: 12-18시간
- **GPU (1x RTX 2080 Ti)**: **4-8시간** ⚡
- **2 samples parallel (2 GPUs)**: 각각 4-8시간

## ✅ 권장 설정 요약

**40 cores + 2× RTX 2080 Ti 서버에 최적화된 설정:**

```json
{
  "humanwgs_singleton.sample_id": "YOUR_SAMPLE",
  "humanwgs_singleton.sex": "MALE",
  "humanwgs_singleton.hifi_reads": ["path/to/reads.bam"],
  "humanwgs_singleton.ref_map_file": "path/to/ref_map.tsv",
  "humanwgs_singleton.tertiary_map_file": "path/to/tertiary_map.tsv",
  "humanwgs_singleton.backend": "HPC",
  "humanwgs_singleton.preemptible": false,
  
  "humanwgs_singleton.gpu": true,                           // ⭐ 필수!
  "humanwgs_singleton.total_deepvariant_tasks": 32,         // 40 cores에 맞춤
  "humanwgs_singleton.deepvariant_tasks_per_shard": 8,      // GPU 메모리 최적화
  "humanwgs_singleton.max_reads_per_alignment_chunk": 100000000
}
```

**이 설정으로 CPU 부족 문제가 완전히 해결되고, 실행 시간도 크게 단축됩니다!** 🚀
