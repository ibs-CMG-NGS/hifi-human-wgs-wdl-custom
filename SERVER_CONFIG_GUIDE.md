# HiFi-human-WGS-WDL 파이프라인 서버 스펙 점검 및 권장 설정

## 📊 현재 서버 스펙
- **CPU**: 40 cores (Intel Xeon E5-2640 v4 @ 2.40GHz)
  - 2 Sockets × 10 Cores × 2 Threads = 40 logical CPUs
- **메모리**: 251 GB (사용 가능: ~220 GB)
- **GPU**: 2× NVIDIA GeForce RTX 2080 Ti (11GB VRAM each) ✅
- **아키텍처**: x86_64
- **NUMA 노드**: 2

## ⚠️ 주요 확인 사항

### 1. **CPU/메모리 요구사항 비교**

| 항목 | 필요 스펙 | 현재 스펙 | 상태 |
|------|----------|----------|------|
| 최소 CPU | 64 cores | 40 cores | ⚠️ **부족** |
| 최소 메모리 | 256 GB | 251 GB | ⚠️ **거의 부족** |

**README.md에 명시된 요구사항:**
> "The most resource-heavy step in the workflow requires **64 cpu cores and 256 GB of RAM**."

### 2. **영향을 받는 주요 작업**

#### DeepVariant Call Variants (CPU 모드)
- **기본 설정**: `cpu = total_deepvariant_tasks` (보통 64)
- **메모리**: `total_deepvariant_tasks * 4 GB` (256 GB)
- **현재 문제**: 40 코어로는 기본 병렬 작업 수행 불가

#### pbmm2 Alignment
- **기본 설정**: 스레드 수에 따라 동적 조정
- **영향**: 상대적으로 적음

## 🔧 필수 설정 변경 사항

### 1. **miniwdl.cfg 수정** (로컬 실행 시)

현재 `/config/miniwdl.cfg`를 복사하여 로컬 설정 생성:

```bash
cp config/miniwdl.cfg config/miniwdl.local.cfg
```

**수정 내용:**
```ini
[scheduler]
container_backend = singularity
# 서버의 40 코어를 고려하여 동시 실행 작업 수 제한
task_concurrency = 2
fail_fast = false

[file_io]
allow_any_input = true

[call_cache]
put = true
get = true
dir = "$PWD/miniwdl_call_cache"

[task_runtime]
command_shell = /bin/bash
defaults = {
        "maxRetries": 2,
        "docker": "ubuntu:20.04",
        "cpu": 16,
        "memory": "64G"
    }

[singularity]
exe = ["/usr/bin/singularity"]
run_options = [
        "--containall",
        "--nv"
    ]
image_cache = "$PWD/miniwdl_singularity_cache"
```

**실행 시:**
```bash
miniwdl run --cfg config/miniwdl.local.cfg workflows/singleton.wdl -i my_sample.inputs.json
```

### 2. **Inputs JSON 파일 수정**

`sample1.inputs.json` 또는 새로운 입력 파일에 다음 설정 추가:

```json
{
  "humanwgs_singleton.sample_id": "KTY9537",
  "humanwgs_singleton.sex": "MALE",
  "humanwgs_singleton.hifi_reads": [
    "/home/ygkim/ngs_pipeline/HiFi-human-WGS-WDL/data/m84285_260108_082608_s1.hifi_reads.bc2016.bam"
  ],
  "humanwgs_singleton.ref_map_file": "/home/ygkim/ngs_pipeline/HiFi-human-WGS-WDL/hifi-wdl-resources/hifi-wdl-resources-v3.1.0/GRCh38.ref_map.v3p1p0.template.tsv",
  "humanwgs_singleton.tertiary_map_file": "/home/ygkim/ngs_pipeline/HiFi-human-WGS-WDL/hifi-wdl-resources/hifi-wdl-resources-v3.1.0/GRCh38.tertiary_map.v3p1p0.template.tsv",
  "humanwgs_singleton.backend": "HPC",
  "humanwgs_singleton.preemptible": false,
  
  "humanwgs_singleton.total_deepvariant_tasks": 32,
  "humanwgs_singleton.deepvariant_tasks_per_shard": 8,
  "humanwgs_singleton.max_reads_per_alignment_chunk": 100000000
}
```

**주요 변경 파라미터:**
- `total_deepvariant_tasks`: 64 → **32** (40코어에 맞춤)
- `deepvariant_tasks_per_shard`: 기본값 → **8** (메모리 사용 최적화)
- `preemptible`: true → **false** (로컬 서버에는 해당 없음)

### 3. **SLURM 사용 시 (HPC 백엔드)**

`backends/hpc/miniwdl.cfg` 참조하여 설정:

```ini
[scheduler]
container_backend = slurm_singularity
task_concurrency = 50
fail_fast = false

[slurm]
# 파티션과 리소스 제한 설정
extra_args = "--partition compute --comment 'HiFi-WGS' --cpus-per-task=16 --mem=64G"
```

## 💡 추가 최적화 권장사항

### 1. **GPU 사용 (강력 권장!)** ⭐
서버에 **NVIDIA RTX 2080 Ti 2개**가 설치되어 있습니다!

**GPU 상태 확인:**
```bash
nvidia-smi
```

**DeepVariant는 GPU를 사용하면 성능이 크게 향상됩니다:**
- CPU 모드: 64 코어 필요, 8-12시간 소요
- **GPU 모드: 1개 GPU로 충분, 2-4시간 소요** 🚀

**inputs.json에서 GPU 활성화 (필수!):**
```json
{
  "humanwgs_singleton.gpu": true
}
```

**주의사항:**
- GPU 드라이버 확인: `nvidia-smi`가 정상 작동해야 함
- CUDA 호환성 확인 (RTX 2080 Ti는 CUDA 7.5 지원)
- Singularity는 `--nv` 옵션으로 GPU 접근 (이미 설정됨)

### 2. **메모리 모니터링**
파이프라인 실행 중 메모리 사용량 모니터링:

```bash
watch -n 5 free -h
```

### 3. **워크플로우 단계별 실행**
리소스가 부족하면 단계별로 나누어 실행 고려:
- Upstream (alignment + variant calling)
- Downstream (phasing + analysis)

### 4. **데이터 청킹 조정**
`max_reads_per_alignment_chunk`를 조정하여 메모리 사용 분산:
- 기본값: 매우 큼 (전체를 한번에)
- 권장값: 50M - 100M reads per chunk

## 📝 설정 파일 생성 스크립트

로컬 환경에 맞는 설정 파일 자동 생성:

```bash
#!/bin/bash
# setup_local_config.sh

# 로컬 miniwdl 설정 생성
cat > config/miniwdl.local.cfg << 'EOF'
[scheduler]
container_backend = singularity
task_concurrency = 2
fail_fast = false

[file_io]
allow_any_input = true

[call_cache]
put = true
get = true
dir = "$PWD/miniwdl_call_cache"

[task_runtime]
command_shell = /bin/bash
defaults = {
        "maxRetries": 2,
        "docker": "ubuntu:20.04",
        "cpu": 16,
        "memory": "64G"
    }

[singularity]
exe = ["/usr/bin/singularity"]
run_options = ["--containall", "--nv"]
image_cache = "$PWD/miniwdl_singularity_cache"
EOF

echo "✅ Local miniwdl configuration created: config/miniwdl.local.cfg"

# 로컬 입력 템플릿 생성
cat > sample.local.inputs.json << 'EOF'
{
  "humanwgs_singleton.sample_id": "YOUR_SAMPLE_ID",
  "humanwgs_singleton.sex": "MALE",
  "humanwgs_singleton.hifi_reads": [
    "/path/to/your/data/sample.hifi_reads.bam"
  ],
  "humanwgs_singleton.ref_map_file": "/home/ygkim/ngs_pipeline/HiFi-human-WGS-WDL/hifi-wdl-resources/hifi-wdl-resources-v3.1.0/GRCh38.ref_map.v3p1p0.template.tsv",
  "humanwgs_singleton.tertiary_map_file": "/home/ygkim/ngs_pipeline/HiFi-human-WGS-WDL/hifi-wdl-resources/hifi-wdl-resources-v3.1.0/GRCh38.tertiary_map.v3p1p0.template.tsv",
  "humanwgs_singleton.backend": "HPC",
  "humanwgs_singleton.preemptible": false,
  "humanwgs_singleton.total_deepvariant_tasks": 32,
  "humanwgs_singleton.deepvariant_tasks_per_shard": 8,
  "humanwgs_singleton.max_reads_per_alignment_chunk": 100000000,
  "humanwgs_singleton.gpu": false
}
EOF

echo "✅ Local input template created: sample.local.inputs.json"
echo ""
echo "📋 Next steps:"
echo "1. Edit sample.local.inputs.json with your sample information"
echo "2. Run: miniwdl run --cfg config/miniwdl.local.cfg workflows/singleton.wdl -i sample.local.inputs.json"
```

## 🚀 실행 전 체크리스트

- [ ] **GPU 확인**: `./check_gpu_setup.sh` 실행하여 GPU 준비 상태 확인
- [ ] `config/miniwdl.local.cfg` 생성 및 `--nv` 플래그 확인
- [ ] inputs.json에 **`"humanwgs_singleton.gpu": true`** 설정 (필수!)
- [ ] inputs.json에 `total_deepvariant_tasks: 32` 설정
- [ ] inputs.json에 `deepvariant_tasks_per_shard: 8` 설정
- [ ] Singularity 이미지 캐시 확인
- [ ] 참조 데이터 다운로드 확인 (`hifi-wdl-resources/`)
- [ ] 충분한 디스크 공간 확인 (최소 500GB 권장)
- [ ] `nvidia-smi`로 GPU 사용 가능 확인

## ⚡ 예상 성능

**40코어 + 2×RTX 2080 Ti 서버에서의 예상 실행 시간:**
- **원본 권장 스펙 (64코어, CPU only)**: ~8-12시간
- **현재 스펙 (40코어, CPU only)**: ~12-18시간
- **현재 스펙 (40코어 + GPU 1개 사용)**: ~3-6시간 ⚡🚀
- **단일 샘플 WGS 30x coverage 기준**

**GPU 사용의 장점:**
- ✅ DeepVariant 단계가 CPU 64코어 대신 GPU 1개로 처리
- ✅ CPU 리소스를 다른 작업에 할당 가능
- ✅ 메모리 부담 대폭 감소 (256GB → ~64GB)
- ✅ 전체 파이프라인 실행 시간 50-70% 단축

**강력 권장: GPU 모드 사용!** 메모리 부족 문제도 해결됩니다.
