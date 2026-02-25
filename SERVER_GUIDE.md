# 서버 설정 가이드

## 서버 스펙

| 항목 | 사양 | 파이프라인 요구사항 |
|------|------|-------------------|
| CPU | 40 cores (Xeon E5-2640 v4) | 64 cores (권장) |
| RAM | 251 GB | 256 GB (권장) |
| GPU | 2× RTX 2080 Ti 11GB | 선택사항 |
| Storage | `/data_4tb` (4TB) | 샘플당 ~500GB |

> CPU·RAM이 권장치보다 약간 부족하나 GPU 사용 시 정상 실행 가능.
> **GPU 0번 고장** — GPU 1번만 사용 (`CUDA_VISIBLE_DEVICES=1`).

## miniwdl 설정

설정 파일: `config/miniwdl.local.cfg` (symlink → `~/.config/miniwdl.cfg`)

```ini
[scheduler]
container_backend = singularity
task_concurrency = 1        # GPU 과부하 방지. 늘리려면 2로 변경
fail_fast = false

[file_io]
allow_any_input = true      # ref_map 파일 사용에 필수
output_hardlinks = false    # cross-device link 오류 방지

[call_cache]
put = true
get = true
dir = "/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_call_cache"

[task_runtime]
command_shell = /bin/bash
defaults = {
    "maxRetries": 2,
    "docker": "ubuntu:20.04",
    "cpu": 16,
    "memory": "64G"
}

[singularity]
exe = ["/usr/bin/apptainer"]
run_options = ["--bind=/etc/hosts:/etc/hosts", "--nv"]
env = ["CUDA_VISIBLE_DEVICES=1", "TF_FORCE_GPU_ALLOW_GROWTH=true", "TF_GPU_THREAD_MODE=gpu_private"]
image_cache = "/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_singularity_cache"
```

### 설정 변경 시 주의사항

- `task_concurrency` 변경: 2로 높이면 메모리 부족 위험. DeepVariant는 ~64GB 사용
- `CUDA_VISIBLE_DEVICES`: GPU 0번 고장이므로 반드시 `1` 유지
- `output_hardlinks`: `/data_4tb`와 `/home`이 다른 파일시스템이므로 `false` 필수

## GPU 설정

### 상태 확인

```bash
nvidia-smi                    # GPU 전체 상태
nvidia-smi -l 5               # 5초마다 갱신
watch -n 2 nvidia-smi         # 별도 터미널에서 모니터링
```

### GPU 사용 활성화 (inputs.json)

```json
{
  "humanwgs_singleton.gpu": true
}
```

GPU 활성화 시 DeepVariant가 GPU 이미지(`google/deepvariant:1.9.0-gpu`)를 사용함.

### GPU vs CPU 성능 비교

| 단계 | CPU (40 cores) | GPU (RTX 2080 Ti) |
|------|---------------|-------------------|
| DeepVariant | ~10-12시간 | ~2-3시간 |
| 전체 파이프라인 | ~12-18시간 | ~4-8시간 |

### GPU 문제 해결

```bash
# GPU 인식 안 됨
nvidia-smi
ls /dev/nvidia*

# Singularity --nv 오류
singularity exec --nv docker://nvidia/cuda:11.0.3-base-ubuntu20.04 nvidia-smi

# GPU 메모리 점유 프로세스 확인
nvidia-smi | grep -E "MiB|PID"
```

## 설정 파일 관리 원칙

### Git에 커밋 O (공유 파일)

- `config/miniwdl.cfg` — 기본 설정 템플릿
- `*.inputs.json.example` — 입력 파일 예제
- `*.ref_map.*.template.tsv` — 레퍼런스 맵 템플릿
- `environment.yml` — conda 환경

### Git에 커밋 X (개인/로컬)

- `config/miniwdl.local.cfg` — 서버 특화 설정 (절대경로 포함)
- `*.inputs.json` — 실제 데이터 경로 포함
- `GRCm39.ref_map.tsv` — 로컬 경로 포함
- `hifi-wdl-resources/` — 대용량 데이터

### 개인 설정 파일 사용 패턴

```bash
# 예제 파일 복사 후 수정
cp sample.inputs.json.example my_sample.inputs.json
# my_sample.inputs.json은 .gitignore에 포함됨
```

## 디스크 공간 관리

```bash
# 전체 사용량
df -h /data_4tb

# 주요 디렉토리별 크기
du -sh /data_4tb/hifi-human-wgs-wdl-custom/miniwdl_call_cache/
du -sh /data_4tb/hifi-human-wgs-wdl-custom/miniwdl_singularity_cache/
du -sh /data_4tb/hifi-human-wgs-wdl-custom/batch_results/
du -sh /data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/
```

## 리소스 모니터링

```bash
htop                          # CPU/메모리
watch -n 2 nvidia-smi         # GPU
iostat -x 5                   # 디스크 I/O
ps aux | grep apptainer       # 실행 중인 컨테이너
```
