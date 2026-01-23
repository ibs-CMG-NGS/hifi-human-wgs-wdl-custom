# Cross-Device Link 에러 해결 방법

## 🔍 문제 원인

```
OSError: [Errno 18] Invalid cross-device link
```

**원인**: Call cache가 `/home` 디렉토리에 있고, 새 결과를 `/data_4tb`에 생성하려고 할 때, 
서로 다른 파일시스템 간에는 hard link를 만들 수 없어서 발생합니다.

- 기존 실행: `/home/ygkim/ngs-pipeline/hifi-human-wgs-wdl-custom/`
- 기존 cache: `/home/ygkim/ngs-pipeline/hifi-human-wgs-wdl-custom/miniwdl_call_cache/`
- 새로운 실행: `/data_4tb/hifi-human-wgs-wdl-custom/batch_results/`
- 새로운 cache 경로: `/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_call_cache/`

## ✅ 해결 방법

### 방법 1: Call Cache 비활성화 (가장 간단, 하지만 느림)

```bash
cd ~/ngs-pipeline/hifi-human-wgs-wdl-custom

# config/miniwdl.local.cfg 파일 백업
cp config/miniwdl.local.cfg config/miniwdl.local.cfg.backup

# call_cache 섹션 수정
cat > config/miniwdl.local.cfg << 'EOF'
[scheduler]
container_backend = singularity
task_concurrency = 1
fail_fast = false

[file_io]
allow_any_input = true
output_hardlinks = true

[call_cache]
# Cross-device link 문제 해결: 캐시 비활성화
put = false
get = false

[task_runtime]
command_shell = /bin/bash
gpu_enabled = true
defaults = {
        "maxRetries": 2,
        "docker": "ubuntu:20.04",
        "cpu": 16,
        "memory": "64G"
    }

[singularity]
exe = ["/usr/bin/apptainer"]
run_options = [
        "--bind", "/etc/hosts:/etc/hosts",
        "--nv"
    ]
env = ["CUDA_VISIBLE_DEVICES=1", "TF_FORCE_GPU_ALLOW_GROWTH=true", "TF_GPU_THREAD_MODE=gpu_private"]
image_cache = "/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_singularity_cache"
EOF
```

**단점**: 캐시를 사용하지 않아서 모든 단계를 처음부터 실행해야 합니다.

---

### 방법 2: Call Cache를 /data_4tb로 이동 (권장)

```bash
cd ~/ngs-pipeline/hifi-human-wgs-wdl-custom

# 1. 기존 캐시를 /data_4tb로 이동
sudo mv /home/ygkim/ngs-pipeline/hifi-human-wgs-wdl-custom/miniwdl_call_cache /data_4tb/hifi-human-wgs-wdl-custom/
sudo chown -R ygkim:ygkim /data_4tb/hifi-human-wgs-wdl-custom/miniwdl_call_cache

# 2. config 파일 수정
cat > config/miniwdl.local.cfg << 'EOF'
[scheduler]
container_backend = singularity
task_concurrency = 1
fail_fast = false

[file_io]
allow_any_input = true
output_hardlinks = true

[call_cache]
put = true
get = true
# 절대 경로로 /data_4tb 사용
dir = "/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_call_cache"

[task_runtime]
command_shell = /bin/bash
gpu_enabled = true
defaults = {
        "maxRetries": 2,
        "docker": "ubuntu:20.04",
        "cpu": 16,
        "memory": "64G"
    }

[singularity]
exe = ["/usr/bin/apptainer"]
run_options = [
        "--bind", "/etc/hosts:/etc/hosts",
        "--nv"
    ]
env = ["CUDA_VISIBLE_DEVICES=1", "TF_FORCE_GPU_ALLOW_GROWTH=true", "TF_GPU_THREAD_MODE=gpu_private"]
# 컨테이너 캐시도 절대 경로로
image_cache = "/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_singularity_cache"
EOF

# 3. 권한 확인
ls -ld /data_4tb/hifi-human-wgs-wdl-custom/miniwdl_call_cache
```

**장점**: 기존 캐시를 활용하여 이미 성공한 단계는 건너뜁니다.

---

### 방법 3: 새로운 Call Cache 디렉토리 생성 (절충안)

```bash
cd ~/ngs-pipeline/hifi-human-wgs-wdl-custom

# 1. config 파일 수정
cat > config/miniwdl.local.cfg << 'EOF'
[scheduler]
container_backend = singularity
task_concurrency = 1
fail_fast = false

[file_io]
allow_any_input = true
output_hardlinks = true

[call_cache]
put = true
get = true
# 새로운 캐시 디렉토리 사용
dir = "/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_call_cache_batch"

[task_runtime]
command_shell = /bin/bash
gpu_enabled = true
defaults = {
        "maxRetries": 2,
        "docker": "ubuntu:20.04",
        "cpu": 16,
        "memory": "64G"
    }

[singularity]
exe = ["/usr/bin/apptainer"]
run_options = [
        "--bind", "/etc/hosts:/etc/hosts",
        "--nv"
    ]
env = ["CUDA_VISIBLE_DEVICES=1", "TF_FORCE_GPU_ALLOW_GROWTH=true", "TF_GPU_THREAD_MODE=gpu_private"]
image_cache = "/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_singularity_cache"
EOF

# 2. 새 캐시 디렉토리 생성
mkdir -p /data_4tb/hifi-human-wgs-wdl-custom/miniwdl_call_cache_batch
```

**장점**: 기존 캐시는 그대로 두고 새로운 캐시 사용
**단점**: 첫 실행은 캐시가 없어서 처음부터 실행

---

## 🎯 권장 사항

**방법 2를 권장합니다** (Call Cache를 /data_4tb로 이동)

이유:
1. ✅ 기존 캐시를 활용하여 시간 절약
2. ✅ 모든 데이터를 /data_4tb에 집중
3. ✅ 디스크 공간 문제 해결

---

## 🚀 실행 순서 (방법 2 기준)

```bash
# 1. 현재 작업 중지 (이미 실패했으면 생략)

# 2. 캐시 이동
cd ~/ngs-pipeline/hifi-human-wgs-wdl-custom
sudo mv miniwdl_call_cache /data_4tb/hifi-human-wgs-wdl-custom/
sudo chown -R ygkim:ygkim /data_4tb/hifi-human-wgs-wdl-custom/miniwdl_call_cache

# 3. config 파일 업데이트 (위의 방법 2 스크립트 실행)

# 4. 실패한 결과 디렉토리 정리
rm -rf /data_4tb/hifi-human-wgs-wdl-custom/batch_results/LDK6217

# 5. 배치 재실행
export CUDA_VISIBLE_DEVICES=1
./batch_run_optimized.sh sequential
```

---

## 📊 확인

```bash
# 캐시 위치 확인
ls -ld /data_4tb/hifi-human-wgs-wdl-custom/miniwdl_call_cache

# config 파일 확인
grep "dir = " config/miniwdl.local.cfg

# 파일시스템 확인
df -h /home/ygkim /data_4tb
```

---

## ⚠️ 주의사항

1. **sudo 권한 필요**: 캐시 디렉토리 이동 시 sudo가 필요할 수 있습니다
2. **소유권 확인**: 이동 후 `chown`으로 소유권 변경 필수
3. **디스크 공간**: `/data_4tb`에 충분한 공간이 있는지 확인
   ```bash
   du -sh /home/ygkim/ngs-pipeline/hifi-human-wgs-wdl-custom/miniwdl_call_cache
   df -h /data_4tb
   ```

---

이 문서를 따라 수정하면 cross-device link 에러가 해결됩니다!
