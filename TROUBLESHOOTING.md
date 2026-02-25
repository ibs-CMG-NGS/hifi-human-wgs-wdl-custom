# 문제 해결 가이드

## Cross-Device Link 오류

### 증상

```
OSError: [Errno 18] Invalid cross-device link
```

### 원인

Call cache(`/home` 파일시스템)와 결과 디렉토리(`/data_4tb` 파일시스템)가 다른 파티션에 있을 때,
miniwdl이 hard link를 만들려다 실패.

### 해결 (현재 설정에 반영됨)

`config/miniwdl.local.cfg`에 두 가지 설정이 적용되어 있음:

```ini
[file_io]
output_hardlinks = false          # hard link 비활성화

[call_cache]
dir = "/data_4tb/.../miniwdl_call_cache"  # /data_4tb에 캐시 위치
```

오류 재발 시 위 두 설정을 확인.

### 기존 캐시가 /home에 있는 경우

```bash
sudo mv /home/ygkim/.../miniwdl_call_cache /data_4tb/hifi-human-wgs-wdl-custom/
sudo chown -R ygkim:ygkim /data_4tb/hifi-human-wgs-wdl-custom/miniwdl_call_cache
```

---

## Singularity/Apptainer 관련 오류

### "singularity: command not found"

```bash
which singularity && which apptainer
grep "exe" ~/.config/miniwdl.cfg
```

### "--nv 옵션 오류" (GPU 관련)

```bash
apptainer --version                # 1.0+ 필요
nvidia-smi                         # GPU 드라이버 확인
singularity exec --nv docker://nvidia/cuda:11.0.3-base-ubuntu20.04 nvidia-smi
```

### "image not found"

```bash
ls /data_4tb/hifi-human-wgs-wdl-custom/miniwdl_singularity_cache/
bash scripts/prefetch_images.sh    # 이미지 재다운로드
```

---

## miniwdl 실행 오류

### "allow_any_input" 오류

```ini
# ~/.config/miniwdl.cfg에 추가
[file_io]
allow_any_input = true
```

### 태스크 실패 후 재실행

`fail_fast = false` 설정 시 다른 태스크는 계속 실행됨.
실패한 태스크만 재실행하려면 해당 call_cache 삭제 후 재실행:

```bash
# 실패한 태스크 캐시 확인 및 삭제
ls miniwdl_call_cache/
rm -rf miniwdl_call_cache/<failed_task_hash>

# 동일 명령으로 재실행 (성공 태스크는 캐시 사용)
miniwdl run workflows/singleton.wdl --input sample.inputs.json --dir batch_results/
```

### 메모리 부족 오류

```bash
free -h
# config/miniwdl.local.cfg에서:
# task_concurrency = 1  (이미 설정됨)
```

---

## GPU 오류

### GPU 인식 안 됨

```bash
nvidia-smi
grep "CUDA_VISIBLE_DEVICES" ~/.config/miniwdl.cfg
# 현재 서버: GPU 0 고장, GPU 1만 사용 → CUDA_VISIBLE_DEVICES=1
```

### GPU 메모리 부족 (OOM)

```bash
nvidia-smi                         # 점유 프로세스 확인
# inputs.json에서 gpu: false로 변경 후 CPU 모드 실행
```

---

## Sawfish 오류

### expected copy number genome regions record has 4 columns where at least 5 are required

```
thread 'main' panicked at src/genome_regions.rs:226:9:
expected copy number genome regions record has 4 columns where at least 5 are required
```

#### 원인

`expected_cn.bed` 파일 형식 오류. Sawfish는 5컬럼을 요구함:

| chrom | start | end | **name** | copy_number |
|-------|-------|-----|----------|-------------|
| chr1  | 0 | 195154279 | chr1 | 2 |

`name` 컬럼(4번째)이 없으면 crash.

#### 해결

```bash
cd /data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39

# name 컬럼(염색체명) 삽입
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $1, $4}' expected_cn.mm39.XX.bed > tmp.bed && mv tmp.bed expected_cn.mm39.XX.bed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $1, $4}' expected_cn.mm39.XY.bed > tmp.bed && mv tmp.bed expected_cn.mm39.XY.bed

# 확인 (5컬럼인지)
head -3 expected_cn.mm39.XX.bed
```

수정 후 파이프라인 재실행 (failed 태스크는 자동으로 캐시 미적용).

---

## 레퍼런스 파일 오류

### ref_map TSV 경로 오류

```bash
# 각 경로 존재 확인
while IFS=$'\t' read -r key val; do
  [ -f "$val" ] && echo "OK: $key" || echo "MISSING: $key -> $val"
done < GRCm39.ref_map.tsv
```

### Mouse 분석에서 PharmCAT/Paraphase 에러

Mouse genome에서 human 전용 도구 실행 시 빈 결과 또는 에러 발생 — 정상적인 동작임.
파이프라인 중단 시:

1. `fail_fast = false` 설정 확인
2. 에러 로그 확인: `batch_results/*/call-*/stderr`
3. 해당 태스크 call_cache 삭제 후 재시도
