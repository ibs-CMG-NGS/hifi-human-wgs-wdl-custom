# 환경 설정 가이드

## 서버 환경 (현재)

| 항목 | 사양 |
|------|------|
| CPU | 40 cores (Intel Xeon E5-2640 v4) |
| RAM | 251 GB |
| GPU | 2× NVIDIA RTX 2080 Ti 11GB — **GPU 0 고장, GPU 1만 사용** |
| Container | Apptainer 1.4.5 (`/usr/bin/singularity`) |
| 파이프라인 | `/data_4tb/hifi-human-wgs-wdl-custom/` |

## Step 1. Conda 환경

```bash
conda env create -f environment.yml
conda activate hifi-human-wgs
miniwdl --version
```

bgzip/tabix가 없으면 추가 설치:

```bash
conda install -c bioconda htslib -y
```

## Step 2. miniwdl 설정 활성화

```bash
mkdir -p ~/.config
ln -sf /data_4tb/hifi-human-wgs-wdl-custom/config/miniwdl.local.cfg ~/.config/miniwdl.cfg
```

현재 `config/miniwdl.local.cfg` 주요 설정:

| 항목 | 값 | 이유 |
|------|-----|------|
| `container_backend` | `singularity` | Apptainer 사용 |
| `task_concurrency` | `1` | GPU 과부하 방지 |
| `output_hardlinks` | `false` | cross-device link 오류 방지 |
| `call_cache dir` | `/data_4tb/.../miniwdl_call_cache` | 절대경로, 재실행 캐시 |
| `image_cache` | `/data_4tb/.../miniwdl_singularity_cache` | 이미지 캐시 |
| `exe` | `/usr/bin/apptainer` | Apptainer 실행 파일 |
| `CUDA_VISIBLE_DEVICES` | `1` | GPU 1번 고정 |

## Step 3. 레퍼런스 데이터

### GRCh38 (Human) — 이미 설치됨

```
hifi-wdl-resources/hifi-wdl-resources-v3.1.0/GRCh38/
├── human_GRCh38_no_alt_analysis_set.fasta(.fai)
├── trgt/adotto_strchive_20250827.hg38.bed.gz
├── sawfish/
├── pharmcat/
└── ...
```

레퍼런스 맵: `backends/hpc/GRCh38.ref_map.v3p1p0.hpc.tsv`

### GRCm39 (Mouse) — 이미 설치됨

```
hifi-wdl-resources/GRCm39/
├── mouse_GRCm39.fasta(.fai)          ← 2.6GB
├── trgt_mm39.bed.gz(.tbi)            ← UCSC simpleRepeat 기반, 1.6M loci
├── sawfish_exclude_mm39.bed.gz(.tbi) ← 빈 파일 (exclude 없음)
├── expected_cn.mm39.XY.bed           ← 수컷 예상 copy number
├── expected_cn.mm39.XX.bed           ← 암컷 예상 copy number
└── cpgIslandExt.sorted.mm39.tsv      ← UCSC CpG islands
```

레퍼런스 맵: `GRCm39.ref_map.tsv`

### Human 레퍼런스 재다운로드 (필요 시)

```bash
wget https://zenodo.org/record/17086906/files/hifi-wdl-resources-v3.1.0.tar
tar -xvf hifi-wdl-resources-v3.1.0.tar -C hifi-wdl-resources/
```

## Step 4. Singularity 이미지 캐시

`miniwdl_singularity_cache/`에 이미 다운로드됨. 신규 이미지 필요 시:

```bash
bash scripts/prefetch_images.sh
```

## Step 5. 실행 확인

```bash
conda activate hifi-human-wgs
cd /data_4tb/hifi-human-wgs-wdl-custom

miniwdl run workflows/singleton.wdl \
  --input BioSample24.inputs.json \
  --dir batch_results/ \
  --verbose
```

## Mouse 데이터 분석 시 주의사항

이 파이프라인은 human 전용으로 설계되었으나 mouse에서도 핵심 단계는 동작함:

| 도구 | Mouse 적용 | 비고 |
|------|-----------|------|
| pbmm2, DeepVariant, Sawfish, HiPhase | ✅ 정상 | Species-agnostic |
| TRGT | ✅ 부분적 | Mouse catalog 사용 (UCSC 기반) |
| PharmCAT, PBstarPhase | ⚠️ 빈 결과 | Human PGx 전용 |
| Paraphase | ⚠️ 빈 결과 | Human HLA 전용 |
| Tertiary 분석 | ❌ 미적용 | Human population DB 필요 |

Mouse 실행 시 `tertiary_map_file`과 `phenotypes`는 inputs.json에서 제외.
