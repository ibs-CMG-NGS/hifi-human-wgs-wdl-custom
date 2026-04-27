# 환경 설정 가이드

## 사전 요구사항

아래 도구들이 시스템에 설치되어 있어야 합니다. 각 명령으로 설치 여부를 확인합니다. 미설치 항목이 있으면 [INSTALL.md](./INSTALL.md)를 먼저 참고합니다.

| 도구 | 확인 명령 | 용도 |
|------|----------|------|
| Conda (Miniconda/Anaconda) | `conda --version` | Python 환경 및 패키지 관리 |
| Apptainer (Singularity) | `apptainer --version` | 컨테이너 실행 (1.0 이상) |
| samtools | `samtools --version` | FASTA 인덱스 생성 |
| git | `git --version` | 저장소 관리 |
| NVIDIA 드라이버 | `nvidia-smi` | GPU 사용 시 필수 |

## 서버 환경 (현재)

| 항목 | 사양 |
|------|------|
| CPU | 40 cores (Intel Xeon E5-2640 v4) |
| RAM | 251 GB |
| GPU | 2× NVIDIA RTX 2080 Ti 11GB — **GPU 0 고장, GPU 1만 사용** |
| Container | Apptainer 1.4.5 (`/usr/bin/singularity`) |
| 파이프라인 | `/data_4tb/hifi-human-wgs-wdl-custom/` |

## Step 0. 저장소 클론

```bash
git clone https://github.com/ibs-CMG-NGS/hifi-human-wgs-wdl-custom.git /data_4tb/hifi-human-wgs-wdl-custom
cd /data_4tb/hifi-human-wgs-wdl-custom
```

## Step 1. Conda 환경

```bash
conda env create -f environment.yml
conda activate hifi-human-wgs
miniwdl --version
```

## Step 2. miniwdl 설정

`config/miniwdl.local.cfg`는 이 서버에 맞는 실행 설정 파일입니다. Git에서 제외되므로 새 서버에서는 직접 생성합니다.

```bash
cat > config/miniwdl.local.cfg << 'EOF'
[scheduler]
container_backend = singularity
task_concurrency = 1
fail_fast = false

[file_io]
allow_any_input = true
output_hardlinks = false

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
EOF
```

> 다른 서버에 설치할 경우 `dir`, `image_cache`의 절대경로와 `CUDA_VISIBLE_DEVICES` 값을 환경에 맞게 수정합니다.

그 다음 miniwdl이 이 설정을 자동으로 읽도록 symlink를 생성합니다:

```bash
mkdir -p ~/.config
ln -sf /data_4tb/hifi-human-wgs-wdl-custom/config/miniwdl.local.cfg ~/.config/miniwdl.cfg
```

주요 설정 항목:

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

### GRCh38 (Human)

```bash
mkdir -p hifi-wdl-resources
wget https://zenodo.org/record/17086906/files/hifi-wdl-resources-v3.1.0.tar
tar -xvf hifi-wdl-resources-v3.1.0.tar -C hifi-wdl-resources/
rm hifi-wdl-resources-v3.1.0.tar
```

다운로드 후 구조:

```
hifi-wdl-resources/hifi-wdl-resources-v3.1.0/GRCh38/
├── human_GRCh38_no_alt_analysis_set.fasta(.fai)
├── trgt/adotto_strchive_20250827.hg38.bed.gz
├── sawfish/
├── pharmcat/
└── ...
```

레퍼런스 맵: `backends/hpc/GRCh38.ref_map.v3p1p0.hpc.tsv` (Git에 포함)

### GRCm39 (Mouse)

UCSC에서 원본 파일을 다운로드한 후 파이프라인 포맷으로 변환합니다.

```bash
mkdir -p hifi-wdl-resources/GRCm39
cd hifi-wdl-resources/GRCm39
```

#### 1) FASTA

UCSC mm39 전체 게놈 서열을 다운로드하고 samtools로 인덱스를 생성합니다. (~2.6 GB)

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
gunzip mm39.fa.gz
mv mm39.fa mouse_GRCm39.fasta
samtools faidx mouse_GRCm39.fasta
```

#### 2) TRGT 탠덤반복 카탈로그

TRGT는 탠덤반복(tandem repeat) 영역을 genotyping하기 위해 반복 카탈로그 BED 파일이 필요합니다. UCSC Genome Browser의 `simpleRepeat` 트랙(TRF 기반, ~1.6M loci)을 TRGT 입력 포맷(`chrom start end id motif (motif)n`)으로 변환합니다.

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/simpleRepeat.txt.gz

zcat simpleRepeat.txt.gz \
  | awk 'BEGIN{OFS="\t"} {print $2, $3, $4, $2"_"$3"_"$4, $17, "("$17")n"}' \
  | bgzip > trgt_mm39.bed.gz
tabix -p bed trgt_mm39.bed.gz
```

#### 3) Sawfish exclude BED

Sawfish SV/CNV 분석에서 제외할 문제 영역(반복 서열이 밀집한 구간 등)을 지정하는 파일입니다. GRCm39에는 별도로 정의된 exclude 영역이 없으므로 빈 파일을 생성합니다.

```bash
echo -n | bgzip > sawfish_exclude_mm39.bed.gz
tabix -p bed sawfish_exclude_mm39.bed.gz
```

#### 4) Expected copy number BED

Sawfish CNV 분석 시 각 염색체의 기준 copy number를 지정합니다. `.fai` 인덱스에서 염색체 목록과 길이를 읽어 sex에 따라 XY/XX 두 파일을 생성합니다. chrX·chrY는 XY에서 1, XX에서는 chrY 계열 전체를 제외합니다.

```bash
# XY (수컷): chrX=1, chrY=1, chrM=1, 나머지=2
awk 'BEGIN{OFS="\t"} {
    if ($1=="chrM") cn=1
    else if ($1=="chrX" || $1=="chrY") cn=1
    else cn=2
    print $1, 0, $2, $1, cn
}' mouse_GRCm39.fasta.fai > expected_cn.mm39.XY.bed

# XX (암컷): chrY 계열 제외, chrX=2, chrM=1, 나머지=2
awk 'BEGIN{OFS="\t"} {
    if ($1 ~ /^chrY/) next
    if ($1=="chrM") cn=1
    else cn=2
    print $1, 0, $2, $1, cn
}' mouse_GRCm39.fasta.fai > expected_cn.mm39.XX.bed
```

#### 5) CpG island TSV

MethBat 메틸화 프로파일링에 사용하는 CpG island 영역 파일입니다. UCSC `cpgIslandExt` 트랙에서 좌표를 추출하고, 정렬 후 순번 기반 레이블(`CpG_1`, `CpG_2`, …)을 붙입니다.

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/cpgIslandExt.txt.gz

zcat cpgIslandExt.txt.gz \
  | awk 'BEGIN{OFS="\t"} {print $2, $3, $4}' \
  | sort -k1,1 -k2,2n \
  | awk 'BEGIN{OFS="\t"; print "chrom\tstart\tend\tcpg_label"} {print $1, $2, $3, "CpG_"NR}' \
  > cpgIslandExt.sorted.mm39.tsv
```

#### 6) ref_map.tsv 생성

위에서 생성한 파일들의 절대경로를 파이프라인에 전달하는 TSV입니다. 절대경로를 포함하므로 Git에서 제외되며, 서버별로 다시 생성해야 합니다. 파이프라인 루트에서 실행합니다:

```bash
cd /data_4tb/hifi-human-wgs-wdl-custom
BASE="$(pwd)/hifi-wdl-resources/GRCm39"

cat > GRCm39.ref_map.tsv << EOF
name	GRCm39
fasta	${BASE}/mouse_GRCm39.fasta
fasta_index	${BASE}/mouse_GRCm39.fasta.fai
trgt_tandem_repeat_bed	${BASE}/trgt_mm39.bed.gz
sawfish_exclude_bed	${BASE}/sawfish_exclude_mm39.bed.gz
sawfish_exclude_bed_index	${BASE}/sawfish_exclude_mm39.bed.gz.tbi
sawfish_expected_bed_male	${BASE}/expected_cn.mm39.XY.bed
sawfish_expected_bed_female	${BASE}/expected_cn.mm39.XX.bed
methbat_region_tsv	${BASE}/cpgIslandExt.sorted.mm39.tsv
EOF
```

> `pharmcat_positions_vcf` 항목은 포함하지 않음 — Mouse 분석 시 pbstarphase가 human PGx DB로 mouse genome 조회를 시도해 범위 초과 panic을 일으킵니다 ([상세](TROUBLESHOOTING.md#mouse-분석에서-pbstarphase-crash)).

생성 후 경로 검증:

```bash
while IFS=$'\t' read -r key val; do
  [ -f "$val" ] && echo "OK: $key" || echo "MISSING: $key"
done < GRCm39.ref_map.tsv
```

## Step 4. Singularity 이미지 캐시

파이프라인에서 사용하는 모든 컨테이너 이미지를 로컬에 미리 다운로드합니다. 실행 중 네트워크 없이도 이미지를 사용할 수 있으며, 이미지 pull로 인한 실행 지연을 방지합니다.

```bash
bash scripts/populate_miniwdl_singularity_cache.sh \
  image_manifest.txt \
  miniwdl_singularity_cache/
```

## Step 5. 실행 확인

inputs.json을 준비하고 테스트 실행으로 환경을 검증합니다.

```bash
# inputs.json 작성 (BATCH_GUIDE.md 참조)
cp sample.inputs.json.example my_sample.inputs.json
# sample_id, hifi_reads, ref_map_file 등을 실제 값으로 수정
```

```bash
conda activate hifi-human-wgs
cd /data_4tb/hifi-human-wgs-wdl-custom

miniwdl run workflows/singleton.wdl \
  --input my_sample.inputs.json \
  --dir batch_results/ \
  --verbose
```

## Mouse 데이터 분석 시 주의사항

이 파이프라인은 human 전용으로 설계되었으나 mouse에서도 핵심 단계는 동작합니다:

| 도구 | Mouse 적용 | 비고 |
|------|-----------|------|
| pbmm2, DeepVariant, Sawfish, HiPhase | ✅ 정상 | Species-agnostic |
| TRGT | ✅ 부분적 | Mouse catalog 사용 (UCSC 기반) |
| pb-cpg-tools, MethBat | ✅ 정상 | CpG island TSV 필요 |
| PharmCAT, PBstarPhase | ❌ 건너뜀 | `run_pgx: false`로 비활성화 |
| Paraphase | ⚠️ 빈 결과 | Human HLA 전용 |
| Tertiary 분석 | ❌ 미적용 | Human population DB 필요 |

Mouse 실행 시 inputs.json에 반드시 포함:

```json
{
  "humanwgs_singleton.run_pgx": false
}
```

`tertiary_map_file`과 `phenotypes`는 지정하지 않습니다 (Mouse population DB 없음).
