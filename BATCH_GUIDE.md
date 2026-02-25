# 샘플 실행 가이드

## 빠른 시작

```bash
conda activate hifi-human-wgs
cd /data_4tb/hifi-human-wgs-wdl-custom

miniwdl run workflows/singleton.wdl \
  --input BioSample24.inputs.json \
  --dir batch_results/ \
  --verbose 2>&1 | tee batch_results/BioSample24.run.log
```

## inputs.json 작성

### 기본 구조

```json
{
  "humanwgs_singleton.sample_id":    "SAMPLE_ID",
  "humanwgs_singleton.sex":          "MALE",
  "humanwgs_singleton.hifi_reads":   ["/path/to/sample.hifi_reads.bam"],
  "humanwgs_singleton.ref_map_file": "/path/to/ref_map.tsv",
  "humanwgs_singleton.backend":      "HPC",
  "humanwgs_singleton.preemptible":  false,
  "humanwgs_singleton.gpu":          false
}
```

### 주요 파라미터

| 파라미터 | 필수 | 설명 |
|----------|------|------|
| `sample_id` | ✅ | 영문자·숫자·`.`·`-`·`_` 허용 |
| `sex` | 선택 | `"MALE"` / `"FEMALE"` (미입력 시 XX 기본값) |
| `hifi_reads` | ✅ | HiFi unaligned BAM 경로 배열 |
| `fail_reads` | 선택 | failed reads BAM (bait capture에 사용) |
| `ref_map_file` | ✅ | 레퍼런스 TSV 경로 |
| `tertiary_map_file` | 선택 | Tertiary 분석 TSV (미입력 시 스킵) |
| `phenotypes` | 선택 | HPO 코드 (미입력 시 tertiary 스킵) |
| `gpu` | 선택 | DeepVariant GPU 가속 (기본값: false) |
| `max_reads_per_alignment_chunk` | 선택 | 기본값: 500000 |

### 레퍼런스 맵 경로

| 종 | ref_map_file |
|----|-------------|
| Human (GRCh38) | `backends/hpc/GRCh38.ref_map.v3p1p0.hpc.tsv` |
| Mouse (GRCm39) | `GRCm39.ref_map.tsv` |

### BAM 파일 위치 패턴 (r84285_20260219 런)

```
/data_4tb/pacbio_rawdata/r84285_20260219_052427/
├── 1_A01/hifi_reads/m84285_260219_053344_s1.hifi_reads.bc2024.bam  → BioSample24
├── 1_B01/hifi_reads/m84285_260219_073640_s2.hifi_reads.bc2025.bam  → BioSample25
├── 1_C01/hifi_reads/m84285_260219_093939_s3.hifi_reads.bc2026.bam  → BioSample26
└── 1_D01/hifi_reads/m84285_260219_114241_s4.hifi_reads.bc2027.bam  → BioSample27
```

> `unassigned.bam` 파일은 사용하지 않음 (바코드 미매칭 reads)

### 새 런 BAM 파일 확인 방법

```bash
# hifi_reads BAM 목록
find /data_4tb/pacbio_rawdata/<RUN_DIR>/ -path "*/hifi_reads/*.bam" | grep -v unassigned | sort

# 샘플 ID 확인 (BAM 헤더)
samtools view -H /path/to/sample.bam | grep "^@RG" | grep -oP "SM:\S+"
```

### inputs.json 일괄 생성 스크립트

```bash
cat > /tmp/make_inputs.py << 'PYEOF'
import json

RAW = "/data_4tb/pacbio_rawdata/RUN_DIR"       # 실제 런 디렉토리로 변경
REF = "/data_4tb/hifi-human-wgs-wdl-custom/GRCm39.ref_map.tsv"  # 또는 GRCh38

samples = [
    {"id": "SAMPLE_ID_1", "dir": "1_A01", "bam": "BAMFILE1.bam"},
    {"id": "SAMPLE_ID_2", "dir": "1_B01", "bam": "BAMFILE2.bam"},
]

for s in samples:
    d = {
        "humanwgs_singleton.sample_id":    s["id"],
        "humanwgs_singleton.hifi_reads":   [f"{RAW}/{s['dir']}/hifi_reads/{s['bam']}"],
        "humanwgs_singleton.ref_map_file": REF,
        "humanwgs_singleton.backend":      "HPC",
        "humanwgs_singleton.preemptible":  False,
        "humanwgs_singleton.gpu":          False,
    }
    fname = f"/data_4tb/hifi-human-wgs-wdl-custom/{s['id']}.inputs.json"
    with open(fname, "w") as f:
        json.dump(d, f, indent=2)
    print(f"생성: {fname}")
PYEOF

python3 /tmp/make_inputs.py
```

## 실행

### 단일 샘플

```bash
cd /data_4tb/hifi-human-wgs-wdl-custom

miniwdl run workflows/singleton.wdl \
  --input BioSample24.inputs.json \
  --dir batch_results/ \
  --verbose \
  2>&1 | tee batch_results/BioSample24.run.log
```

### 다중 샘플 순차 실행

```bash
cd /data_4tb/hifi-human-wgs-wdl-custom

for sample in BioSample24 BioSample25 BioSample26 BioSample27; do
  echo "=== 시작: $sample $(date) ==="
  miniwdl run workflows/singleton.wdl \
    --input ${sample}.inputs.json \
    --dir batch_results/ \
    --verbose \
    2>&1 | tee batch_results/${sample}.run.log
  echo "=== 완료: $sample $(date) ==="
done
```

> `task_concurrency = 1` 설정이므로 동시 실행보다 순차 실행이 효율적.

### 백그라운드 실행

```bash
nohup bash -c '
for sample in BioSample24 BioSample25 BioSample26 BioSample27; do
  miniwdl run workflows/singleton.wdl \
    --input ${sample}.inputs.json \
    --dir batch_results/ \
    --verbose \
    2>&1 | tee batch_results/${sample}.run.log
done' > batch_results/batch_master.log 2>&1 &

echo "PID: $!"
```

## 진행 상황 모니터링

별도 터미널에서:

```bash
# monitor_progress.sh 사용
./monitor_progress.sh batch_results/BioSample24.run.log

# 간단 확인
watch -n 10 'grep -c "task done" batch_results/BioSample24.run.log'
```

예상 태스크 수: ~18단계 (singleton, tertiary 제외 시)

## 결과 구조

```
batch_results/
└── _YYYY-MM-DD_HH-MM-SS_humanwgs_singleton/
    └── out/
        ├── merged_haplotagged_bam/        ← 정렬 BAM (haplotagged)
        ├── phased_small_variant_vcf/      ← 소변이 VCF (SNV/Indel)
        ├── small_variant_gvcf/            ← gVCF (joint calling용)
        ├── phased_sv_vcf/                 ← 구조변이 VCF
        ├── phased_trgt_vcf/               ← 탠덤반복 VCF
        ├── cpg_combined_bed/              ← 메틸화 BED
        ├── cpg_hap1_bed/ cpg_hap2_bed/   ← 해플로타입별 메틸화
        ├── pbstarphase_json/              ← PGx 다이플로타입
        ├── pharmcat_report_html/          ← 약물유전체 보고서
        ├── mosdepth_summary/              ← 커버리지 요약
        ├── phase_stats/                   ← 위상결정 통계
        └── stats_file/                    ← 전체 통계 TSV
```

## 다중 BAM 파일 처리

한 샘플에 여러 BAM 파일이 있을 때 (여러 cell 또는 run):

```json
{
  "humanwgs_singleton.hifi_reads": [
    "/path/to/sample_cell1.bam",
    "/path/to/sample_cell2.bam",
    "/path/to/sample_cell3.bam"
  ]
}
```

파이프라인이 내부에서 자동으로 merge함.

## Call Cache 활용

재실행 시 성공한 태스크는 자동으로 캐시에서 불러옴:

```
miniwdl_call_cache/   ← 태스크별 결과 캐시
```

특정 태스크만 재실행하려면 해당 캐시 디렉토리 삭제 후 재실행.
