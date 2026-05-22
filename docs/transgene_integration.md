# Transgene Integration Site Analysis

트랜스진(transgene)이 게놈에 삽입된 위치를 de novo assembly 기반으로 분석하는 파이프라인입니다.

## 워크플로우 개요

```
workflows/transgene_integration.wdl
```

### 분석 단계

```
[De novo 기반 — 항상 실행]
  map_transgene_to_assembly   트랜스진 → hap1/hap2 contig 매핑 (minimap2 asm20)
    → extract_chimeric_contigs  트랜스진 서열을 포함하는 chimeric contig 추출
      → align_chimeric_to_ref   chimeric contig → 레퍼런스 정렬 (minimap2 asm5)
        → extract_integration_coords  대략적 삽입 위치(chr/pos) 추출

[Hybrid Reference 정밀 분석 — haplotagged_bam 제공 시]
  build_hybrid_ref            삽입 염색체 + 트랜스진 FASTA로 hybrid reference 구축
    → extract_region_reads    삽입 위치 ±window bp 범위 HiFi reads 추출
      → align_to_hybrid_ref   reads → hybrid reference 재정렬 (minimap2 map-hifi)
        → detect_chimeric_reads  chimeric reads 검출 및 정확한 breakpoint 결정

[유전자 어노테이션 — annotation_gtf_bgz 제공 시]
  query_gene_annotation       breakpoint 주변 유전자 조회 (tabix + Gencode GTF)

[리포트 — 항상 실행]
  integration_report          txt / TSV / HTML 통합 리포트 생성
```

### 핵심 개념

- **Chimeric contig**: 트랜스진 서열과 게놈 서열이 하나의 contig에 함께 조립된 것. 삽입이 일어난 위치의 직접적 증거
- **Hybrid reference 분석**: 삽입 염색체 + 트랜스진을 이어붙인 레퍼런스에 HiFi reads를 재정렬 → split alignment로 정확한 breakpoint 결정
- **Breakpoint**: chimeric reads의 host-TG junction 좌표 중앙값. 수십 개 reads로 지지될 때 수 bp 수준의 정밀도

## 입력 파라미터

### 필수

| 파라미터 | 설명 |
|---|---|
| `sample_id` | 샘플 ID (예: `TG-6283`) |
| `transgene_name` | 트랜스진 구축명 (예: `gfa2_CreERT2`) |
| `transgene_fasta` | 트랜스진 서열 FASTA |
| `hap1_fasta` | Haplotype 1 assembly FASTA (hifiasm 출력) |
| `hap2_fasta` | Haplotype 2 assembly FASTA (hifiasm 출력) |
| `ref_fasta` | 레퍼런스 게놈 FASTA (예: GRCm39) |
| `default_runtime_attributes` | 런타임 설정 |

### 선택 — Hybrid Reference 정밀 분석

| 파라미터 | 기본값 | 설명 |
|---|---|---|
| `haplotagged_bam` | — | HiPhase 해플로태깅 BAM (humanwgs 파이프라인 출력) |
| `haplotagged_bai` | — | BAM 인덱스 (.bai) |
| `extract_window_bp` | 300000 | 삽입 위치 기준 read 추출 범위 (bp) |

### 선택 — 유전자 어노테이션

| 파라미터 | 기본값 | 설명 |
|---|---|---|
| `annotation_gtf_bgz` | — | Sorted + bgzipped Gencode GTF |
| `annotation_gtf_tbi` | — | Tabix 인덱스 (.tbi) |
| `annotation_window_bp` | 200000 | Breakpoint 기준 유전자 조회 범위 (bp) |

### 기타

| 파라미터 | 기본값 | 설명 |
|---|---|---|
| `min_match_bp` | 500 | Chimeric contig 판정 최소 매칭 bp |

## 출력 파일

| 파일 | 설명 |
|---|---|
| `report_html` | 통합 리포트 (HTML, 주요 결과) |
| `report_txt` | 통합 리포트 (텍스트) |
| `integration_tsv` | 삽입 위치 상세 테이블 |
| `hap1_chimeric_fasta` / `hap2_chimeric_fasta` | Chimeric contig 서열 |
| `chimeric_reads_tsv` | Chimeric reads 상세 (hybrid ref 분석 시) |
| `annotation_tsv` | 삽입 위치 주변 유전자 목록 (어노테이션 분석 시) |

## 실행 방법

### inputs.json 작성

`{SAMPLE}.tg_integration.inputs.json` 형식으로 작성합니다.

**기본 (de novo 분석만)**

```json
{
  "transgene_integration.sample_id":       "TG-XXXX",
  "transgene_integration.transgene_name":  "transgene_construct_name",
  "transgene_integration.transgene_fasta": "/path/to/transgene.fa",
  "transgene_integration.hap1_fasta":      "/path/to/hap1.p_ctg.fa",
  "transgene_integration.hap2_fasta":      "/path/to/hap2.p_ctg.fa",
  "transgene_integration.ref_fasta":       "/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/mouse_GRCm39.fasta",
  "transgene_integration.default_runtime_attributes": {
    "backend": "HPC", "preemptible_tries": 0, "max_retries": 2,
    "zones": "", "cpuPlatform": "", "gpuType": "",
    "container_registry": "quay.io/biocontainers", "container_namespace": null
  }
}
```

**전체 분석 (hybrid ref + 어노테이션 포함)**

위에 아래 항목을 추가합니다.

```json
{
  "transgene_integration.haplotagged_bam": "/path/to/TG-XXXX.GRCm39.haplotagged.bam",
  "transgene_integration.haplotagged_bai": "/path/to/TG-XXXX.GRCm39.haplotagged.bam.bai",
  "transgene_integration.annotation_gtf_bgz": "/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/gencode.vM36.annotation.sorted.gtf.bgz",
  "transgene_integration.annotation_gtf_tbi": "/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/gencode.vM36.annotation.sorted.gtf.bgz.tbi"
}
```

### 실행

```bash
conda activate hifi-human-wgs
cd /data_4tb/hifi-human-wgs-wdl-custom

bash run_tg_integration.sh {SAMPLE_ID} {SAMPLE}.tg_integration.inputs.json \
  /mnt/JJ_dis_8tb/tg-integration-denovo-assembly/{SAMPLE}_integration_analysis
```

또는 직접 실행:

```bash
miniwdl run workflows/transgene_integration.wdl \
  --input {SAMPLE}.tg_integration.inputs.json \
  --dir /mnt/JJ_dis_8tb/tg-integration-denovo-assembly/{SAMPLE}_integration_analysis \
  --verbose 2>&1 | tee {SAMPLE}.tg_integration.log
```

## 레퍼런스 데이터 (공용)

| 파일 | 경로 |
|---|---|
| GRCm39 레퍼런스 | `/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/mouse_GRCm39.fasta` |
| Gencode vM36 GTF (bgzip) | `/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/gencode.vM36.annotation.sorted.gtf.bgz` |
| Gencode vM36 GTF index | `/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/gencode.vM36.annotation.sorted.gtf.bgz.tbi` |
| 트랜스진 FASTA 디렉토리 | `/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG FA files/` |

## 샘플별 입력 파일

| 샘플 | inputs.json | 트랜스진 FASTA |
|---|---|---|
| TG-6283 B01 | `TG6283.tg_integration.inputs.json` | `gfa2_CreERT2_TG-6283.fa` |
| TG-6102 A01 | `TG6102.tg_integration.inputs.json` | `Aldh1l1_CreERT2_TG-6102.fa` |

## 분석 결과 경로

```
/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/
  {SAMPLE}_integration_analysis/_LAST/out/
    report_html/          ← HTML 리포트
    report_txt/           ← 텍스트 리포트
    integration_tsv/      ← 삽입 위치 테이블
    chimeric_reads_tsv/   ← Chimeric reads 상세
    annotation_tsv/       ← 유전자 어노테이션
```

## 선행 조건

이 파이프라인은 humanwgs_singleton 파이프라인의 출력을 입력으로 사용합니다.

```
humanwgs_singleton (run_assembly=true)
  ├── hap1_fasta  →  transgene_integration.hap1_fasta
  ├── hap2_fasta  →  transgene_integration.hap2_fasta
  └── haplotagged_bam  →  transgene_integration.haplotagged_bam (선택)
```

haplotagged BAM 경로는 다음 위치에서 찾을 수 있습니다.

```
{SAMPLE_OUTPUT_DIR}/_LAST/out/merged_haplotagged_bam/TG-XXXX.GRCm39.haplotagged.bam
```
