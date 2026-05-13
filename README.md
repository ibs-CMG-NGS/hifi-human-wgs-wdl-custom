# PacBio WGS Variant Pipeline

PacBio HiFi whole-genome sequencing 데이터 분석을 위한 WDL 기반 파이프라인 (v3.1.1).

## 문서 안내

### 운영 가이드

이 서버 환경에 맞춘 실무 가이드입니다.

| 문서 | 읽어야 할 때 |
|------|------------|
| [INSTALL.md](./INSTALL.md) | 처음 서버에 필수 도구(Conda, Apptainer, NVIDIA 드라이버 등) 설치 시 |
| [SETUP.md](./SETUP.md) | 파이프라인 환경 구축 시 (miniwdl 설정, 레퍼런스 데이터, 이미지 캐시) |
| [BATCH_GUIDE.md](./BATCH_GUIDE.md) | 샘플 실행 시 — inputs.json 작성, 단일/배치 실행, 결과 구조, QC 리포트 |
| [SERVER_GUIDE.md](./SERVER_GUIDE.md) | 서버 스펙·miniwdl 설정·GPU 설정·디스크 관리 |
| [TROUBLESHOOTING.md](./TROUBLESHOOTING.md) | 에러 발생 시 — Cross-device link, Singularity, GPU, Sawfish, Mouse 분석 오류 |
| [GIT_WORKFLOW.md](./GIT_WORKFLOW.md) | 코드/설정 버전 관리, Windows ↔ 서버 동기화 |

### 기술 레퍼런스

워크플로우 입출력 스펙, 파일 포맷 명세, 도구별 상세 문서입니다.

**워크플로우 I/O**

| 문서 | 내용 |
|------|------|
| [docs/singleton.md](./docs/singleton.md) | singleton.wdl 전체 입력/출력 파라미터 및 DAG |
| [docs/family.md](./docs/family.md) | family.wdl 전체 입력/출력 파라미터 및 DAG |

**파일 포맷 명세**

| 문서 | 내용 |
|------|------|
| [docs/ref_map.md](./docs/ref_map.md) | ref_map TSV 키 목록 및 설명 |
| [docs/tertiary_map.md](./docs/tertiary_map.md) | tertiary_map TSV 키 목록 및 설명 |
| [docs/tertiary.md](./docs/tertiary.md) | Tertiary 분석 (slivar/svpack) 상세 |

**도구별 상세**

| 문서 | 내용 |
|------|------|
| [docs/assembly.md](./docs/assembly.md) | De novo assembly (hifiasm + BUSCO + Merqury + SyRI + QUAST) |
| [docs/deepvariant.md](./docs/deepvariant.md) | DeepVariant 설정 |
| [docs/pbmm2.md](./docs/pbmm2.md) | pbmm2 정렬 옵션 |
| [docs/trgt.md](./docs/trgt.md) | TRGT 탠덤반복 genotyping |
| [docs/pharmcat.md](./docs/pharmcat.md) | PharmCAT 약물유전체 분석 |
| [docs/gpu.md](./docs/gpu.md) | GPU 지원 설정 및 클라우드별 GPU 타입 |
| [docs/bam_statistics.md](./docs/bam_statistics.md) | BAM 통계 항목 설명 |
| [docs/tools_containers.md](./docs/tools_containers.md) | 도구 버전 및 컨테이너 이미지 목록 |

**백엔드 설정**

| 문서 | 내용 |
|------|------|
| [docs/backend-hpc.md](./docs/backend-hpc.md) | HPC (이 서버) 설정 |
| [docs/backend-gcp.md](./docs/backend-gcp.md) | Google Cloud 설정 |
| [docs/backend-azure.md](./docs/backend-azure.md) | Azure 설정 |
| [docs/backend-dnanexus.md](./docs/backend-dnanexus.md) | DNAnexus 설정 |
| [docs/backend-aws-healthomics.md](./docs/backend-aws-healthomics.md) | AWS HealthOmics 설정 |

## 워크플로우 구성

두 가지 진입점:

- `workflows/singleton.wdl` — 단일 샘플 분석
- `workflows/family.wdl` — 가족/다중 샘플 joint calling

### 분석 단계 (Singleton 기준)

```
pbmm2 (정렬)
  → mosdepth (커버리지)
  → DeepVariant (소변이 SNV/Indel)
  → Sawfish (구조변이 SV)
  → TRGT (탠덤반복)
  → Paraphase (HLA) / MitorSaw (mtDNA)
  → HiPhase (위상결정·해플로태깅)
  → pb-cpg-tools + MethBat (메틸화)
  → PBstarPhase + PharmCAT (약물유전체, human only)
  → [선택] Tertiary 분석 (slivar + svpack)

[선택, run_assembly=true]
  → hifiasm (de novo assembly)
  → BUSCO + Merqury (assembly QC)
  → minimap2 + SyRI (assembly-to-reference 구조변이)
  → QUAST (reference 대비 assembly 품질)
```

> De novo assembly는 `run_assembly=true` 플래그로 활성화. 기본값 `false`. → [상세 가이드](./docs/assembly.md)

## 레퍼런스 데이터

| 종 | 경로 | 레퍼런스 맵 |
|----|------|------------|
| Human (GRCh38) | `hifi-wdl-resources/hifi-wdl-resources-v3.1.0/GRCh38/` | `backends/hpc/GRCh38.ref_map.v3p1p0.hpc.tsv` |
| Mouse (GRCm39) | `hifi-wdl-resources/GRCm39/` | `GRCm39.ref_map.tsv` |

## 주요 경로 (이 서버)

```
/data_4tb/hifi-human-wgs-wdl-custom/           ← 파이프라인 루트
/data_4tb/pacbio_rawdata/                       ← Raw HiFi BAM
/data_4tb/hifi-human-wgs-wdl-custom/batch_results/  ← 분석 결과
```

## 빠른 실행

```bash
conda activate hifi-human-wgs
cd /data_4tb/hifi-human-wgs-wdl-custom

miniwdl run workflows/singleton.wdl \
  --input BioSample24.inputs.json \
  --dir batch_results/ \
  --verbose 2>&1 | tee batch_results/BioSample24.run.log
```

---

*Based on [PacBio HiFi-human-WGS-WDL v3.1.1](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL)*
