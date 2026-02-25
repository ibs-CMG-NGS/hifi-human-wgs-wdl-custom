# PacBio WGS Variant Pipeline

PacBio HiFi whole-genome sequencing 데이터 분석을 위한 WDL 기반 파이프라인 (v3.1.1).

## 문서 안내

| 문서 | 읽어야 할 때 |
|------|------------|
| [SETUP.md](./SETUP.md) | 처음 환경 구축 시 |
| [BATCH_GUIDE.md](./BATCH_GUIDE.md) | 샘플 실행 시 (단일/배치) |
| [SERVER_GUIDE.md](./SERVER_GUIDE.md) | 서버 설정·GPU·리소스 변경 시 |
| [GIT_WORKFLOW.md](./GIT_WORKFLOW.md) | 코드/설정 버전 관리 시 |
| [TROUBLESHOOTING.md](./TROUBLESHOOTING.md) | 에러 발생 시 |

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
```

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
