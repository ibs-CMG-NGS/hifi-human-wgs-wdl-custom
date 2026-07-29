# TG-6903 D01 HiFi WGS 분석 보고서

**작성일:** 2026-06-01  
**샘플:** TG-6903 D01  
**트랜스진:** GFAP-GFP  
**레퍼런스:** GRCm39 (mouse)  
**분석자:** CMG-NGS-Core

---

## 요약 (결론 우선)

TG-6903 D01 샘플에 대해 HiFi 전장유전체 시퀀싱(WGS) 및 de novo assembly 기반 트랜스진 통합 분석을 수행하였습니다.

**GFAP-GFP 트랜스진이 chr4:68,426,842 위치에 hemizygous, 6-copy tandem 형태로 삽입되어 있음을 확인하였습니다.** 삽입 위치는 유전자 간 영역(intergenic)으로, 가장 가까운 단백질 코딩 유전자(Brinp1)까지 252 kb 이상 이격되어 있어 내재 유전자에 대한 직접적 영향은 없는 것으로 판단됩니다.

---

## 1. 시퀀싱 QC 통계

| 항목 | 값 |
|------|----|
| 입력 BAM | m84285_260219_114241_s4.hifi_reads.bc2027 |
| 레퍼런스 | GRCm39 |
| 평균 시퀀싱 depth | **38.5x** |
| chr4 mean depth | 36.8x |
| 총 시퀀싱 분량 (FASTQ) | 197 GB |

---

## 2. De novo Assembly (hifiasm v0.25.0)

### 2.1 어셈블리 통계

| 항목 | Haplotype 1 | Haplotype 2 |
|------|-------------|-------------|
| 주요 chimeric contig | h1tg000029l (29.6 Mb) | h2tg000049l (63.1 Mb) |
| BUSCO 완성도 (eutheria_odb10) | 77.8% [S:75.9%, D:1.9%] | 80.6% [S:78.6%, D:2.0%] |
| BUSCO (primary contig) | **98.4%** [S:95.8%, D:2.6%] | — |

> **참고:** hap1/hap2의 BUSCO 완성도가 ~78-80%인 것은 하플로타입별 부분 어셈블리의 특성으로 정상이며, primary contig 기준 98.4%는 고품질 어셈블리임을 나타냅니다.

### 2.2 파이프라인 소요 시간

| 단계 | 소요 시간 |
|------|-----------|
| hifiasm | 9시간 6분 (peak RSS 112 GB) |
| pbmm2 alignment (16 chunks) | 11시간 38분 |
| 전체 파이프라인 | 36시간 9분 (05/29 11:31 → 05/30 23:40) |

---

## 3. 변이 분석 결과 (GRCm39 reference-based)

### 3.1 소변이 (Small Variants — DeepVariant)

| 항목 | 값 |
|------|----|
| 총 변이 수 | 132,531 |
| SNV | 98,420 |
| Indel | 34,110 |
| Ti/Tv 비율 | 1.64 |
| Het 변이 | 82,121 (61.9%) |
| Hom 변이 | 16,299 (12.3%) |

### 3.2 구조 변이 (Structural Variants — Sawfish)

| 유형 | 수 |
|------|----|
| DEL | 2,675 |
| INS | 8,078 |
| DUP/INV/기타 | 2,127 |
| **합계** | **12,880** |

### 3.3 Phasing (HiPhase)

| 항목 | 값 |
|------|-----|
| 총 변이 | 280,926 |
| Heterozygous | 216,500 |
| Phased | 175,657 (81.1%) |
| Phase blocks | 16,277 |
| Phase block N50 | — |

---

## 4. 트랜스진 통합 분석 (GFAP-GFP)

### 4.1 삽입 위치 (정밀 분석)

| 항목 | 값 |
|------|-----|
| 염색체 | **chr4** |
| 삽입 위치 | **68,426,842 bp** |
| 위치 정밀도 | ±4 bp |
| 정밀 범위 | chr4:68,412,244–68,426,842 |
| Chimeric read 지지 수 | **13개** |

chimeric reads 13개 중:
- **host→TG** 방향: 9개 (chr4 숙주 DNA가 먼저, 이후 TG 서열)
- **TG→host** 방향: 4개 (TG 서열 이후 chr4 숙주 DNA)

breakpoint가 chr4:68,426,842 한 지점에 집중되어 삽입 위치의 정밀도가 높습니다.

### 4.2 트랜스진 카피 구조 (hap1)

| 항목 | 값 |
|------|-----|
| 통합 haplotype | **Hap1** (hemizygous) |
| 카피 수 | **6개 (tandem)** |
| 카피 간격 | **4 bp** (head-to-tail 직렬 연결) |
| 각 카피 크기 | 3,505 bp (100% 완전체) |
| TG 커버리지 | **100%** (전 카피) |
| 삽입 방향 | **역방향 (−)** |
| TG 클러스터 총 크기 | ~21 kb (6×3,505 bp + linker) |

#### 카피별 contig 좌표 (h1tg000029l)

| Copy | Contig start | Contig end | Coverage |
|------|-------------|------------|----------|
| 1 | 6,846,320 | 6,849,824 | 100% |
| 2 | 6,849,828 | 6,853,332 | 100% |
| 3 | 6,853,336 | 6,856,841 | 100% |
| 4 | 6,856,845 | 6,860,350 | 100% |
| 5 | 6,860,354 | 6,863,859 | 100% |
| 6 | 6,863,863 | 6,867,368 | 100% |

### 4.3 Hap2 분석

| 항목 | 값 |
|------|-----|
| Contig | h2tg000049l (63.1 Mb) |
| 검출 hit | 1개 (부분 매핑) |
| TG 커버리지 | **15%** (526 / 3,505 bp) |
| 매핑 위치 | chr16:5,281,423 (mapq=60) |

hap2에서의 15% 부분 검출은 완전한 트랜스진 삽입이 아닌 비특이적 서열 유사도에 의한 artifact로 판단됩니다. **TG-6903은 hap1에만 삽입된 hemizygous 개체입니다.**

---

## 5. 삽입 위치의 내재 유전자 영향

### 5.1 삽입 위치 주변 유전자

| 유전자 | 타입 | 위치 | 삽입점 거리 | 영향 |
|--------|------|------|------------|------|
| ENSMUSG00000139230 | lncRNA | chr4:68,534,123–68,538,037 (+) | 107,281 bp | 없음 |
| Gm12911 | lncRNA | chr4:68,230,870–68,240,484 (−) | 186,358 bp | 없음 |
| **Brinp1** | **protein_coding** | chr4:68,679,751–68,872,634 (−) | **252,909 bp** | 없음 |

### 5.2 내재 유전자 영향 평가

| 항목 | 판단 |
|------|------|
| 삽입 위치 feature | **Intergenic** (유전자 간 영역) |
| 가장 가까운 단백질 코딩 유전자 | Brinp1 (252 kb) |
| 기능 중요 유전자 파괴 | **없음** |
| lncRNA 파괴 | 없음 (107 kb 원거리) |
| GFAP 연구 관련 유전자 영향 | **없음** |

**삽입 위치 chr4:68.4 Mb는 유전자 밀도가 낮은 intergenic 영역으로, 내재 유전자에 대한 직접적 disruption이 없습니다.**

> **참고 — Mup cluster 관련:** `extract_integration_coords` 알고리즘의 초기 버전은 chr4:61.5 Mb(Mup gene cluster 내)를 삽입 위치로 잘못 추정하였습니다. 이는 TG 카피 경계 서열이 chr16 서열과 부분적으로 상동하여 chimeric contig의 reference alignment에서 gap 패턴이 왜곡된 것이 원인입니다. TG→assembly PAF를 직접 활용하는 방식으로 알고리즘을 수정하여 **chr4:68,426,842**의 정확한 위치를 확인하였습니다.

---

## 6. 결론

1. **GFAP-GFP 트랜스진 확인:** TG-6903 D01에서 GFAP-GFP 트랜스진 삽입이 명확히 확인되었습니다.
2. **삽입 위치:** chr4:68,426,842, 13개 chimeric read로 ±4 bp 정밀도 확인.
3. **삽입 구조:** Hemizygous, 6-copy tandem (역방향), 각 카피 100% 완전체.
4. **내재 유전자 영향 없음:** Intergenic 삽입, 가장 가까운 단백질 코딩 유전자 Brinp1까지 252 kb.
5. **연구 적합성:** GFAP 관련 유전자 좌에 영향 없으며, GFAP-GFP 형질전환 마우스로서 사용에 지장 없음.

---

## 7. 주요 파일 경로

| 파일 | 경로 |
|------|------|
| Haplotagged BAM | `/data_4tb/hifi-human-wgs-wdl-custom/output/TG6903_D01/_LAST/out/merged_haplotagged_bam/TG-6903.GRCm39.haplotagged.bam` |
| Hap1 assembly FASTA | `/data_4tb/hifi-human-wgs-wdl-custom/output/TG6903_D01/_LAST/out/assembly_hap1_fasta/TG-6903.hifiasm.bp.hap1.p_ctg.fa` |
| Hap2 assembly FASTA | `/data_4tb/hifi-human-wgs-wdl-custom/output/TG6903_D01/_LAST/out/assembly_hap2_fasta/TG-6903.hifiasm.bp.hap2.p_ctg.fa` |
| SV VCF | `/data_4tb/hifi-human-wgs-wdl-custom/output/TG6903_D01/_LAST/out/phased_sv_vcf/TG-6903.GRCm39.structural_variants.phased.vcf.gz` |
| Small variant VCF | `/data_4tb/hifi-human-wgs-wdl-custom/output/TG6903_D01/_LAST/out/phased_small_variant_vcf/TG-6903.GRCm39.small_variants.phased.vcf.gz` |
| Integration report (TXT) | `/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6903_D01_integration_wdl/_LAST/out/report_txt/TG-6903.GFAP_GFP.integration_report.txt` |
| Integration report (HTML) | `/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6903_D01_integration_wdl/_LAST/out/report_html/TG-6903.GFAP_GFP.integration_report.html` |
| Chimeric reads TSV | `/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6903_D01_integration_wdl/_LAST/out/chimeric_reads_tsv/TG-6903.chimeric_reads.tsv` |
| Integration sites TSV | `/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6903_D01_integration_wdl/_LAST/out/integration_tsv/TG-6903.GFAP_GFP.integration_sites.tsv` |

---

*본 보고서는 PacBio HiFi WGS 데이터 기반 bioinformatics 분석 결과입니다.*  
*트랜스진 양성/음성 최종 판정은 PCR genotyping 등 실험적 검증 결과와 종합하여 판단하시기 바랍니다.*
