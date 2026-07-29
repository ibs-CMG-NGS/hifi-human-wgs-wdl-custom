# TG-6425 C01 HiFi WGS 분석 보고서

**작성일:** 2026-05-29  
**샘플:** TG-6425 C01  
**트랜스진:** Aldh1l1-EGFP (예상)  
**레퍼런스:** GRCm39 (mouse)  
**분석자:** CMG-NGS-Core

---

## 요약 (결론 우선)

TG-6425 C01 샘플에 대해 HiFi 전장유전체 시퀀싱 분석을 수행하였습니다. 시퀀싱 품질과 어셈블리는 모두 정상적으로 완료되었으나, **Aldh1l1-EGFP 트랜스진 서열이 raw reads, de novo assembly, reference-based SV 분석 어느 단계에서도 검출되지 않았습니다.**

이는 분석 파이프라인의 문제가 아니며, **raw sequencing data 자체에 EGFP 서열이 존재하지 않는다는 것이 직접 확인되었습니다.**

> ⚠️ **연구팀 확인 필요:** TG-6425 C01 개체의 PCR genotyping 결과 재확인 및 샘플 동일성 검증을 권고합니다.

---

## 1. 시퀀싱 QC 통계

| 항목 | 값 |
|------|----|
| 총 reads 수 | 6,601,600 |
| 평균 read 길이 | 14,344 bp |
| Read N50 | 15,738 bp |
| 평균 시퀀싱 품질 | Q34 |
| 레퍼런스 매핑률 | 99.92% |
| 평균 시퀀싱 depth | 34.5x |
| 추정 성별 | FEMALE |

시퀀싱 품질은 정상 범위이며, 트랜스진 탐지에 충분한 coverage(34.5x)가 확보되었습니다.

---

## 2. De novo assembly 결과

| 항목 | Haplotype 1 | Haplotype 2 |
|------|-------------|-------------|
| 총 크기 | 2.38 Gb | 1.92 Gb |
| Contig 수 | 514 | 553 |
| BUSCO 완성도 | 98.3% | - |
| 어셈블리 품질 | 정상 | 정상 |

De novo assembly 자체는 정상적으로 수행되었습니다.

---

## 3. Aldh1l1-EGFP 트랜스진 탐지 분석

### 3.1 분석 방법

트랜스진 탐지는 다음 4가지 독립적인 방법으로 수행하였습니다.

| 방법 | 탐색 범위 | 탐색 대상 |
|------|-----------|-----------|
| ① Raw reads 전수 검색 | FASTQ 177G (660만 reads) | EGFP 특이적 20-mer |
| ② De novo assembly 검색 | hap1, hap2, unphased, raw unitig | EGFP + 트랜스진 FASTA |
| ③ Reference-based SV 분석 | 전체 게놈 | 3,000~5,500 bp novel insertion |
| ④ Haplotagged BAM 직접 검색 | chr6:86~91 Mb (Aldh1l1 locus) | EGFP 20-mer |

### 3.2 결과

| 방법 | 결과 |
|------|------|
| ① Raw reads (FASTQ) EGFP 20-mer 검색 (fwd) | **0 hits** |
| ① Raw reads (FASTQ) EGFP 20-mer 검색 (rev) | **0 hits** |
| ② hap1 / hap2 어셈블리 minimap2 정렬 | **alignment 없음** |
| ② p_ctg (unphased) / r_utg (raw unitig) | **alignment 없음** |
| ③ SV VCF — chr6 전체 novel INS ≥ 3 kb | **트랜스진 크기 삽입 없음** |
| ④ haplotagged BAM chr6:86~91 Mb EGFP 검색 | **0 hits** |

탐지 기준으로 사용한 EGFP 특이 서열:
- Forward: `ATGGTGAGCAAGGGCGAGGA` (EGFP 개시코돈 포함 20-mer)
- Reverse: `TCCTCGCCCTTGCTCACCAT`

### 3.3 Minimap2 PAF 분석 결과

트랜스진 FASTA(`Aldh1l1-EGFP_cassette`, 3,994 bp)를 de novo assembly에 정렬한 결과, hap1·hap2 모두에서 극히 작은 단편 alignment만 검출되었으며, 이는 트랜스진 내 Aldh1l1 5'UTR 부분(~200 bp)이 내재성(endogenous) Aldh1l1 유전자좌와 부분 상동성을 갖는 것으로 해석됩니다.

| 어셈블리 | 정렬 단편 | 위치 | 크기 | 비고 |
|----------|-----------|------|------|------|
| hap1 (h1tg000021l) | #1 | 20,022,127–20,022,413 | 286 bp | 임계값(500 bp) 미만 |
| hap1 (h1tg000021l) | #2 | 20,061,825–20,062,062 | 237 bp | 임계값(500 bp) 미만 |
| hap2 (h2tg000015l) | #1 | 20,022,127–20,022,413 | 286 bp | 동일 패턴 |
| hap2 (h2tg000015l) | #2 | 20,061,825–20,062,062 | 237 bp | 동일 패턴 |

두 단편 모두 reverse orientation, 서로 ~39 kb 이격 — 트랜스진 삽입 패턴이 아닌 숙주 유전체 내 부분 상동 서열로 판단됩니다.

---

## 4. 정상 양성 대조군과의 비교

동일 파이프라인으로 분석한 확인된 양성 개체(TG-6283 B01)와의 비교:

| 비교 항목 | TG-6283 B01 (양성 대조) | TG-6425 C01 (이번 분석) |
|-----------|------------------------|------------------------|
| 트랜스진 | gfa2-CreERT2 (4,309 bp) | Aldh1l1-EGFP (3,994 bp) |
| Raw reads 내 트랜스진 특이 서열 검출 | 다수 검출 (CreERT2 kmer) | **0** (EGFP 20-mer fwd/rev 모두) |
| De novo assembly chimeric contig | 검출됨 (h1tg000090l) | **검출 안 됨** |
| SV VCF novel INS | chr2:5,559,779 (5,315 bp) | **없음** |
| 2차 파이프라인 통과 여부 | 통과, 보고서 생성 | **실패 (삽입 위치 미검출)** |

---

## 5. 결론 및 권고사항

### 결론

1. TG-6425 C01의 시퀀싱 데이터 품질은 정상입니다.
2. **Aldh1l1-EGFP 트랜스진 서열이 raw sequencing data에 존재하지 않습니다.**
3. 이는 분석 파이프라인의 오류나 한계가 아닌, 데이터 자체의 특성입니다.
4. 2차 파이프라인(트랜스진 삽입 위치 탐지)은 이 결과로 인해 정상적으로 완료할 수 없습니다.

### 권고사항

| 우선순위 | 항목 |
|----------|------|
| ★★★ | TG-6425 C01 개체의 **PCR genotyping 원본 결과** 재확인 |
| ★★★ | 시퀀싱 라이브러리 준비 기록에서 **샘플 동일성(sample identity)** 확인 |
| ★★☆ | 가능하다면 **동일 개체로부터 재채혈 후 PCR genotyping** 재수행 |
| ★☆☆ | 트랜스진 FASTA(`Aldh1l1-EGFP_TG-6425.fa`)의 서열 정확성 검토 |

---

## 6. Misidentification 배제 분석 (Cross-check)

TG-6425 C01이 다른 TG line으로 mis-label되었을 가능성을 배제하기 위해, **TG-6425의 de novo assembly에 다른 3개 line의 transgene 서열을 매핑**하는 cross-check 분석을 수행하였습니다 (2026-06-04).

### 분석 방법

TG-6425의 HAP1/HAP2 assembly (hifiasm)를 대상으로, 아래 3개 transgene FASTA를 사용하여 동일한 TG integration detection 파이프라인을 실행하였습니다.

| 테스트 | Transgene FASTA | 대상 line |
|--------|----------------|----------|
| Cross-check #1 | `Aldh1l1_CreERT2_TG-6102.fa` | TG-6102 |
| Cross-check #2 | `gfa2_CreERT2_TG-6283.fa` | TG-6283 |
| Cross-check #3 | `GFAP_GFP_TG-6903.fa` | TG-6903 |

### 결과

| 테스트 | HAP1 chimeric contigs | HAP2 chimeric contigs | 판정 |
|--------|----------------------|----------------------|------|
| vs TG-6102 (Aldh1l1-CreERT2) | 0 | 0 | ✅ 음성 |
| vs TG-6283 (gfa2-CreERT2) | 0 | 0 | ✅ 음성 |
| vs TG-6903 (GFAP-GFP) | 1 (※) | 0 | ✅ 음성 (false positive) |

**※ GFAP-GFP 부분 hit 해석:**
- 검출된 contig: `h1tg000040l` (63 Mb 대형 배경 contig)
- TG coverage: **15%** (3505 bp 중 말단 526 bp, 위치 2979–3505만 매핑)
- 매핑 위치: chr16:5,281,198
- **False positive 근거:** 동일한 패턴(말단 15%, chr16:5,281 위치)이 TG-6903 자체 분석에서도 비특이적 hit로 확인됨. 이는 GFAP_GFP 말단부(SV40 polyA 추정)가 마우스 chr16 서열과 부분 유사성을 가지기 때문으로 판단되며, transgene 삽입과는 무관한 genome 배경 신호임.

### 결론

**TG-6425 C01이 TG-6102, TG-6283, TG-6903 중 어느 line으로도 mis-label되었을 가능성은 배제됩니다.**

본 샘플의 assembly에는 세 line의 transgene 서열(Aldh1l1 promoter, gfa2 promoter, GFAP-GFP)이 존재하지 않습니다. Aldh1l1-EGFP transgene 미검출 원인은 mis-labeling이 아닌, 해당 개체의 genotyping 결과 재확인이 필요한 상황입니다.

---

## 8. 주요 파일 경로

| 파일 | 경로 |
|------|------|
| HiFi reads (FASTQ) | `/data_4tb/.../call-bam_to_fastq-0/work/m84285_260219_093939_s3.hifi_reads.bc2026.fastq` |
| Haplotagged BAM | `/data_4tb/hifi-human-wgs-wdl-custom/output/TG6425_C01/_LAST/out/merged_haplotagged_bam/TG-6425.GRCm39.haplotagged.bam` |
| Hap1 assembly FASTA | `/data_4tb/.../out/assembly_hap1_fasta/TG-6425.hifiasm.bp.hap1.p_ctg.fa` |
| Hap2 assembly FASTA | `/data_4tb/.../out/assembly_hap2_fasta/TG-6425.hifiasm.bp.hap2.p_ctg.fa` |
| SV VCF | `/data_4tb/.../out/phased_sv_vcf/TG-6425.GRCm39.structural_variants.phased.vcf.gz` |
| 분석 통계 | `/data_4tb/.../out/stats_file/TG-6425.stats.txt` |

---

*본 보고서는 PacBio HiFi WGS 데이터 기반 bioinformatics 분석 결과입니다.  
트랜스진 양성/음성 최종 판정은 PCR genotyping 등 실험적 검증 결과와 종합하여 판단하시기 바랍니다.*
