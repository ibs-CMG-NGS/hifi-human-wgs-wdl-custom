# Base Modification Discovery Pipeline — 설계 문서 (계획 v3)

> 이 문서는 파이프라인 구축 전 작성된 **설계/계획 문서**입니다. 문헌 근거와 설계 의도를
> 보존하기 위해 원문 그대로 남겨두었습니다. **실제 입력/출력 파라미터, 실행 방법은
> [docs/base_modification_discovery.md](./base_modification_discovery.md) 를 기준으로 확인하세요.**
> 구현 과정에서 일부 도구/출력이 계획과 달라졌으며, 차이점은 문서 맨 아래
> [구현 시 변경된 사항](#구현-시-변경된-사항) 절에 정리했습니다.

## Context

PacBio HiFi kinetics(IPD/PW) 데이터를 이용하여 oxidative stress 약물이 마우스 게놈에 유발하는
**알려지지 않은 DNA 염기 변형(unknown base modification)** 을 genome-wide로 탐색하고,
알려진 변형(5mC, **5hmC**, 6mA)과 통합하여 포괄적인 에피게놈 변화 지도를 구축하는 WDL 파이프라인.

### 실험 설계

| 그룹 | 역할 | 예상 변형 신호 |
|------|------|---------------|
| **Control** | 기준선 | baseline |
| **Drug A** | ROS 유도제 | 8-oxoG↑, FapyG/FapyA↑, 5hmC↑(TET 활성), 5mC→5hmC 전환↑ |
| **Drug B** | 항산화제 | baseline 또는 소폭 ↓ |
| **Drug A+B** | 복합 처치 | Drug A 효과 부분 구제(rescue) 예상 |

- 각 그룹 3마리 → pooling 후 ~90x 유효 커버리지
- 종: 마우스 (GRCm39)
- 데이터 위치: `/mnt/` 하위 외장 드라이브

### 핵심 비교 축

```
① Drug A vs Control          → ROS 유도 변형 (induction)
② Drug A+B vs Drug A         → Rescue 변형 ← 가장 중요
③ Drug B vs Control          → 항산화 단독 효과
④ Drug A+B vs Control        → 복합 처치 순효과
```

### 변형 종류별 탐지 전략

| 변형 | 탐지 방법 | 도구 |
|------|-----------|------|
| **5mC** | MM/ML 태그 + ML 모델 | jasmine → pb-CpG-tools |
| **5hmC** | ip/pw 키네틱스 + ML 모델 | jasmine (5mC와 구분) |
| **6mA** | ip/pw 키네틱스 | jasmine 또는 ipdSummary |
| **Unknown** | ip/pw IPD ratio (vs control) | ipdSummary (two-sample) |

---

## 기존 리소스 재사용

| 리소스 | 경로 | 재사용 Task |
|--------|------|-------------|
| RuntimeAttributes struct | `workflows/wdl-common/wdl/structs.wdl` | 전체 |
| GRCm39 FASTA | `hifi-wdl-resources/GRCm39/mouse_GRCm39.fasta` | 1,3,6 |
| Gencode vM36 GTF | `hifi-wdl-resources/GRCm39/gencode.vM36.annotation.sorted.gtf.bgz` | 8,10 |
| CpG island TSV | `hifi-wdl-resources/GRCm39/cpgIslandExt.sorted.mm39.tsv` | 8 |
| `cpg_pileup` task | `workflows/wdl-common/wdl/tasks/cpg_pileup.wdl` | Task 3c |

---

## 신규 파일 목록

| 파일 | 설명 |
|------|------|
| `workflows/base_modification_discovery.wdl` | 메인 워크플로우 |
| `workflows/wdl-common/wdl/tasks/base_modification.wdl` | 신규 Task 정의 |
| `base_mod_discovery.inputs.template.json` | 입력 템플릿 |
| `run_base_mod_discovery.sh` | 실행 스크립트 |

---

## 전체 워크플로우 구조

```
12개 raw HiFi BAM (/mnt/, kinetics ip/pw + MM/ML 태그 포함)
    │
    ▼ [scatter: 12 병렬]
Task 1. pbmm2_align_kinetics
  → 12개 aligned BAM (ip/pw + MM/ML 모두 보존)
    │
    ▼ [scatter: 12 병렬 — 정렬 직후 per-sample ML 실행]
Task 2a. jasmine_modification_calling
  → 12개 BAM (MM/ML 업데이트: 5mC + 5hmC + 6mA 확률 추가)
    │
    ▼ [4 그룹별 병렬]
Task 2b. merge_group_bams
  → 4개 pooled BAM (~90x, jasmine 태그 포함)
    │
    ├──────────────────────────────────────────────────────┐
    │                                                      │
    ▼ [3 비교 병렬]                    ▼ [4 그룹 병렬]   │
Task 3a. ipd_summary               Task 3b. jasmine_pileup │
  (treated vs control, IPD 기반)     (5mC / 5hmC / 6mA    │
  → unknown 포함 전체 변형 탐지       per-site 확률 BED)   │
  → modifications.gff/.csv          + Task 3c. cpg_pileup  │
                                       (5mC 기존 task 비교용)
    │                                          │
    ▼ [3 병렬]                               ▼
Task 4a. filter_and_classify         Task 4b. compare_known_mods
  → unknown_candidates.tsv            → 5mC/5hmC/6mA 그룹 비교
  → known_mod_summary.tsv             → differential_mods.tsv

    │                                          │
    └──────────────┬───────────────────────────┘
                   ▼
Task 5. multigroup_comparison
  → categorized_modifications.tsv
  → rescue_candidates.tsv (A→A+B에서 복구된 위치)
  → 5hmC_rescue.tsv (Drug A로 증가한 5hmC가 Drug B로 복구되는지)
  → Venn 요약
    │
    ├──────────┬──────────┬──────────┬──────────┐
    ▼          ▼          ▼          ▼          ▼
Task 6.   Task 7.    Task 8.    Task 9.    Task 10.
bigwig    motif      annotate   pathway    report
(IGV용)   (MEME)     (GTF+CpGi) (GO/KEGG)  (HTML)
```

---

## Task 상세 설계

### Task 1: `pbmm2_align_kinetics`
```
목적: ip/pw(kinetics) + MM/ML(base mod) 태그 모두 보존하며 정렬
입력: raw_hifi_bam, ref_fasta, sample_id
출력: aligned.bam, aligned.bam.bai
컨테이너: quay.io/pacbio/pbmm2:1.13.0 (기존 파이프라인과 동일 버전)
핵심: strip_kinetics=false (기존 pbmm2.wdl 기본값 override)
스레드: 16, 메모리: 48 GiB
```
> **문헌**: Flusberg et al. (2010) *Nature Methods* 7:461-465 — SMRT sequencing kinetics로
> DNA 변형 직접 탐지 최초 증명. 기존 싱글톤 파이프라인은 `pbmm2.wdl` L.269에서
> `strip_kinetics=true`로 이 정보를 소실시킴.

### Task 2a: `jasmine_modification_calling`
```
목적: kinetics(ip/pw)를 ML 모델에 적용 → 5mC / 5hmC / 6mA 확률을 MM/ML 태그로 기록
입력: aligned.bam (kinetics 보존)
출력: modified.bam (MM/ML 태그 업데이트 — 5mC, 5hmC, 6mA 구분)
컨테이너: quay.io/pacbio/jasmine:2.0.0
명령:
  jasmine kinetics \
    --bam aligned.bam \
    --output-bam modified.bam \
    --log-level INFO
스레드: 8, 메모리: 32 GiB
```
> **5hmC 탐지 근거**:
> - PacBio (2023) Jasmine Technical Note — kinetics-based ML 모델로 5mC와 5hmC의 서로
>   다른 IPD/PW 패턴을 구분. 5hmC는 5mC보다 약한 IPD 효과이나 SMRT kinetics로 식별 가능.
> - Eid et al. (2009) *Science* 323:133 — SMRT 방법론 기반 논문: 각 염기의 통과 시간이
>   화학적 구조에 고유하게 의존함.
> - Laszlo et al. (2013) *Nature Biotechnology* — 단일분자 수준에서 5mC vs 5hmC 구분 가능성 보고.
> - **5hmC in oxidative stress**: Ito et al. (2011) *Science* 333:1300 — TET 효소가 ROS
>   환경에서 5mC→5hmC 전환을 가속. Drug A(ROS 유도) 처치 시 5hmC 증가 예상.

### Task 2b: `merge_group_bams`
```
목적: 동일 그룹 3 replicate 병합 → ~90x 유효 커버리지 확보
입력: Array[File] modified_bams (jasmine 처리된 3개), out_prefix
출력: merged.bam, merged.bam.bai
컨테이너: quay.io/biocontainers/samtools:1.21--h50ea8bc_0
명령: samtools merge → samtools sort -@ 8 → samtools index
```
> **문헌**: Schadt et al. (2013) *Genome Research* 23:129-141 — 커버리지가 kinetics 기반
> 변형 탐지 통계력의 핵심 결정인자. 25-35x 개별 샘플 pooling → ~90x로 false positive rate 감소.

### Task 3a: `ipd_summary`
```
목적: per-base IPD ratio로 unknown modification 포함 모든 변형 위치 탐지
입력: treated_bam (pooled), control_bam (pooled), ref_fasta, out_prefix
출력: modifications.gff, modifications.csv
컨테이너: quay.io/biocontainers/kineticstools:2.5.0--py311h4f5f2e0_3
명령:
  ipdSummary treated.bam \
    --reference ref.fasta \
    --control control.bam \
    --gff modifications.gff \
    --csv modifications.csv \
    --pvalue 0.001 \
    --minCoverage 25 \
    --numWorkers 16
스레드: 16, 메모리: 48 GiB
```
> **예상 unknown modification 후보 및 문헌**:
> - **8-oxoG**: Cadet & Wagner (2013) *Cold Spring Harb Perspect Biol* 5:a012559 — 산화 스트레스의
>   가장 빈번한 DNA 손상. GGG/G4 서열에 집중. OGG1 BER 기질.
>   Lujan et al. (2012) *Nucleic Acids Res* — PacBio로 8-oxoG 간접 탐지 가능성 보고.
> - **FapyG / FapyA**: 동일 Cadet 리뷰 — 퓨린 고리 개열로 중합효소 속도 현저 저하 → 높은 IPD.
> - **Thymine glycol**: Cooke et al. (2003) *FASEB J* 17:1195 — 피리미딘 산화 손상.
> - **5fC / 5caC**: 5hmC의 추가 산화 산물. Ito et al. (2011) *Science* 333:1300.

### [hifimeth 검토 결과]
> - **hifimeth** (GitHub: xiaochuanle/hifimeth, bioRxiv 2024.08.14.607879): deep graph CNN 기반 5mC 단일분자 calling. 5mC 정확도 ~95%, AUC 98% 이상으로 jasmine보다 5mC 단독 성능 우수.
> - **그러나 5hmC 미지원** — jasmine만이 5mC/5hmC/6mA 동시 탐지 가능.
> - conda/bioconda 패키지 없음 → 별도 컨테이너 빌드 필요.
> - **결론**: jasmine을 주 도구로 사용. hifimeth는 선택적 Task로 5mC 결과 cross-validation에만 활용 (기본 파이프라인에서 제외, 입력 파라미터 `run_hifimeth_qc: Boolean = false`로 조건부 실행).

### Task 3b: `jasmine_pileup`
```
목적: jasmine 처리된 BAM에서 5mC / 5hmC / 6mA per-site 확률 집계
입력: merged.bam (jasmine MM/ML 업데이트된), ref_fasta, out_prefix
출력:
  - 5mC_pileup.bed (CpG site별 5mC 확률)
  - 5hmC_pileup.bed (CpG site별 5hmC 확률) ← 핵심 신규 출력
  - 6mA_pileup.bed (non-CpG 6mA 위치)
컨테이너: quay.io/pacbio/pb-cpg-tools:2.3.2
명령: aligned_bam_to_cpg_scores --bam merged.bam --ref ref.fasta --output-prefix out
```
> **5hmC 생물학적 맥락**:
> - Pastor et al. (2013) *Nature Reviews Genetics* 14:341 — 5hmC는 뇌, 간, 배아에서 높은
>   농도로 존재. DNA demethylation의 중간 산물이자 독립 에피마크.
> - Kriaucionis & Heintz (2009) *Science* 324:929 — 뉴런에서 5hmC 발견.
> - Drug A+B vs Drug A 비교에서 5hmC 감소 → Drug B가 TET 활성(및 ROS) 억제 경로를
>   통해 산화 관련 5mC→5hmC 전환을 차단하는지 확인 가능.

### Task 3c: `cpg_pileup` (기존 task 재사용 — 비교용)
```
재사용: workflows/wdl-common/wdl/tasks/cpg_pileup.wdl
목적: 기존 MM/ML 기반 5mC 결과와 jasmine 결과 cross-validation
비고: jasmine_pileup(3b)와 동일 BAM 사용 → 5mC 수치 일치 여부 QC
```

### Task 4a: `filter_and_classify`
```
목적: ipdSummary 결과에서 변형 분류 및 unknown 후보 추출
입력: modifications.csv, group_label
출력:
  - unknown_candidates.tsv (ipdRatio≥2.0, modScore≥30, coverage≥25,
    5mC@CpG 및 m6A canonical context 제외)
  - known_mod_summary.tsv (참고용)
  - GGG_flagged.tsv (8-oxoG 우선 후보 — GGG context 표시)
컨테이너: busco (python3)
```

### Task 4b: `compare_known_mods`
```
목적: 4개 그룹 5mC / 5hmC / 6mA 수준 비교
입력: 4개 그룹 jasmine_pileup BED
출력: differential_known_mods.tsv
  (per-site: 5mC_delta, 5hmC_delta, 5hmC_5mC_ratio 변화 — Drug A에서 ratio 증가 예상)
컨테이너: busco (python3)
핵심 지표: 5hmC/(5mC+5hmC) 비율 — TET 활성도 대리 지표
```

### Task 4c: `cluster_unknown_modifications` (신규)
```
목적: unknown modification 후보를 IPD 프로파일 + 서열 context로 클러스터링
      → 몇 종류의 구조적으로 다른 변형이 존재하는지 추정
입력: 4개 그룹 unknown_candidates.tsv (통합), modifications.csv (IPD 상세값)
출력:
  - clustered_modifications.tsv (후보 × 클러스터 ID)
  - cluster_summary.tsv (클러스터별: 크기, 평균 ipdRatio, 대표 context, 예상 변형 유형)
  - umap_plot.pdf (UMAP 시각화 — 클러스터 분포)
  - cluster_ipd_heatmap.pdf (클러스터별 ±3bp IPD 프로파일 히트맵)
컨테이너: quay.io/biocontainers/busco:5.7.1--pyhdfd78af_1
  (python3 + sklearn + umap-learn + seaborn 포함)
특징 벡터 구성:
  - IPD ratio: position-2, -1, 0, +1, +2 (5-dimensional local profile)
  - 서열 context: 3-mer at modification position (4-hot encoding)
  - strand symmetry score: |forward_ipdRatio - reverse_ipdRatio|
  - modificationScore (정규화)
알고리즘:
  1. UMAP (n_components=2) → 2D embedding
  2. HDBSCAN (min_cluster_size=50) → 클러스터 할당 (노이즈 포인트 별도 표시)
  3. 클러스터 대표 k-mer 추출 → 알려진 변형 IPD 패턴과 비교
     * 높은 ipdRatio + GGG context → 8-oxoG 후보 클러스터
     * 대칭 + CpG context → 5hmC 누락 클러스터
     * 비대칭 + 퓨린 → FapyG/FapyA 후보 클러스터
```
> **clustering 근거**:
> - Schadt et al. (2013) *Genome Research* — 각 변형 유형은 수정 위치와 ±수 염기에 걸쳐
>   고유한 IPD kinetic 지문(fingerprint)을 가짐. 이를 클러스터링하면 modification type 구분 가능.
> - Beaulaurier et al. (2019) *Nature Reviews Microbiology* 17:157 — 세균 에피게놈에서
>   IPD 클러스터링으로 m6A, m4C, m5C를 구분한 선례.
> - HDBSCAN 선택 이유: unknown modification의 클러스터 수가 미지이므로 k 사전 지정 불필요.

### Task 5: `multigroup_comparison`
```
목적: unknown + known 변형 통합 4-group 비교
입력: 4개 그룹 unknown_candidates + differential_known_mods + clustered_modifications.tsv
출력:
  - categorized_modifications.tsv (A_only / rescued_by_B / B_protective / synergistic)
  - rescue_candidates.tsv ← Drug A↑ → Drug A+B에서 복구된 unknown 변형
  - 5hmC_rescue.tsv ← Drug A에서 증가한 5hmC가 Drug B로 감소하는 위치
  - venn_summary.tsv
  - cluster_group_distribution.tsv ← 클러스터별로 4개 그룹에서의 분포 (특정 변형 유형이 Drug A에서만 나타나는지)
컨테이너: busco (python3)
```

### Task 6: `ipd_to_bigwig`
```
목적: per-base ipdRatio → BigWig (IGV / UCSC Genome Browser 시각화)
입력: modifications.csv, chrom_sizes
출력: {group}_{comparison}_ipdRatio.bw (각 비교별 BigWig)
추가: 5hmC_pileup.bw, 5mC_pileup.bw (jasmine_pileup에서 변환)
컨테이너: quay.io/biocontainers/ucsc-bedgraphtobigwig:469--h2b7bef9_0
추가 출력: IGV 세션 파일(.xml) — hybrid_bam + ipdRatio.bw + 5hmC.bw 통합 뷰
```

### Task 7: `motif_analysis`
```
목적: rescue 후보 및 A_only 후보 주변 서열 motif 탐색
입력: rescue_candidates.tsv, ref_fasta
출력: meme_rescued/, streme_out/
컨테이너: quay.io/biocontainers/meme:5.5.5--py311pl5321h6f6fdc4_1
명령:
  meme candidate_seqs.fa -dna -oc meme_out -nmotifs 10 -minw 6 -maxw 20 -mod zoops
  streme --p rescued.fa --n Aonly.fa --dna --oc streme_out
```
> 예상 motif: GGG/GGGG (8-oxoG), CpG (5hmC), ARE (Antioxidant Response Element: TGA[CG]NNNGC)

### Task 8: `annotate_modifications`
```
목적: 후보 위치에 유전체 context 부착
입력: categorized_modifications.tsv, gencode.vM36 GTF, CpG island TSV
출력: annotated_modifications.tsv (gene_name, feature_type, cpg_island, dist_to_TSS, is_promoter)
컨테이너: quay.io/biocontainers/bedtools:2.31.1--h1365805_1
재사용: hifi-wdl-resources/GRCm39/gencode.vM36.annotation.sorted.gtf.bgz
재사용: hifi-wdl-resources/GRCm39/cpgIslandExt.sorted.mm39.tsv
```

### Task 9: `pathway_enrichment`
```
목적: 변형 유전자 GO/KEGG enrichment
입력: annotated_modifications.tsv (promoter/exon 위치)
출력: go_enrichment.tsv, kegg_enrichment.tsv, enrichment_plot.pdf
컨테이너: quay.io/biocontainers/bioconductor-clusterprofiler (R)
우선 확인 pathway: BER (base excision repair), Nrf2/KEAP1, NF-κB, TET/DNA demethylation
```
> **문헌**: Huang et al. (2009) *Nature Protocols* — clusterProfiler GO/KEGG enrichment 방법론

### Task 10: `modification_report`
```
목적: 전체 분석 통합 HTML 보고서
입력: 모든 Task 출력
출력:
  - base_mod_report.html (대화형)
  - final_candidates.tsv
  - rescue_priority_list.tsv (생물학적 검증 우선순위)
보고 내용:
  - 요약 카드: unknown 후보 수 / rescue 후보 수 / 5hmC 변화 유전자 수 / 추정 변형 종류 수(클러스터)
  - Venn 다이어그램 (4 그룹 비교)
  - **UMAP 클러스터 plot** (unknown modification 종류 시각화)
  - 클러스터별 IPD 프로파일 히트맵 + 예상 변형 유형 레이블
  - 상위 rescue 후보 테이블 (ipdRatio, 5hmC_delta, cluster_id, gene, feature, motif)
  - 5hmC / 5mC 변화량 산포도 (Drug A vs Drug A+B)
  - MEME 최상위 motif
  - Pathway enrichment dot plot
  - 검증 제안 (DIP-seq / hMeDIP-seq / LC-MS/MS / dot blot 타겟)
  - hifimeth QC: run_hifimeth_qc=true일 경우 5mC concordance 테이블 추가
```

---

## 워크플로우 입력

```wdl
input {
  String experiment_id

  Array[File] control_bams      # 3개 (/mnt/ 외장 드라이브)
  Array[File] drug_a_bams       # 3개 (ROS 유도제)
  Array[File] drug_b_bams       # 3개 (항산화제)
  Array[File] drug_ab_bams      # 3개 (복합)

  File ref_fasta                # mouse_GRCm39.fasta
  File ref_fai
  File annotation_gtf_bgz       # gencode.vM36 (기존 리소스)
  File annotation_gtf_tbi
  File cpg_island_tsv           # cpgIslandExt.sorted.mm39.tsv (기존)

  Float   min_ipd_ratio    = 2.0
  Int     min_mod_score    = 30
  Int     min_coverage     = 25
  Float   ipd_pvalue       = 0.001

  # 선택적: hifimeth QC (5mC cross-validation)
  # hifimeth는 conda 패키지 없음 → 별도 이미지 빌드 필요 시만 활성화
  Boolean run_hifimeth_qc  = false

  RuntimeAttributes default_runtime_attributes
}

# 비고: Task 4c (cluster_unknown_modifications)는 python3 + umap-learn + hdbscan 필요
# busco 컨테이너에 없을 경우 별도 Dockerfile 준비:
#   FROM quay.io/biocontainers/busco:5.7.1--pyhdfd78af_1
#   RUN pip install umap-learn hdbscan seaborn scikit-learn
```

---

## 실행 스크립트 (기존 run_*.sh 패턴 동일)

```bash
#!/bin/bash
source /home/ygkim/program/anaconda3/etc/profile.d/conda.sh
conda activate hifi-human-wgs
OUTPUT_DIR="/data_4tb/hifi-human-wgs-wdl-custom/output/base_mod_discovery"
cd /data_4tb/hifi-human-wgs-wdl-custom || exit 1
miniwdl run workflows/base_modification_discovery.wdl \
    --input base_mod_discovery.inputs.json \
    --dir "${OUTPUT_DIR}" --verbose \
    2>&1 | tee "${OUTPUT_DIR}.log"
```

---

## 검증 방법

1. `miniwdl check workflows/base_modification_discovery.wdl`
2. Task 1→2a: 1개 샘플로 jasmine 실행 후 BAM MM/ML 태그에 `C+m` (5mC), `C+h` (5hmC) 존재 확인
   ```bash
   samtools view modified.bam | head -5 | tr '\t' '\n' | grep -E '^(MM|ML):'
   ```
3. Task 3b: 5hmC_pileup.bed에 0이 아닌 5hmC 확률 존재 여부 확인
4. Task 5 rescue 후보 수 sanity check (Drug A에서 ≥100개 없으면 파라미터 재검토)
5. 최종: `rescue_priority_list.tsv` 상위 20개 IGV 수동 확인

---

## 구현 순서

1. `base_modification.wdl` (tasks) — Task 1, 2a, 2b, 3a, 3b, 3c(재사용), 4a, 4b, 5~10
2. `base_modification_discovery.wdl` (workflow) — scatter/call/output
3. `base_mod_discovery.inputs.template.json`
4. `run_base_mod_discovery.sh`
5. miniwdl check
6. BAM 경로 확정 후 inputs.json 채우기

---

## 구현 시 변경된 사항

실제 `base_modification.wdl` / `base_modification_discovery.wdl` 구현 시 위 계획에서 다음 부분이 변경되었습니다. 문헌 근거와 설계 의도(위 본문)는 그대로 유효하지만, **실행 가능한 구체적 스펙은 아래 변경사항 + [docs/base_modification_discovery.md](./base_modification_discovery.md) 가 최종**입니다.

| 항목 | 계획 (이 문서) | 실제 구현 |
|---|---|---|
| Task 1 컨테이너 | `quay.io/pacbio/pbmm2:1.13.0` | `~{runtime_attributes.container_registry}/pbmm2@sha256:...` (기존 humanwgs pbmm2.wdl과 동일 SHA 고정 이미지, `container_registry`를 `quay.io/pacbio`로 설정 필요) |
| Task 3b 도구 | `pb-cpg-tools aligned_bam_to_cpg_scores` | **`modkit pileup`** — `aligned_bam_to_cpg_scores`는 5mC만 지원하고 5hmC를 지원하지 않는다는 사실을 구현 중 확인하여 변경 |
| Task 4c 출력 | `cluster_ipd_heatmap.pdf` 포함 | heatmap 미구현 (UMAP plot만 생성). IPD 5-dim 위치별 프로파일 대신 ipdRatio + modScore + GGG flag + 3-mer encoding을 특징 벡터로 사용 |
| Task 5 출력 | `cluster_group_distribution.tsv` 포함 | 미구현 (categorized/rescue/5hmC_rescue/venn 4종만 생성) |
| Task 6 | 5hmC/5mC BigWig + IGV 세션 xml 포함 | ipdRatio BigWig만 구현 (5hmC/5mC BigWig, IGV xml 미구현) |
| Task 9 도구 | `bioconductor-clusterprofiler` (R) | **`gseapy`** (Python) — busco/python3 기반 task 스타일과 통일하기 위해 R 대신 채택 |
| Task 10 출력 | `final_candidates.tsv` 포함 | 미구현 (`report_html`, `rescue_priority_list.tsv` 2종만 생성) |
| hifimeth QC | `run_hifimeth_qc` 플래그로 조건부 실행 예정 | 파라미터만 예약되어 있고 실제 cross-validation 로직은 미구현 |
| kinetics 태그 표기 | `ip/pw` (정렬 전후 구분 없이 표기) | raw HiFi BAM은 `fi`/`fp`/`ri`/`rp`, pbmm2 정렬 후에야 `ip`/`pw`로 재정리됨 — 시점에 따라 태그명이 다름 (상세: [base_modification_discovery.md](./base_modification_discovery.md) 실행 전 준비 절) |

추가로, `jasmine:2.0.0`, `kineticstools:2.5.0`, `modkit:0.3.3` 컨테이너 태그의 실제 publish 여부는 구현 시점에 검증하지 못했습니다. 최초 실행 전 1개 샘플로 스모크 테스트를 권장합니다.
