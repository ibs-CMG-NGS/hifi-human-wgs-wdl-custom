# Base Modification Discovery

PacBio HiFi kinetics(IPD/PW) 데이터를 이용해 oxidative stress 약물이 마우스 게놈에 유발하는
**알려진 변형(5mC, 5hmC, 6mA)** 과 **알려지지 않은 DNA 염기 변형(unknown modification)** 을
genome-wide로 탐색하고, 약물 간 비교를 통해 에피게놈 변화 지도를 구축하는 파이프라인입니다.

## 실험 설계

| 그룹 | 역할 | 예상 변형 신호 |
|---|---|---|
| Control | 기준선 | baseline |
| Drug A | ROS(oxidative stress) 유도제 | 8-oxoG↑, FapyG/FapyA↑, 5hmC↑ (TET 활성), 5mC→5hmC 전환↑ |
| Drug B | 항산화제 | baseline 또는 소폭 ↓ |
| Drug A+B | 복합 처치 | Drug A 효과 부분 구제(rescue) 예상 |

- 각 그룹 3마리 → pooling 후 ~90x 유효 커버리지 (개별 샘플은 25~35x)
- 종: 마우스 (GRCm39)
- 시퀀싱 시 kinetics 정보 포함 필수 (raw HiFi BAM 기준 `fi`/`fp`/`ri`/`rp` 태그 — 정렬 전 PacBio HiFi reads의 표준 kinetics 태그명. 정렬 전 BAM은 `samtools view`로 직접 확인 가능)

### 핵심 비교 축

```
① Drug A vs Control          → ROS 유도 변형 (induction)
② Drug A+B vs Drug A         → Rescue 변형 ← 가장 중요
③ Drug B vs Control          → 항산화 단독 효과
④ Drug A+B vs Control        → 복합 처치 순효과
```

## 워크플로우 개요

```
workflows/base_modification_discovery.wdl
workflows/wdl-common/wdl/tasks/base_modification.wdl   (task 정의)
```

### 분석 단계

```
[1. 정렬 — kinetics 보존, 12샘플 병렬]
  pbmm2_align_kinetics        raw HiFi BAM(fi/fp/ri/rp) → GRCm39 정렬 (정렬 시 게놈 strand 기준 ip/pw로 재정리, MM/ML 태그도 유지)

[2. 변형 calling + pooling]
  jasmine_modification_calling  kinetics → ML 모델로 5mC/5hmC/6mA 확률 산출 (MM/ML 갱신)
    → merge_group_bams          그룹별 3 replicate 병합 (~90x)

[3. 변형 탐지 — known + unknown]
  ipd_summary                 treated vs control two-sample IPD 비교 → unknown 포함 전체 변형
  jasmine_pileup (modkit)      5mC / 5hmC / 6mA per-site 확률 BED
  cpg_pileup (기존 task 재사용) 5mC 결과 cross-validation용 QC

[4. 필터링 / 비교 / 클러스터링]
  filter_and_classify          unknown 후보 추출 (known CpG 변형 제외, 8-oxoG GGG context 표시)
  compare_known_mods           4개 그룹 5mC/5hmC/6mA 수준 비교, 5hmC/(5mC+5hmC) 비율 계산
  cluster_unknown_modifications  UMAP + HDBSCAN으로 unknown 변형을 종류별로 클러스터링

[5. 통합 비교 — rescue 분석]
  multigroup_comparison        Drug A에서 증가 → Drug A+B에서 복구된 변형(rescue) 추출
                                5hmC_rescue (TET 활성 관련 5hmC 복구) 별도 추출

[6~9. Downstream 분석]
  ipd_to_bigwig                IPD ratio → BigWig (IGV 시각화)
  motif_analysis (MEME/STREME) rescue 후보 주변 서열 motif 탐색
  annotate_modifications        Gencode GTF + CpG island로 유전체 context 부착
  pathway_enrichment (gseapy)   변형 유전자 GO/KEGG enrichment

[10. 리포트 — 항상 실행]
  modification_report          통합 HTML 리포트 + rescue 우선순위 목록 생성
```

### 핵심 개념

- **Two-sample ipdSummary**: 약물처치 BAM과 control BAM을 직접 비교하여 IPD ratio 산출. In-silico 모델 대신 생물학적 control을 사용해 배경 노이즈를 줄임
- **jasmine**: kinetics(ip/pw)를 ML 모델에 적용해 5mC(`C+m`)와 5hmC(`C+h`)를 구분하는 PacBio 공식 도구. `aligned_bam_to_cpg_scores`는 5mC만 지원하므로 5hmC에는 `modkit pileup` 사용
- **Rescue 분석**: Drug A에서 증가한 변형이 Drug A+B(항산화제 동시 투여)에서 사라지거나 감소하는지 확인 → 항산화제의 보호 효과를 변형 수준에서 직접 검증
- **Unknown modification 클러스터링**: ipdRatio 프로파일 + 3-mer 서열 context를 특징 벡터로 UMAP 임베딩 후 HDBSCAN으로 그룹화. 클러스터 수(k)를 사전에 지정할 필요가 없어 알려지지 않은 변형 종류 수 추정에 적합. GGG context + 높은 IPD인 클러스터는 8-oxoG 후보로 분류
- **hifimeth QC (선택)**: 5mC만 지원하는 별도 deep-learning 도구로, jasmine 5mC 결과와 cross-validation하고 싶을 때만 `run_hifimeth_qc=true`로 활성화 (현재 미구현, 파라미터만 예약됨)

## 입력 파라미터

### 필수

| 파라미터 | 설명 |
|---|---|
| `experiment_id` | 실험 식별자 (보고서 헤더에 사용) |
| `control_bams` | Control 그룹 raw HiFi BAM 3개 (Array, kinetics 태그 `fi`/`fp`/`ri`/`rp` 포함 — 정렬 전 원본) |
| `drug_a_bams` | Drug A(ROS 유도제) 그룹 raw HiFi BAM 3개 |
| `drug_b_bams` | Drug B(항산화제) 그룹 raw HiFi BAM 3개 |
| `drug_ab_bams` | Drug A+B 복합 처치 그룹 raw HiFi BAM 3개 |
| `control_ids` / `drug_a_ids` / `drug_b_ids` / `drug_ab_ids` | 각 그룹 replicate 샘플 ID (Array[String], BAM 순서와 일치해야 함) |
| `ref_fasta` | GRCm39 레퍼런스 게놈 FASTA |
| `ref_fai` | `ref_fasta`의 FAI 인덱스 |
| `annotation_gtf_bgz` | Sorted + bgzipped Gencode GTF (vM36) |
| `annotation_gtf_tbi` | Tabix 인덱스 (.tbi) |
| `cpg_island_tsv` | CpG island 좌표 TSV |
| `default_runtime_attributes` | 런타임 설정 (`container_registry`는 `quay.io/pacbio` 필요 — pbmm2 이미지 위치) |

### 선택 — 파라미터 튜닝

| 파라미터 | 기본값 | 설명 |
|---|---|---|
| `min_ipd_ratio` | 2.0 | unknown 후보 판정 최소 IPD ratio |
| `min_mod_score` | 30 | unknown 후보 판정 최소 modificationScore |
| `min_coverage` | 25 | ipdSummary 최소 site coverage |
| `ipd_pvalue` | 0.001 | ipdSummary p-value cutoff |
| `extract_window_bp` | 250 | motif 분석 시 변형 위치 좌우 서열 추출 범위 (bp) |
| `annotation_window_bp` | 2000 | promoter 분류용 TSS 기준 window (bp) |
| `min_cluster_size` | 50 | HDBSCAN 최소 클러스터 크기 |
| `run_hifimeth_qc` | false | hifimeth 5mC QC cross-validation 활성화 여부 (예약, 미구현) |

## 출력 파일

| 출력 | 설명 |
|---|---|
| `report_html` | 통합 HTML 리포트 (요약 통계, rescue 후보, 클러스터, enrichment) |
| `rescue_priority_list` | Drug A→A+B rescue 후보 우선순위 TSV (검증 실험 대상 선정용) |
| `rescue_candidates` | unknown modification rescue 후보 전체 목록 |
| `fivehmC_rescue` | 5hmC 증가→복구 위치 (TET 활성 변화 추정) |
| `categorized_modifications` | 모든 unknown 후보의 카테고리 분류 (A_only_persistent / rescued_by_B) |
| `differential_known_mods` | 4그룹 5mC/5hmC/6mA per-site 비교 테이블 |
| `clustered_modifications` / `cluster_summary` / `umap_plot` | unknown 변형 클러스터링 결과 (TSV + UMAP PDF) |
| `bw_drug_a_vs_ctrl` 등 | IPD ratio BigWig (IGV 시각화용) |
| `meme_results` / `streme_results` | motif 분석 결과 (tar.gz) |
| `annotated_tsv` | 변형 위치별 유전자/CpG island/promoter 어노테이션 |
| `go_enrichment` / `kegg_enrichment` | 변형 유전자 pathway enrichment 결과 + plot |
| `ctrl_merged_bam` 등 | 그룹별 pooled BAM (jasmine MM/ML 태그 포함, IGV 확인용) |
| `ctrl_cpg_qc_bed` / `drug_a_cpg_qc_bed` | 기존 `cpg_pileup` task 기반 5mC QC (jasmine 결과와 대조용) |

## 실행 방법

### 1. inputs.json 작성

`base_mod_discovery.inputs.template.json`을 복사해 `base_mod_discovery.inputs.json`을 만들고
`/mnt/` 하위 실제 BAM 경로와 샘플 ID를 채웁니다.

```bash
cp base_mod_discovery.inputs.template.json base_mod_discovery.inputs.json
# base_mod_discovery.inputs.json 을 열어 control_bams / drug_a_bams / ... 경로 수정
```

### 2. 실행

```bash
conda activate hifi-human-wgs
cd /data_4tb/hifi-human-wgs-wdl-custom
bash run_base_mod_discovery.sh
```

또는 직접 실행:

```bash
miniwdl run workflows/base_modification_discovery.wdl \
  --input base_mod_discovery.inputs.json \
  --dir output/base_mod_discovery --verbose \
  2>&1 | tee base_mod_discovery.log
```

### 3. 결과 확인

```
output/base_mod_discovery/_LAST/out/
  report_html/                  ← 통합 HTML 리포트 (가장 먼저 볼 파일)
  rescue_priority_list/         ← 검증 실험 우선순위 목록
  bw_drug_a_vs_ctrl/             ← IGV에 로드해서 IPD ratio 확인
  umap_plot/                    ← unknown 변형 클러스터 시각화
```

## 선행 조건 / 실행 전 준비

1. **kinetics 포함 시퀀싱**: 원본(정렬 전) HiFi BAM에 kinetics 태그가 존재해야 합니다. 정렬 전 BAM에서는 `fi`(forward IPD) / `fp`(forward PulseWidth) / `ri`(reverse IPD) / `rp`(reverse PulseWidth) 형태로 존재하며, `pbmm2_align_kinetics` task가 정렬하면서 게놈 strand 기준으로 재정리되어 `ip`/`pw` 태그가 됩니다. 즉 raw BAM에서 `ip`/`pw`를 찾으면 안 되고 `fi`/`fp`/`ri`/`rp` 존재 여부로 확인해야 합니다.
   ```bash
   samtools view raw_hifi.bam | head -1 | tr '\t' '\n' | grep -E '^(fi|fp|ri|rp):'
   ```
   일반 humanwgs 파이프라인의 `pbmm2.wdl`은 기본적으로 `strip_kinetics=true`로 이 태그를 정렬 시 제거하므로, 이 파이프라인은 별도의 `pbmm2_align_kinetics` task로 직접 정렬합니다 — **반드시 raw HiFi BAM(정렬 전, kinetics 포함)을 입력으로 사용해야 합니다.**
2. **그룹/replicate 순서 일치**: `control_bams`와 `control_ids` 등 BAM 배열과 ID 배열의 순서가 반드시 일치해야 합니다.
3. **컨테이너 검증 필요**: `jasmine:2.0.0`, `kineticstools:2.5.0`, `modkit:0.3.3` 등 일부 컨테이너 태그는 실제 publish 여부를 사전에 확인하지 못했습니다. 전체 12샘플 실행 전에 **그룹당 1개씩(control/drug_a/drug_b/drug_ab 각 1개, 총 4개) BAM**으로 먼저 파일럿 실행을 권장합니다. 샘플 1개만으로는 `pbmm2_align_kinetics` → `jasmine_modification_calling`까지만 검증되고 `ipd_summary`(treated vs control 비교), `merge_group_bams` 로직은 확인되지 않으므로, 4그룹 각 1개 BAM으로 전체 워크플로우를 한 번 통과시켜보는 것이 더 안전합니다.
   ```bash
   # 1) 파일럿용 inputs.json: 각 그룹 배열에 BAM 1개 / ID 1개만 채워서 실행
   # 2) 정렬+jasmine 출력 확인
   samtools view modified.bam | head -5 | tr '\t' '\n' | grep -E '^(MM|ML):'
   # 출력에 C+m (5mC), C+h (5hmC) 가 보이면 정상
   # 3) ipd_summary, multigroup_comparison까지 에러 없이 끝까지 도는지 확인
   ```
4. **디스크 용량**: 샘플당 raw BAM, 정렬 BAM, jasmine BAM이 모두 보존되며 그룹별 pooled BAM(~90x)도 추가로 생성되므로 충분한 디스크 여유가 필요합니다.
5. **`container_registry` 설정**: `pbmm2_align_kinetics` task는 `~{runtime_attributes.container_registry}/pbmm2@sha256:...` 형태로 이미지를 가져오므로 `default_runtime_attributes.container_registry`를 `quay.io/pacbio`로 설정해야 합니다 (다른 task들은 `quay.io/biocontainers` 공개 이미지를 직접 사용).

## 레퍼런스 데이터 (공용)

| 파일 | 경로 |
|---|---|
| GRCm39 레퍼런스 | `/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/mouse_GRCm39.fasta` |
| GRCm39 FAI | `/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/mouse_GRCm39.fasta.fai` |
| Gencode vM36 GTF (bgzip) | `/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/gencode.vM36.annotation.sorted.gtf.bgz` |
| Gencode vM36 GTF index | `/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/gencode.vM36.annotation.sorted.gtf.bgz.tbi` |
| CpG island TSV | `/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/cpgIslandExt.sorted.mm39.tsv` |

## 검증 가이드 (분석 후 sanity check)

1. `cluster_unknown_modifications` 결과에서 클러스터 수가 0이면 `min_cluster_size`를 낮춰 재실행
2. `rescue_candidates` 수가 너무 적으면(<100) `min_ipd_ratio` / `min_mod_score`를 완화해 재검토
3. `5hmC_rescue.tsv`가 비어 있으면 `compare_known_mods`의 5hmC 임계값(코드 내 0.1/-0.05)을 조정
4. 최종적으로 `rescue_priority_list` 상위 20개를 IGV에서 `bw_drug_a_vs_ctrl` / `ctrl_merged_bam` 등과 함께 수동 확인
5. 실험적 검증 권장: 8-oxoG 후보 → 항-8-oxoG dot blot, 5hmC rescue → hMeDIP-seq/TAB-seq, FapyG/FapyA 후보 → LC-MS/MS
