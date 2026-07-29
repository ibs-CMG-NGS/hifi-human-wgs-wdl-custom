# H2O2 성상세포 WGS 발표 — 그림(Plot) 배치 가이드

**작성일:** 2026-07-26
**짝 문서:** [`H2O2_PPT_source_material_2026-07-25.md`](./H2O2_PPT_source_material_2026-07-25.md)(슬라이드 콘텐츠), [`H2O2_PPT_presenter_notes_2026-07-25.md`](./H2O2_PPT_presenter_notes_2026-07-25.md)(발표자 노트)
**용도:** 어느 슬라이드에 **실제 결과 그림(PDF)** 을 넣을지, 각 그림이 **네거티브 데이터(신호 없음)인지 포지티브 데이터(그라디언트 있음)인지**, 그리고 발표 상 어떻게 배치·강조할지 지정한다. 모든 경로는 `/data_4tb/shared/wgs-tertiary-analysis/` 기준.

- 기본 그림 폴더: `h2o2_analysis_results/cohort/`
- joint 버전 폴더: `h2o2_analysis_results_joint/cohort/`
- 샘플별 ASM: `h2o2_analysis_results/asm/`

---

## 0. 한눈에 보는 배치 표 (슬라이드 → 그림 → 데이터 성격)

| 슬라이드 | 넣을 그림 파일 | 데이터 성격 | 발표 상 역할 |
|---|---|---|---|
| 1-2 파이프라인 개요 | (그림 없음 — AI 모식도/플로우차트 생성) | — | 개념도 |
| 2-1 VCF 포맷 | (그림 없음 — 텍스트 레코드 예시 캡처) | — | 포맷 설명 |
| **2-5b Burden 결과** | `variant_burden_plots.pdf` | **네거티브 (핵심)** | "개별 유전자엔 신호 없음" 증명 |
| **2-5b Burden 결과** | `sv_burden_plots.pdf` (단일샘플) | **네거티브** | SV도 신호 없음 |
| **2-5b Burden 결과 (함정)** | `h2o2_analysis_results_joint/cohort/sv_burden_plots.pdf` | **네거티브(아티팩트)** | joint 112개는 다카피 아티팩트임을 대조로 보여줌 |
| 2-5b Burden 결과(보조) | `geneset_burden_variant_plots.pdf` / `geneset_burden_sv_plots.pdf` | 네거티브 | GO 세트도 무의미 |
| 2-4 DMR(선택) | `dmr_A_vs_D_plots.pdf` (대표 1개) | 약한/보조 | 메틸화 차이 예시(p값 한계 명시) |
| **3-4 het_hom_peak_frac** | `summary_compare_plots.pdf` (4패널 중 hom-peak% 패널) | **포지티브 (핵심)** | 그라디언트 첫 등장 |
| 3-3 VAF 개념(보조) | `summary_compare_plots.pdf` (median_het_vaf 패널) | **주의: 아티팩트 있음** | median VAF는 참고만, 핵심 아님 |
| **3-5 5개 검증 종합** | `heterogeneity_consensus_summary.pdf` | **포지티브 (최종 결론)** | 발표의 클라이맥스 |
| 3-5 방법 2 | `summary_compare_plots_discordance_rank.pdf` | 포지티브 | 불일치율 + 순위 일관성 |
| 3-5 방법 4 | `hap_imbalance_plots.pdf` | **포지티브 (유일한 완전 독립)** | haplotype 방향성 |
| 3-5 참고지표 | `hp_consistency_plots.pdf` | 포지티브(단, 비독립) | n_het_snv 상관 산점도로 재해석 근거 |
| 3-5 방법 5 (선택) | `asm/{sample}.asm_plots.pdf` | 포지티브(비독립) | 메틸화 차원 예시 |

---

## 1. 네거티브 데이터 배치 — "범인 유전자는 없다" (Slide 2-5b)

**배치 전략**: 이 슬라이드는 발표의 **전환점**이다. "개별 유전자를 봤지만 신호가 없어서 게놈 전역 판세로 갔다"는 논리를 그림으로 증명. 한 슬라이드에 그림 2~3개를 병치.

### 왼쪽/위: `variant_burden_plots.pdf` — SNV+indel (네거티브 핵심)
- **무엇이 보여야 하나**: 유전자별 p-value 분포 또는 volcano 형태. 대부분이 유의선(예: p=0.05, FDR=0.05) 위에 안 걸림. 최상위 B4galt3(p=0.039)조차 FDR=1.0.
- **캡션(권장)**: "80개 유전자 중 명목상 p<0.05는 1개(B4galt3)뿐, FDR 보정 후 유의 유전자 0개 — 개별 SNV/indel 층위엔 군을 가르는 신호가 없음"
- **발표 멘트**: "보시면 유의선을 넘는 게 사실상 없습니다. 하나 걸린 것도 다중검정 보정하면 사라집니다."

### 오른쪽/아래 1: `sv_burden_plots.pdf` (단일샘플) — SV (네거티브)
- **캡션**: "SV도 동일 — 278개 중 FDR<0.05 0개"

### 오른쪽/아래 2 (함정 대조): `h2o2_analysis_results_joint/cohort/sv_burden_plots.pdf` — joint SV
- **무엇이 보여야 하나**: joint에서는 112개가 유의선 위로 튀어 보임. 단 최상위 라벨이 Scgb1b10/Sirpb1b/Scgb2b11 등 **다카피 유전자족**.
- **캡션(중요)**: "joint calling에선 112개가 '유의'하게 보이지만 상위가 전부 다카피 유전자족(Scgb/Sirpb) = 정렬 아티팩트. 단일샘플(0개)을 신뢰"
- **발표 멘트**: "이건 함정입니다. 유의하게 보이지만 열어보면 전부 게놈에 사본이 여러 개인 유전자족이라 진짜 신호가 아닙니다."

### 보조(공간 있으면): `geneset_burden_variant_plots.pdf` / `geneset_burden_sv_plots.pdf`
- **캡션**: "산화스트레스 반응·DNA repair 유전자 세트로 묶어도 유의하지 않음(최소 FDR 0.22)"

> **배치 팁**: 네거티브 3개(variant/SV단일/SV joint)를 가로로 나란히 두고, joint만 빨간 테두리로 "아티팩트" 표시하면 대조가 극적. 또는 단일샘플 2개를 같이 두고 joint는 별도 "함정" 콜아웃 박스로.

---

## 2. 포지티브 데이터 배치 — 이질성 그라디언트 (Part 3)

### 2-1. het_hom_peak_frac & median_het_vaf → Slide 3-4 / 3-3

**파일**: `summary_compare_plots.pdf` (2×2 패널: total_pass_variants / total_pass_sv / **het_hom_peak_frac** / median_het_vaf)

- **het_hom_peak_frac 패널 (포지티브 핵심)** → **Slide 3-4**에 배치
  - 이 패널이 A<C<B<D 그라디언트를 boxplot으로 처음 보여주는 곳. 각 그룹 평균값 라벨 포함.
  - **캡션**: "이형접합 SNV 중 VAF>0.9 비율 — A<C<B<D로 단조 증가(Kruskal-Wallis H≈10). 처치가 셀수록 클론 우세↑"
  - **강조**: 이 패널만 잘라내 크게 쓰는 것을 권장(4패널 중 주인공).

- **median_het_vaf 패널 (주의: 아티팩트)** → **Slide 3-3**에 "참고"로만
  - **캡션(정직)**: "이형접합 VAF 중앙값 — 방향은 그라디언트와 같으나 B의 정확한 0.350000은 depth 양자화 수치 아티팩트로 확인됨. 핵심 지표 아님, 참고용"
  - **발표 멘트**: "VAF 중앙값도 비슷한 방향이지만, 이건 시퀀싱 depth가 만든 수치 아티팩트가 섞여 있어서 우리는 주력 지표로 쓰지 않습니다. 정직하게 빼둔 겁니다."

- **total_pass_variants / total_pass_sv 패널** → 사용 선택. 기술적 confound(배치·커버리지) 가능성 있어 "참고"로만. Slide 3-4 하단 작게 또는 생략.

### 2-2. 방법 2 (불일치율 + 순위 일관성) → Slide 3-5

**파일**: `summary_compare_plots_discordance_rank.pdf` (좌: 웰별 불일치율 barplot, 우: 지표별 순위 heatmap)
- **캡션**: "3반복 간 joint genotype 불일치율(좌)과 지표 간 순위 일관성(우) — 균질할수록 불일치↓, 여러 지표가 같은 군 순서를 가리킴(Kendall's W)"

### 2-3. 방법 4 (haplotype 방향성) → Slide 3-5 ★

**파일**: `hap_imbalance_plots.pdf`
- **데이터 성격**: **포지티브 + 유일한 완전 독립 증거** — 가장 강조할 그림.
- **캡션**: "PS(phase block) 단위 phase-0 allele 편향 — frac_biased가 A 0.12 < C 0.35 < B 0.63 < D 0.71. 개별 SNP가 아닌 haplotype 구간 일관성을 본 독립 증거(롱리드 고유)"

### 2-4. 참고지표 (HP 미배정률) → Slide 3-5

**파일**: `hp_consistency_plots.pdf` (미배정률 + n_het_snv 상관 산점도)
- **데이터 성격**: 포지티브 방향이나 **비독립**(이형접합 밀도와 rho=−0.993).
- **캡션(정직)**: "HP 미배정률도 같은 방향이나, n_het_snv와 거의 완벽 상관(rho=−0.993) — 독립 증거가 아니라 이형접합 밀도 감소의 재관측. 산점도로 그 근거를 명시"

### 2-5. 종합 (클라이맥스) → Slide 3-5 결론 ★★

**파일**: `heterogeneity_consensus_summary.pdf`
- **데이터 성격**: **포지티브 최종 결론** — 5개 방법 순위를 한 장 heatmap으로.
- **배치**: 발표의 마지막 결과 슬라이드에 가장 크게. 앞의 개별 그림들이 이 한 장으로 수렴함을 보여줌.
- **캡션**: "원리가 다른 5개 방법 모두 A<C<B<D를 가리킴 — 이 연구의 핵심 결론. (단 완전 독립은 방법 4뿐, 나머지는 방향 재현으로 신뢰 보강)"

### 2-6. ASM (선택) → Slide 3-5 방법 5 또는 부록
**파일**: `asm/H2O2-A01-ctrl1.asm_plots.pdf` 등 대표 1~2개
- **캡션**: "haplotype별 메틸화 차이 영역 예시 — n_asm_regions도 A>C>B>D(비독립, 이형접합 밀도와 연동)"

---

## 3. 발표 흐름 상 네거티브→포지티브 배치 논리 (스토리라인)

발표에서 그림 순서는 **네거티브 먼저, 포지티브 나중**이 설득력이 크다:

1. **(2-5b) 네거티브**: "개별 유전자를 봤다 → burden test 그림 → 신호 없음(FDR 소멸), joint는 아티팩트" → 청중에게 "그럼 뭐가 있나?"라는 긴장 유발.
2. **(3-4) 첫 포지티브**: "het_hom_peak_frac boxplot → 처음으로 A<C<B<D 그라디언트 등장" → 긴장 해소 시작.
3. **(3-5) 포지티브 누적**: 불일치율 → haplotype 방향성(독립!) → 종합 heatmap → "여러 창으로 봐도 같은 그림".
4. **정직성 배치**: median VAF·HP미배정·ASM은 각각 "아티팩트/비독립"임을 그림 캡션에 명시해, 포지티브를 부풀리지 않는 신뢰 확보.

> **핵심 원칙**: 포지티브 그림(그라디언트)마다 "이게 독립 증거인가"를 캡션에 표시. 방법 4(hap_imbalance)만 "완전 독립"으로, 나머지는 "방향 재현/보강"으로 라벨링해야 리뷰어 신뢰를 얻는다.

---

## 4. 그림 없이 텍스트/모식도로 처리할 슬라이드

| 슬라이드 | 처리 방법 |
|---|---|
| 1-1 HiFi 배경 | short vs long read 비교 **모식도**(AI 생성) |
| 1-2 파이프라인 개요 | **플로우차트 모식도**(소스 문서의 ASCII 흐름 기반 AI 생성) |
| 1-5b joint calling | 웰당 3반복 병합 **개념도** |
| 2-1 / 2-1b VCF 포맷 | 실제 레코드 한 줄 **텍스트 캡처** + 컬럼 주석 |
| 3-2/3-3 VAF 역추론 | 구슬통 비유 **개념도** + 역추론 수식 |
| 3-8 용어집 | 표(그림 없음) |

---

## 5. 파일 실경로 빠른 참조

```
# 네거티브 (Slide 2-5b)
h2o2_analysis_results/cohort/variant_burden_plots.pdf
h2o2_analysis_results/cohort/sv_burden_plots.pdf
h2o2_analysis_results_joint/cohort/sv_burden_plots.pdf      # joint (아티팩트 대조)
h2o2_analysis_results/cohort/geneset_burden_variant_plots.pdf
h2o2_analysis_results/cohort/geneset_burden_sv_plots.pdf

# 포지티브 (Part 3)
h2o2_analysis_results/cohort/summary_compare_plots.pdf                    # het_hom_peak_frac(3-4) + median_vaf(3-3)
h2o2_analysis_results/cohort/summary_compare_plots_discordance_rank.pdf   # 방법 2
h2o2_analysis_results/cohort/hap_imbalance_plots.pdf                      # 방법 4 (독립)
h2o2_analysis_results/cohort/hp_consistency_plots.pdf                     # 참고지표
h2o2_analysis_results/cohort/heterogeneity_consensus_summary.pdf          # 종합 결론 ★
h2o2_analysis_results/asm/<sample>.asm_plots.pdf                          # 방법 5 예시

# 보조
h2o2_analysis_results/cohort/dmr_A_vs_D_plots.pdf   # DMR 대표
h2o2_analysis_results/cohort/trgt_group_compare_plots.pdf
```
