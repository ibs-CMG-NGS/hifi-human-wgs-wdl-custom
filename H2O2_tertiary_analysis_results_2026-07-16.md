# H2O2 astrocyte WGS — wgs-tertiary-analysis 파이프라인 실행 결과 종합

**작성일:** 2026-07-16
**선행 문서:** [`H2O2_variant_count_qc_analysis_2026-07-14.md`](./H2O2_variant_count_qc_analysis_2026-07-14.md)(애드혹 분석, QC버그·이질성 가설 최초 발견), [`H2O2_key_findings_summary_2026-07-14.md`](./H2O2_key_findings_summary_2026-07-14.md)(그 요약)
**계획 문서:** [`wgs-tertiary-analysis/docs/H2O2_setup_plan_2026-07-14.md`](../shared/wgs-tertiary-analysis/docs/H2O2_setup_plan_2026-07-14.md), [`H2O2_group_comparison_v2_plan_2026-07-15.md`](../shared/wgs-tertiary-analysis/docs/H2O2_group_comparison_v2_plan_2026-07-15.md)
**작업 디렉토리:** `/data_4tb/shared/wgs-tertiary-analysis/` (아래 모든 상대경로의 기준)

이 문서는 애드혹 스크립트로 찾았던 이질성 신호를 정식 Snakemake 파이프라인(`wgs-tertiary-analysis`)으로 재현·검증한 결과를 정리한다. **가장 중요한 발견은 §0의 성별 교란**이며, 이후 모든 결과 해석에 영향을 준다. 각 분석 항목마다 핵심 도구/입출력/읽는 법/해석을 표준 형식으로 기재했다.

---

## 0. 핵심 발견: A/D군=수컷, B/C군=암컷 (완벽한 성별 교란)

**도구:** `mosdepth`(WDL 파이프라인 내 이미 실행됨, 별도 실행 불필요) — 염색체별 평균 depth 요약.

**입력:** 각 샘플 `{run_dir}/out/mosdepth_summary/*.txt` (예: `/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-A01-ctrl1/20260517_163205_humanwgs_singleton/out/mosdepth_summary/H2O2-A01-ctrl1.mosdepth.summary.txt`)

**출력:** 별도 저장 파일 없음(이 세션에서 1회성 python 스크립트로 계산, 문서에 결과만 기록). 재현하려면: 각 샘플의 `mosdepth_summary/*.txt`에서 `chr1`~`chr19` 평균(`mean` 컬럼)의 평균을 상염색체 기준값으로 하고, `chrX`/`chrY` 행의 `mean`을 이 값으로 나누면 된다.

**읽는 법:** `mosdepth.summary.txt`는 tab-delimited, 컬럼 `chrom, length, bases, mean, min, max`. `_region` 접미사가 붙은 중복 행은 무시(동일 값).

**해석:**

| 군 | 샘플 | chrX/상염색체 평균 | chrY/상염색체 평균 | 추정 성별 |
|---|---|---|---|---|
| A | A01-ctrl1, 2-A01, 3-A01 | 0.77~0.79 | 0.13~0.14 | **수컷(XY)** |
| B | B01, 2-B01, 3-B01 | 0.96~0.97 | 0.05 | **암컷(XX)** |
| C | C01, 2-C01, 3-C01 | 0.97~0.98 | 0.05 | **암컷(XX)** |
| D | D01, 2-D01, 3-D01 | 0.77~0.79 | 0.13~0.14 | **수컷(XY)** |

XY 개체는 chrX가 반수체라 상염색체 대비 ~0.5배가 기대치지만 실측은 ~0.78배(PAR 영역·매핑 특성 등으로 완전한 0.5는 아님) — 그래도 XX군(~0.97)과는 뚜렷이 구분된다. chrY는 XX 개체에서 거의 0(잔여 신호는 반복서열 교차매핑), XY에서 뚜렷이 검출됨. 3배치 전부 완벽하게 일관 — 우연이 아니다. 원본 WDL 입력(`inputs.json`)엔 성별 파라미터가 아예 없어 사전에 인지되지 않은 채 실행됐다.

**영향:** 성염색체(chrX/Y)를 포함한 모든 "군간비교"는 처치효과가 아니라 성별차를 반영할 위험이 있다. 아래 §1(SV burden)과 §5(DMR)에서 실제로 이 교란이 결과를 오염시킨 것을 확인했다. §2(샘플단위 요약통계)는 상염색체만으로 재확인해도 결과가 유지되어 이 교란의 영향을 받지 않았다.

---

## 1. SV burden (`cohort_sv_burden`) — 3라운드 조사 끝에 "single-sample이 신뢰할 결과"로 결론

**핵심 도구:** `svpack consequence`(SV의 유전자 영향 주석) → `scripts/sv_burden.R`(Fisher's exact/χ² + BH FDR, 효과크기 계산)

**목적:** 구조변이(SV, 결실/중복/삽입 등)가 특정 유전자를 파괴하는 빈도가 4개 처치군 사이에 다른지 유전자 단위로 검정 — H2O2 처치강도에 비례해 유전자 파괴성 SV 부담이 증가/감소하는 유전자가 있는지 탐색.

**입력:**
- VCF: 4군×3배치=12샘플의 `phased_sv_vcf`(single-sample, `filelist_h2o2.csv` 경유) 또는 joint-called `split_joint_structural_variant_vcfs`(`filelist_h2o2_joint.csv` 경유)
- 유전자 주석: `/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/Mus_musculus.GRCm39.113.chr.gff3`
- 중간산출물: `h2o2_analysis_results{,_joint}/structural_variants/{sample}.annotated_sv.tsv` (CHROM, POS, END, SVTYPE, SVLEN, BCSQ, GT)

**출력:**
- `h2o2_analysis_results/cohort/sv_burden_stats.csv` (single-sample, **최종 채택**)
- `h2o2_analysis_results/cohort/sv_burden_plots.pdf` (6페이지: SVTYPE 분포, SVLEN 분포, top20 raw-p 순위+효과크기, volcano, bubble)
- `h2o2_analysis_results_joint/cohort/sv_burden_stats.csv` (joint, **참고용 — 신뢰 불가**, §1.2 참고)

**출력 데이터 읽는 법 (`sv_burden_stats.csv` 컬럼):**

| 컬럼 | 의미 |
|---|---|
| `gene` | 유전자 심볼 |
| `pvalue` | 보정 전(raw) 검정 p-value (4군이라 χ² 검정, 셀 카운트<5면 시뮬레이션 p-value) |
| `odds_ratio` | 2군 비교시에만 값 있음(4군은 NA) |
| `n_A`..`n_D` | 해당 군에서 이 유전자에 SV(sv:cds/sv:utr consequence, ≤1Mb)를 가진 샘플 수 (군당 최대 3) |
| `max_prop_diff` | 군간 침투율(n/3) 최대-최소 차, 0~1 |
| `group_max`/`group_min` | 침투율 최고/최저 군 |
| `odds_ratio_extreme` | 최고/최저 군 사이 Haldane-Anscombe(+0.5) 보정 오즈비 — 카운트 0이어도 항상 유한값, fold-change처럼 해석 |
| `cramers_v` | 전체 군에 대한 연관성 강도(효과크기), 0~1 |
| `fdr` | BH 보정 후 p-value — **이 값으로 유의성 판단** |

**결과 해석:**

| 라운드 | 조치 | joint-called VCF | single-sample VCF |
|---|---|---|---|
| 0 (원본) | — | 582/1694 유의(FDR<0.05), 84%가 단일 9.8Mb duplication 아티팩트 | 최저 FDR 0.29~1.0 (0개 유의) |
| 1 | SVLEN>1Mb 제외 (`--max-svlen`) | 704/993 유의, 64%가 chrX/Y | 0/278 |
| 2 | chrX/Y 제외 (`--exclude-chroms`) | 218/479 유의, 다중카피 유전자군(Vmn1r/Mup/Pramel/Ear/Gm 등)이 다수 | 0/278 |
| 3 | 위 유전자군 접두어 제외 (`--exclude-gene-prefixes`) | **112/277 유의** — Scgb/Sirpb/Or(후각수용체)/Ang 등 또 다른 다중카피군으로 이동(두더지잡기) | **0/278 (일관됨)** |

joint calling은 마우스 게놈에 광범위하게 퍼진 segmental duplication/다중카피 유전자군의 **개체별 진짜 카피수 다형성**을 민감하게(정확하게) 잡아내지만, n=3/군에서는 이것이 "처치효과"와 "어느 개체가 어느 군에 배정됐는지"를 구분할 수 없게 만든다. mm39용 UCSC segmental-duplication 트랙이 없어 위치 기반 체계적 제외도 불가능했다.

**→ SV burden에 한해 single-sample(개별 콜) 결과를 채택: 유의한 유전자 없음(깨끗한 negative).**

---

## 2. 샘플단위 요약통계 (`cohort_summary_compare`) — 가장 견고한 결과

**개념 설명(VAF/hom-peak% 심화)**: [`VAF_het_hom_peak_frac_interpretation_2026-07-22.md`](./VAF_het_hom_peak_frac_interpretation_2026-07-22.md) — bulk VAF가 왜 "세포집단 구성비"를 반영하는지, het_hom_peak_frac이 왜 "우세 계통 존재"의 신호인지를 개념·수식·실제 수치로 정리.

**핵심 도구:** `bcftools`(query/view) + `scripts/sample_level_summary.py`(샘플별 지표 계산), `scripts/joint_genotype_discordance.py`(웰별 불일치율), `scripts/group_summary_compare.R`(Kruskal-Wallis + Friedman/Kendall's W)

**목적:** §1/§3/§4처럼 유전자마다 검정하면 다중검정 부담(n=3/군 대비 검정 수백~수천 개)으로 검정력이 사실상 0에 가까워지는 문제를 우회 — 유전자 단위가 아니라 **샘플/웰 단위로 지표 하나당 검정 1회**만 수행해서, 세포집단 이질성/클론성(hom-peak%, joint genotype 불일치율 등)이 처치강도(A<B<C<D 등)에 따라 일관되게 변하는지 확인.

**입력:**
- 지표 계산용: 12샘플의 `phased_small_variant_vcf`/`phased_sv_vcf`(single-sample)
- 불일치율 계산용: 4웰의 joint calling 출력(`/mnt/hdd1_1tb/h2o2-wgs/joint-H2O2-{A,B}01/*/out/split_joint_small_variant_vcfs/`, `/mnt/hdd2_1tb/h2o2-wgs/joint-H2O2-{C,D}01/*/out/split_joint_small_variant_vcfs/`)

**출력:**
- `h2o2_analysis_results/cohort/sample_metrics_combined.tsv` — 샘플별 원자료(12행)
- `h2o2_analysis_results/cohort/discordance_combined.tsv` — 웰별 원자료(4행)
- `h2o2_analysis_results/cohort/summary_compare_stats.csv` — 지표별 Kruskal-Wallis 결과
- `h2o2_analysis_results/cohort/summary_compare_stats_rank_matrix.csv` — 지표×군 순위행렬
- `h2o2_analysis_results/cohort/summary_compare_stats_friedman.csv` — Friedman/Kendall's W 결과
- `h2o2_analysis_results/cohort/summary_compare_plots.pdf` — 지표별 boxplot 2×2 패널(통계주석·평균값 라벨 포함)
- `h2o2_analysis_results/cohort/summary_compare_plots_discordance_rank.pdf` — 불일치율 barplot + 순위 히트맵 2패널

**출력 데이터 읽는 법:**

- `sample_metrics_combined.tsv` 컬럼: `sample, group, total_pass_variants`(PASS 소변이 수), `total_pass_sv`(PASS SV 수), `n_het_snv`(이형접합 PASS SNV 수), `het_hom_peak_frac`(이형접합인데 VAF>0.9인 비율 — 클론성 지표, 클수록 균질), `median_het_vaf`(이형접합 VAF 중앙값)
- `discordance_combined.tsv` 컬럼: `well, group, n_sites_compared`(3배치 모두 유효 GT인 비교대상 사이트 수), `n_discordant`(3배치 GT가 불일치한 사이트 수), `discordance_rate`(불일치 비율 — 클수록 이질적)
- `summary_compare_stats.csv` 컬럼: `metric, statistic(χ²), df, pvalue, fdr` — **fdr<0.05면 그 지표가 4군 간 유의하게 다름**
- `rank_matrix.csv`: 각 지표(행)마다 A/B/C/D의 순위(1=최소~4=최대) — `het_hom_peak_frac`은 코드에서 이미 "이질성" 방향으로 반전 저장됨(값 자체는 클론성 방향이므로 주의)
- `friedman.csv`: `kendalls_w`(0~1, 1=완전 일치), `pvalue_asymptotic`(지표 수 적어 근사 부정확, 참고용)

**결과 해석:**

| 지표 | H (χ²) | p-value | FDR | 비고 |
|---|---|---|---|---|
| `total_pass_variants` | 9.97 | 0.0188 | **0.046** | 유의 |
| `het_hom_peak_frac` | 9.51 | 0.0232 | **0.046** | 유의 |
| `total_pass_sv` | 4.74 | 0.192 | 0.239 | — |
| `median_het_vaf` | 4.22 | 0.239 | 0.239 | — |

(상염색체만 사용 — 성별교란 제외 확인 완료)

`het_hom_peak_frac`(단일샘플 VAF 분포 기반)과 `discordance_rate`(3배치 joint genotype 비교 기반)는 서로 완전히 다른 원리로 계산되는데도, 방향을 통일해서 보면 **군 순위가 정확히 일치**한다(둘 다 A=4위, B=2위, C=3위, D=1위 — "이질성 정도" 기준):

| 지표 | A | B | C | D |
|---|---|---|---|---|
| het_hom_peak_frac (반전) | 4 | 2 | 3 | 1 |
| discordance_rate | 4 | 2 | 3 | 1 |

Kendall's W = 1.0. 4개 항목의 두 무작위 순위가 우연히 완전히 일치할 확률은 1/4! = 1/24 ≈ 4.2%로, 근사 카이제곱 p-value(0.112, m=2라 근사 부정확)보다 이 사실 자체가 더 의미 있다. **A(처치없음)가 가장 이질적, D(최강처치 추정)가 가장 클론성** — 기존 애드혹 분석(2026-07-14 문서 §3)의 A<C<B<D 그라디언트를 정식 파이프라인으로 재현.

**참고**: 웰당 3개 replicate(예: H2O2-A01-ctrl1, H2O2-2-A01, H2O2-3-A01)는 동일 gDNA의 재시퀀싱이 아니라 **독립적인 생물학적 재현(triplicate)**임이 확인됨 — 즉 배양·처치부터 별도로 수행된 3개 웰이라, 아래 §2.3의 replicate간 변동은 시퀀싱 기술 잡음이 아니라 진짜 생물학적 재현성을 반영한다.

### 2.3 세 번째 독립 증거: 다중카피 유전자군 카피수 안정성 (HiFi 고유 depth 기반 카피수 활용)

**핵심 도구:** Sawfish(`sawfish_call`, 원본 `humanwgs_singleton` WDL 실행에 이미 포함 — 별도 실행 불필요) — depth 기반 전체 게놈 카피수 세그멘테이션

**목적:** §1에서 문제가 됐던 다중카피 유전자군(Vmn1r/Mup/Pramel/Ear 등) 로커스의 **실제 카피수가 웰 내 3개 독립 배양 replicate 사이에서 얼마나 안정적인지** 확인 — 안정적이면 germline(개체 고유) 다형성, replicate마다 들쭉날쭉하면 배양·처치 중 발생하는 진짜 클론 동역학(서브클론 비율 변화)을 시사. §2.1/2.2와 완전히 다른 원리(depth 기반 카피수 vs VAF/genotype)로 같은 "이질성 그라디언트" 가설을 검증하는 세 번째 방법.

**입력:** 12샘플의 `{run_dir}/out/sv_copynum_bedgraph/*.bedgraph` (Sawfish가 계산한 구간별 정수 카피수, GC bias 보정됨) — 예: `/mnt/JJ_dis_8tb/h2o2-wgs/H2O2-A01-ctrl1/20260517_163205_humanwgs_singleton/out/sv_copynum_bedgraph/H2O2-A01-ctrl1.GRCm39.structural_variants.copynum.bedgraph`

**출력:** 이 세션에서 1회성으로 계산(별도 저장 파일 없음, 문서에 결과만 기록). 재현 방법: 각 샘플의 `copynum.bedgraph`(컬럼: chrom, start, end, copy_number)에서 관심 로커스와 겹치는 구간들을 길이 가중 최빈값으로 요약하면 됨.

**출력 데이터 읽는 법:** `copynum.bedgraph`는 tab-delimited, 컬럼 `chrom, start, end, copy_number`(정수). 인접한 동일 카피수 구간은 이미 병합돼 있음. `*.copynum.summary.json`은 염색체별 카피수 분포(`bases_per_copy_number`)와 `most_common_copy_number`를 제공 — 전체 게놈 개요 파악에 유용하나 이번엔 특정 로커스만 봐서 bedgraph를 직접 사용.

**결과 해석:**

7개 다중카피 로커스에서 웰 내 3-replicate 카피수 range(최대-최소)를 계산하고, 로커스별로 4개 군의 range에 순위(불안정할수록 높은 순위)를 매겨 평균한 결과:

| 군 | 평균 순위(카피수 불안정성) |
|---|---|
| A | 2.93 (가장 불안정) |
| C | 2.64 |
| B | 2.36 |
| D | 2.07 (가장 안정) |

**A>C>B>D** — §2.2의 hom-peak%·joint discordance와 정확히 같은 순서. 서로 다른 3가지 방법(VAF 분포, joint genotype 일치도, depth 기반 카피수 안정성)이 모두 같은 이질성 그라디언트를 가리킨다.

**한계 (반드시 함께 고려):**
1. 로커스 7개 중 Vmn1r-a/b/c는 같은 chr7 슈퍼클러스터의 인접 서브영역이라 서로 통계적으로 독립적인 증거가 아님 — 정식 검정(p-value)은 계산하지 않음, 기술적(descriptive) 순위 비교로만 제시.
2. 반복서열/다중카피 영역은 원래 depth 기반 카피수 콜링 자체의 노이즈가 큰 구간이라, 관찰된 "불안정성"이 순수 생물학적 신호인지 콜링 노이즈인지 완전히 분리되지 않음.
3. 그래도 위 §2.1/2.2와 원리가 전혀 다른 세 번째 방법에서 같은 방향이 재현된 것은 이질성 그라디언트 가설을 뒷받침하는 추가 근거로 볼 수 있음.

---

## 3. 유전자셋(GO pathway) 집계 검정 (`cohort_geneset_burden_{sv,variant}`)

**핵심 도구:** MGI GO annotation(GAF) → `scripts/fetch_go_genesets.py`(유전자셋 추출) → `scripts/geneset_burden.R`(SV/변이 burden을 유전자셋 단위로 집계 후 Fisher/χ²+Kruskal-Wallis)

**목적:** §1/§4의 개별 유전자 단위 검정이 다중검정 부담으로 검정력이 낮은 문제를, **생물학적으로 관련된 유전자 묶음(산화스트레스 반응, DNA손상복구) 전체를 하나의 검정 단위**로 묶어 우회 — H2O2가 일으키는 산화스트레스/DNA손상이라는 가설과 직접 관련된 경로 수준에서 처치군간 차이가 있는지 확인 (검정 수를 유전자셋 개수(2개)로 줄여 다중검정 부담 최소화).

**입력:**
- `resources/mgi.gaf.gz` (Mouse Genome Informatics GO 주석, https://current.geneontology.org/annotations/mgi.gaf.gz 에서 다운로드)
- `h2o2_analysis_results/structural_variants/{sample}.annotated_sv.tsv` (SV 모드) 또는 `small_variants/{sample}.canonical_impacts.tsv` (variant 모드)

**출력:**
- `resources/genesets_oxidative_dna_repair.json` — 유전자셋 정의(`oxidative_stress_response` 303개, `dna_repair` 457개 유전자 심볼 목록)
- `h2o2_analysis_results/cohort/geneset_burden_sv.csv`, `geneset_burden_sv_plots.pdf`
- `h2o2_analysis_results/cohort/geneset_burden_variant.csv`, `geneset_burden_variant_plots.pdf`

**출력 데이터 읽는 법 (`geneset_burden_*.csv` 컬럼):**

| 컬럼 | 의미 |
|---|---|
| `geneset` | 유전자셋 이름 |
| `n_genes` | 해당 유전자셋의 후보 유전자 수 |
| `pvalue_presence` | "샘플이 이 유전자셋 내 유전자 중 하나라도 이벤트를 가지는가"(이진) 기준 검정 p-value |
| `pvalue_count` | 샘플당 유전자셋 내 총 이벤트 수(연속값) 기준 Kruskal-Wallis p-value — 검정력 더 좋음 |
| `max_prop_diff`, `cramers_v` | §1과 동일한 효과크기 지표 |
| `n_presence_{A..D}` | 해당 군에서 유전자셋 내 이벤트를 가진 샘플 수 |
| `mean_count_{A..D}` | 해당 군의 샘플당 평균 이벤트 수 |
| `fdr_presence`, `fdr_count` | BH 보정 p-value (유전자셋이 2개뿐이라 보정 영향 작음) |

**결과 해석:**
- **SV 기반**: 둘 다 유의하지 않음 (`dna_repair` p_count=0.108, `oxidative_stress_response` p_presence=0.48/p_count=0.31)
- **VEP(소변이) 기반**: 12샘플 전체에서 HIGH/MODERATE 변이가 걸린 유전자가 80개뿐이고 이 중 유전자셋과 **겹치는 게 0개** (직접 확인 결과 버그 아님, 우연히 안 겹친 것 — 그래서 `pvalue_*`가 전부 NA)
- 주의(한계): 유전자셋은 GO term 직접 주석만 사용(하위 term 계층 전개 없음) — 완전한 GO closure가 아닌 근사치

---

## 4. VEP variant burden (`cohort_variant_burden`) — 깨끗한 negative

**핵심 도구:** Ensembl VEP(`docker://ensemblorg/ensembl-vep:release_110.1`, singularity) + mouse cache(release 110, GRCm39) → `bcftools +split-vep` → `scripts/variant_burden.R`

**목적:** 소변이(SNV/indel) 중 단백질 기능에 영향을 줄 가능성이 큰(HIGH/MODERATE impact) 변이만 골라, 그 유전자별 부담이 4개 처치군 사이에 다른지 검정 — §1의 SV 대신 점변이 관점에서 같은 질문(처치강도에 비례한 유전자 손상 부담 차이)을 확인.

**입력:** 12샘플의 `phased_small_variant_vcf` → `filter_small_variants`(slivar 필터) → `annotate_vep` → `extract_vep_canonical`

**출력:**
- 중간산출물: `h2o2_analysis_results/small_variants/{sample}.canonical_impacts.tsv` (CHROM, POS, REF, ALT, SYMBOL, IMPACT, Consequence, GT)
- `h2o2_analysis_results/cohort/variant_burden_stats.csv`, `variant_burden_plots.pdf` (6페이지: raw-p 순위+효과크기, volcano, bubble, heatmap)

**출력 데이터 읽는 법:** 컬럼은 §1 `sv_burden_stats.csv`와 동일 구조(gene, pvalue, odds_ratio, max_prop_diff, group_max/min, odds_ratio_extreme, cramers_v, n_A..n_D, fdr) — 대상이 SV가 아니라 HIGH/MODERATE impact 소변이라는 점만 다름.

**결과 해석:** 80개 유전자(HIGH/MODERATE impact) 테스트, **FDR<0.05 유의 유전자 0개**. 최상위 B4galt3도 raw p=0.043, FDR 보정 후 유의성 없음. 성염색체 편중은 SV burden만큼 심하지 않아 별도 재확인은 하지 않았음(주의: 완전히 배제 검증된 것은 아님).

---

## 5. DMR pairwise (`dmr_pairwise`) — 상위 후보가 성별 교란이었음, 상염색체 재확인 필요했음

**핵심 도구:** `DSS::DMLtest`+`DSS::callDMR`(R/Bioconductor, `docker://bioconductor/bioconductor_docker:RELEASE_3_18` 환경) via `scripts/dmr_analysis.R`

**목적:** 유전(변이) 수준이 아니라 **후성유전학(CpG 메틸화)** 수준에서 군간 차이가 나는 영역(differentially methylated region, DMR)을 찾아 — H2O2 처치가 유전체 서열 변화 없이도 메틸화 패턴을 바꾸는지, 그리고 §2에서 발견한 "이질성 그라디언트"가 메틸화 패턴에도 나타나는지 확인.

**입력:** 4군×3배치 `cpg_combined_bed`(CpG 메틸화 pileup) — `h2o2_analysis_results/cohort/methylation_{A,B,C,D}.txt`(그룹별 샘플 목록)에서 참조

**출력:** 6개 군쌍(`A_vs_B, A_vs_C, A_vs_D, B_vs_C, B_vs_D, C_vs_D`) 각각:
- `h2o2_analysis_results/cohort/dmr_{group1}_vs_{group2}.csv`
- `h2o2_analysis_results/cohort/dmr_{group1}_vs_{group2}_plots.pdf` (2페이지: 메틸화차 분포, top DMR 등)

**출력 데이터 읽는 법 (`dmr_*.csv` 컬럼):**

| 컬럼 | 의미 |
|---|---|
| `chr, start, end, length` | DMR 후보 영역 위치 |
| `nCG` | 영역 내 CpG 사이트 수 |
| `meanMethy1`, `meanMethy2` | group1/group2의 평균 메틸화 비율(0~1) |
| `diff_methy` | meanMethy1 - meanMethy2 (효과크기, +면 group1이 더 높은 메틸화) |
| `areaStat` | 영역 내 개별 CpG Wald 통계량의 합 — **유일하게 유효한 순위 지표**, `\|areaStat\|` 내림차순 정렬됨 |
| `pvalue`, `fdr` | **항상 NA** — DSS::callDMR()이 영역 단위 p-value를 제공하지 않음(아래 §5.1) |

**결과 해석:**

**§5.1 pvalue/fdr이 NA인 이유(설계상 정상, 버그 아님):** `DSS::callDMR()`은 영역 단위 p-value를 제공하지 않는다. areaStat/√nCG를 표준정규분포 근사로 처리하는 결합검정(Stouffer's method)을 시도했으나, 인접 CpG의 공간적 상관관계(smoothing으로 더 심함)를 무시해 반코저버티브였다(전체 영역의 99%+ 가 FDR<0.05로 나오는 등 비현실적) — 폐기하고 정직하게 NA로 복원, 대신 `areaStat` 기준 정렬로 대체(`scripts/dmr_analysis.R` 주석 참고). 엄밀한 검정을 하려면 permutation test(그룹 라벨을 섞어 areaStat의 귀무분포를 영역 길이/nCG별로 추정)가 필요 — 별도 작업으로 남음.

**§5.2 성별 교란 발견 및 상염색체 재확인:** 전체 영역 기준 6개 군쌍 전부 최상위 후보가 chrX(A_vs_B/C/D, C_vs_D) 또는 chrY(B_vs_C)였다 — 처치효과가 아니라 **성별에 따른 X-inactivation/도스보상 등 성염색체 고유의 생물학적 메틸화 차이**로 봐야 한다. 상염색체만으로(csv를 `chr` 컬럼 기준 사후 필터링) 재확인한 최상위 후보:

| 비교 | 위치 | diff_methy | areaStat |
|---|---|---|---|
| A_vs_B | chr11:104,705,502-104,706,036 | +0.213 | 1300.3 |
| A_vs_C | chr18:7,608,927-7,609,733 | -0.230 | -1393.7 |
| A_vs_D | chr4:59,034,926-59,035,774 | -0.136 | -1110.8 |
| B_vs_C | chr11:11,975,832-11,976,648 | -0.135 | -1334.6 |
| B_vs_D | chr6:28,927,796-28,929,707 | -0.233 | -1001.6 |
| C_vs_D | chr6:28,927,796-28,929,822 | -0.308 | -1446.4 |

**주목**: B_vs_D와 C_vs_D가 **같은 chr6:28.9Mb 위치**를 최상위로 가리킴 — 둘 다 D가 관여하는 비교라 D 특이적 상염색체 메틸화 후보일 가능성. 다만 areaStat 기반 순위일 뿐 통계적 유의성 검증은 안 됨 — 확정 아님, IGV 확인 등 후속 검증 필요.

---

## 6. 종합 결론

1. **가장 신뢰할 수 있는 발견은 §2 (샘플단위 요약통계)**: 서로 다른 원리의 세 지표(hom-peak%, joint discordance, 다중카피 유전자군 카피수 안정성)가 모두 일치하는 A<C<B<D 이질성 그라디언트. 성별 교란 배제 확인 완료, replicate는 독립 생물학적 재현(triplicate)으로 확인됨.
2. **SV burden은 single-sample 결과가 유효, joint는 다중카피 유전자군에 오염돼 신뢰 불가.**
3. **VEP variant burden, 유전자셋 검정은 깨끗한 negative** (유의 신호 없음).
4. **DMR은 성별 교란을 반드시 제거하고 봐야 함** — chr6:28.9Mb(D 특이적 추정) 등 상염색체 후보는 있으나 미확정.
5. **성별이 실험군과 완전히 교란돼 있다는 사실 자체가 이 데이터셋의 가장 중요한 구조적 한계**다. 향후 유사 실험은 반드시 웰별 성별을 사전에 기록·균형 배치해야 한다.
6. **HiFi 롱리드 고유 기능(phasing, ASM, read-level HP 태그) 3가지를 추가 검증**한 결과도 전부 A<C<B<D 방향과 일치 — 그중 haplotype 방향성 allelic imbalance만 완전히 독립적인 새 원리이고, 나머지 둘은 이형접합 밀도 감소라는 이미 아는 원인의 재확인. 상세: [`wgs-tertiary-analysis/docs/H2O2_hifi_native_analysis_plan_2026-07-16.md`](../shared/wgs-tertiary-analysis/docs/H2O2_hifi_native_analysis_plan_2026-07-16.md).
7. **LOH 쏠림 지점의 유전자를 4군 교차검증**해 Hsf2bp/Tusc3(약물 보호 실패 조건 B+D에서만 공통), Piwil1(D 특이적) 등 탐색적 후보 발굴 — 확정 아님, 상세: [`H2O2_LOH_hotspot_cross_validation_2026-07-22.md`](./H2O2_LOH_hotspot_cross_validation_2026-07-22.md).

## 7. 관련 파일 (모두 `/data_4tb/shared/wgs-tertiary-analysis/` 기준 상대경로)

- `config/config.yaml`, `config/config_joint.yaml` — 실행 설정(`exclude_chroms`, `exclude_gene_prefixes` 포함)
- `scripts/sv_burden.R`, `scripts/variant_burden.R`, `scripts/dmr_analysis.R`, `scripts/sample_level_summary.py`, `scripts/joint_genotype_discordance.py`, `scripts/group_summary_compare.R`, `scripts/geneset_burden.R`, `scripts/fetch_go_genesets.py`
- `h2o2_analysis_results/`, `h2o2_analysis_results_joint/` — 실행 산출물 (각 절의 "출력" 참고)
- `docs/H2O2_group_comparison_v2_plan_2026-07-15.md` — 1~3단계 실행 계획
