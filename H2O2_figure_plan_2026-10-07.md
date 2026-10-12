# H2O2 astrocyte HiFi WGS — Figure 구성안 (캡션 수준)

**작성일:** 2026-10-07 · **개정:** 2026-10-12 (품질 교란 검증 결과 반영, 스토리라인 전면 개정)
**용도:** publication 수준 figure 제작 파이프라인에 넘길 패널별 캡션·데이터 출처·상태 정의
**스토리라인:** 처치는 검출 가능한 유전체 변이를 늘리지 않음(Fig 2) → 클론 구성도 바뀌지 않음(Fig 3) → 처음 보였던 조건 간 그라디언트는 라이브러리·시료 수준의 기술적 교란(Fig 4) → 메틸화도 전역·국소 모두 처치 효과 없음, 대조군 대 처치군 비교는 설계상 분리 불가(Fig 5) → 5hmC·6mA(Fig 6) → 종합과 설계 교훈(Fig 7)
**경로 기준:** `TA=/data_4tb/shared/wgs-tertiary-analysis/h2o2_analysis_results`, `BM=/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample`, `FW=$BM/fig_atlas_work`, `QM=$FW/qmatch`

패널 상태 표기: **[ready]** 데이터 있음, 그림만 다시 그리면 됨 · **[analysis]** 추가 분석 필요(짧음) · **[pending]** 실행 중인 계산 결과 대기 · **[decide]** 결과 보고 넣을지 결정 · **[withdrawn]** 이전 판의 주장 철회(근거는 Fig 4)

---

## 개정 요약 (2026-10-12)

| 이전 판의 주장 | 판정 | 근거 |
|---|---|---|
| Fig 3: 처치에 따른 클론 균질화(A<C<B<D) | **철회** | 지표들이 indel을 포함했고 STR 수축 artifact·호출 민감도 차이를 따라감. SNV만 쓴 클론성 지표와 read 수준 VAF 비교에서는 조건 차이 없음 |
| Fig 2e: H2O2에 따른 C>A(8-oxoG) 증가 | **철회** | 무처치(A)의 비-C>A 오류 호출이 많아 생긴 비율 효과. 품질 맞춘 재호출에서 C>A 비율 p = 0.30 |
| Fig 4c·4d·4f: 무처치 반복 샘플의 메틸화 이질성, 유전 이질성과의 상관 | **철회/Supplementary** | 같은 교란. DMR 상위 영역도 라이브러리 순서를 따름 |
| Fig 2b–d, 2f: 변이 부담 음성 + 검출력 | 유지 | 교란이 A를 더 "변이가 많아 보이게" 했으므로 "처치로 증가 없음"은 더 단단해짐 |
| Fig 4a: 전역 5mC 동일 | 유지 | |
| 6mA: 생물학적 신호 없음 | 유지(Supplementary) | 5-mer 문맥 PC1이 read QV를 따름(ρ = 0.86) |

---

## 0. 확정된 설계 정보

- **그룹:** A = 무처치(control), B = 급성 H2O2, C = 급성 H2O2 + 약물, D = 만성 H2O2 + 약물. 2×2 요인 설계가 아님.
- **비교 대비(contrast):** H2O2 = B − A, 약물(급성 H2O2 하) = C − B, 만성 vs 급성(약물 하) = D − C, 그리고 D − A. 통계 모델은 `condition + run`, 유의성은 run 안 정확 permutation(13,824 배치, 최소 p = 0.0017) 또는 무작위 permutation.
- **배양:** 여러 개체가 섞인 primary culture. 공유 이형접합 SNV의 약 60 %가 Mouse Genomes Project(REL2021) 계통 SNP와 일치 → **배양 개체가 순수 C57BL/6가 아니거나 비-B6 배경이 섞여 있음**(실험팀 확인 필요).
- **구조적 교란:** 조건 A는 세 run 모두 SMRT cell 1번 슬롯이고 read quality가 가장 낮음(A Q29–32, B/C Q33–35, D Q34–36). 라이브러리 수준 artifact 지표(STR 수축 지수)도 A > C > B > D 순서. **"A 대 처치군(B·C·D)" 차이는 H2O2 공통 효과와 A 시료·라이브러리 특성을 이 설계로 분리할 수 없음.** 처치군끼리의 대비(C − B, D − C)가 상대적으로 깨끗한 비교.
- **실험팀 확인 요청 사항:** 라이브러리 제작 방식(PCR-free 여부, 증폭 cycle), 샘플별 DNA 투입량, DIN/GQN, 추출 날짜·배치, 마우스 계통, 네 웰이 하나의 혼합 세포 준비물에서 나뉘었는지.

---

## 공통 스타일 규칙

- **그룹 색:** 4색 고정 팔레트(A는 회색 계열 control, B–D는 처치 강도 순으로 진해지는 단일 hue 계열 또는 colorblind-safe 4색). 모든 figure에서 동일.
- **기술적 지표 색:** Fig 4·5에서 라이브러리 지표(STR 수축 지수, read QV)를 표시할 때는 조건 색과 다른 연속 색척(예: viridis)을 써서 "조건"과 "기술"을 시각적으로 구분.
- **run(배치) 표시:** 점 모양(●▲■)으로 run 1/2/3 구분. 같은 run의 샘플은 얇은 선으로 연결해 run 내 순위 재현성을 보여줌.
- **n이 작으므로:** 모든 boxplot은 개별 점(n=3)을 함께 표시. 막대그래프만 단독 사용 금지.
- **통계 표기:** 검정 이름·p·n을 캡션에 명시. permutation 결과는 null 분포 histogram + 관측값 세로선으로.
- **음성 결과 표기:** 효과 없음은 "p > 0.05"가 아니라 검출 가능한 최소 효과 크기(또는 신뢰구간)로 표현.
- **염색체:** 군간 비교는 상염색체만(성 교란 때문). 성염색체는 Fig 1d/Supp에서만.
- **출력:** 벡터(SVG/PDF), 패널 문자 a, b, c… 소문자 bold, 폰트 Arial/Helvetica 6–8 pt.

---

## Figure 1. Experimental design and a single long-read dataset that reads genome and epigenome together

**a. Study design (schematic).** Primary astrocytes pooled from multiple mice were cultured under four conditions — A, untreated; B, acute H2O2; C, acute H2O2 + drug; D, chronic H2O2 + drug —, each in three independent biological replicates processed in separate sequencing runs (run 1–3; n = 12). Genomic DNA from each well was sequenced on PacBio Revio (HiFi). Schematic also shows the SMRT cell slot of each sample (A always slot 1). **[ready — 모식도]**

**b. One dataset, multiple layers (schematic).** From the same HiFi reads: small variants (DeepVariant), structural variants and copy number (Sawfish), haplotype phasing (HiPhase), and base modifications from polymerase kinetics (5mC, 5hmC, 6mA; jasmine/pb-CpG-tools/modkit). Arrows indicate which layer feeds which figure. **[ready — 모식도]**

**c. Sequencing quality per sample.** HiFi read count (5.6–9.3 M), read-length N50 (14.3–16.9 kb), median read quality (Q29–Q36), mapping rate (≥ 99.9 %), mean autosomal depth (25–42×), and fraction of heterozygous variants phased (0.61–0.81) for the 12 samples (colour, condition; shape, run). Untreated samples (A) have the lowest read quality in every run (Q29–32 vs Q33–36), a technical difference that Fig. 4 shows to drive the apparent between-condition gradients. **[ready — `$BM/figure_data/sequencing_qc.tsv`]**
  - phase_block_ng50 값은 stats_file에서 대부분 0으로 기록돼 쓰지 않음.

**d. Composition of the starting cell populations.** Per-sample chrX/autosome and chrY/autosome depth ratios: A and D cultures contain male-derived cells (chrY/autosome ≈ 0.13), whereas B and C are female-only. About 60 % of heterozygous SNVs shared between samples are known Mouse Genomes Project strain polymorphisms, indicating non-C57BL/6 genetic background in the pooled animals. **[ready — mosdepth summary; MGP 비율 `$QM/snv_origin/mgp.log`]** *(성 구성 차이를 선택으로 쓸지 출발 집단 차이로 쓸지는 웰 준비 방식 답에 따라 결정.)*

---

## Figure 2. Treatment does not increase detectable mutational burden

**a. Variant counts per sample.** Autosomal PASS SNVs, indels, and SVs per sample. Between-condition differences in small-variant counts are technical (Fig. 4); SV counts do not differ. **[ready — `sample_metrics_combined.tsv`, 상염색체 재집계 확인]**

**b. Gene-level burden of protein-altering small variants.** Rank plot (or volcano) of HIGH/MODERATE-impact SNV/indel burden across conditions for 80 genes; no gene reaches FDR < 0.05 (top hit nominal p = 0.04, FDR = 1). **[ready — `TA/cohort/variant_burden_stats.csv`; fig-atlas 번들 fig2b 준비됨]**

**c. Gene-level burden of gene-disrupting SVs.** Same as b for SVs from single-sample calls (278 genes, 0 at FDR < 0.05). **[ready — `TA/cohort/sv_burden_stats.csv`; 번들 fig2c]**

**d. Pathway-level burden.** Forest plot of 20 curated gene sets (oxidative stress, BER/NER/MMR/HR/NHEJ, apoptosis, ferroptosis, mitochondrion, …) for small variants and SVs; minimum FDR 0.33–0.39. **[ready — `geneset_burden_{variant,sv}_curated.csv`; 번들 fig2d(상염색체 재계산)]**

**e. Detection power.** With the observed between-replicate variation (residual SD of log private-SNV burden 0.18; df = 6), a 3-vs-3 run-paired comparison has 80 % power at α = 0.05 to detect a 1.65-fold change in total burden. No treated condition exceeds the untreated burden; the only between-condition differences run in the opposite direction (untreated highest) and are technical (Fig. 4). **[ready]**

**f. Spectrum of sample-private SNVs after quality matching (Supplementary 후보).** Six-class spectrum of private SNVs from quality- and depth-matched re-calls (chr1–3, 22.8×): C>A fraction does not differ between conditions (A 8.6 %, B 8.3 %, C 9.8 %, D 11.1 %; condition + run F = 1.68, exact permutation p = 0.30), nor does C>A per heterozygous site (p = 0.99). The previously reported C>A rise in the original calls was a ratio effect of excess non-C>A error calls in the low-quality untreated libraries. **[ready — `$QM/qm_compare.summary.txt`, `$QM/qm_compare.samples.tsv`]**
  - 이전 판 2e(C>A 증가, SBS18 유사도)와 2g(전사 가닥 비대칭)는 **[withdrawn]**. 전사 가닥 분석에서 C>A가 대칭이었던 것은 이 결론과 일관됨 — Supplementary 각주로만.

*Supplementary:* joint-called SV burden giving 112 "significant" genes that are all multi-copy gene families (alignment artefact) — 왜 single-sample을 채택했는지.

---

## Figure 3. Clonal composition of the cultures is unchanged by treatment

**a. Principle (schematic).** In a polyclonal bulk population, expansion of one lineage shifts bulk allele fractions: alleles of the expanding lineage rise and those of receding lineages fall. Two read-outs: (i) the fraction of heterozygous SNVs drifting to VAF > 0.9 (loss of heterozygosity in a dominant lineage) and (ii) the spread of per-site VAF changes between two samples, which widens in both directions under a composition change. **[ready — 모식도]**

**b. Homozygous-peak fraction from SNVs only.** Fraction of heterozygous SNVs (indels excluded) with VAF > 0.9 per sample: no condition order in the original calls (condition + run F = 0.46, exact permutation p = 0.67) or in quality-matched re-calls (F = 0.10, p = 0.96). **[ready — `$QM/qm_compare.samples.tsv` (chr1–3), 전 게놈 값 `$FW/vcfpass/`]**

**c. Haplotype-directional allelic imbalance from SNVs only.** Fraction of phase blocks with consistent haplotype skew, SNVs only: 0.65–0.75 in all samples with no condition order (A 0.69–0.75, B 0.67–0.70, C 0.69–0.73, D 0.65–0.75). Shown next to the previously reported all-variant version (A 0.09–0.14 → D 0.67–0.74), whose gradient disappears once indels are removed. **[ready — `$FW/robust/*.robust.tsv`]**

**d. Pairwise read-level VAF shifts (core panel).** For 21,813 heterozygous SNV sites, standardized VAF differences between every pair of samples (z = ΔVAF / binomial SE; read counts from quality-matched BAMs). Variance of z for untreated–chronic (A–D) pairs (1.30, 1.43, 1.50) equals that of same-condition pairs from different runs (mean 1.32; e.g. A1–A2 1.43, D1–D2 1.19), with no directional shift (mean z ≈ 0) and equal fractions of sites shifting up or down by > 3 SE (0.6–0.7 % each). Non-repeat sites give the same result. Given ~30× depth, composition changes exceeding ~0.03–0.04 VAF beyond well-to-well variation would have been detected. **[ready — `$QM/vaf_shift/vs.pairs.tsv`, `vs.sites.tsv.gz`]**
  - 그림 형태: 왼쪽 — 조건 쌍 유형별(같은 조건·다른 run / 같은 run·다른 조건) var_z 점 그림; 오른쪽 — A1–D1과 A1–A2의 z 분포 겹침 histogram.

**e. Replicate genotype discordance after quality matching (Supplementary 후보).** SNV pairwise discordance among replicates: A 0.171, B 0.165, C 0.168, D 0.141 (original calls 0.196 / 0.173 / 0.182 / 0.146). The untreated excess largely disappears with quality matching. **[ready — `$QM/qm_compare.pairs.tsv`]**

*이전 판 3c(hom-peak, indel 포함), 3d(불일치율), 3e(전체 변이 haplotype 불균형), 3g(방법 간 일치), 3h(LOH hotspot) — **[withdrawn]**. 3f(copy-number instability)는 기술 교란과 독립성이 확인되지 않아 제외.*

---

## Figure 4. The apparent condition gradient is a library-level technical artefact

**a. Microsatellite contraction.** Per sample, 2-bp deletions versus 2-bp insertions at dinucleotide repeats of 10–30 units (chr1): untreated samples carry 10,428–13,194 deletions but only 201–274 insertions (contraction index 0.98), chronic H2O2 + drug samples are nearly balanced (670–1,063 vs 503–571; 0.57–0.66). The index orders samples D < B < C < A in all three runs and correlates with read quality (Spearman ρ = −0.86; 1-bp: −0.95). Median VAF of the excess deletions ~0.38, i.e. present in many molecules — a library/DNA-level property (polymerase slippage-like), not random read error. **[ready — `$QM/str_contraction_chr1.tsv`, `$QM/deletion_character_run1.{summary,dist}.tsv`]**

**b. Quality- and depth-matched re-calling.** Reads were subsampled per sample to an identical read-quality distribution and depth (22.8×, chr1–3) and re-called with the original DeepVariant/GLnexus settings. Heterozygous deletions halve but keep the D < B < C < A order (p = 0.0017); low-VAF SNV fraction (p = 0.014 → 0.44), private-SNV count (p = 0.0035 → 0.064) and replicate discordance lose their gradient; heterozygous SNV count keeps it (p = 0.0017). Read quality per se therefore explains the low-VAF/private differences, while the deletion and SNV-count differences persist. **[ready — `$QM/qm_compare.summary.txt`, `qm_compare.samples.tsv`; 방법 `fig_atlas_export/quality_matched_subsample.py`, `run_qm_deepvariant.sh`]**
  - 그림 형태: 지표별로 original vs matched를 나란히(run 연결선), 패널 위에 permutation p.

**c. The residual SNV excess is a calling-sensitivity difference.** At sites with equal read-level VAF, DeepVariant calls a heterozygous genotype more often in untreated libraries: VAF 0.25–0.35, A 0.53–0.59 vs D 0.30–0.42; VAF 0.35–0.45, A 0.58–0.71 vs D 0.36–0.47, while the number of sites per VAF bin is the same across conditions. **[ready — `$QM/vaf_shift/call_rate_by_vaf.tsv`]**

**d. The extra SNVs are real alleles present in all cultures.** SNVs called only in A (vs D of the same run; 3,981–4,837 per run) are mostly supported by reads in D (69 % with ≥ 5 alt reads; median VAF in D 0.27), recur in other A replicates (55–65 %), match known strain SNPs at the same rate as shared heterozygous sites (51–57 % vs 60 % baseline), and have the same substitution spectrum as the reverse set. **[ready — `$QM/snv_origin/mgp.{log,sites.tsv.gz,summary.tsv}`]**

**e. Kinetic signature follows the same axis (Supplementary 후보).** PC1 (77 % of variance) of per-sample 6mA 5-mer context profiles orders samples A < C < B < D in every run and correlates with read quality (ρ = 0.86). **[ready — `$BM/6mA/stats/`]**

> Fig 4 메시지: 한 축의 라이브러리 차이(A > C > B > D)가 STR 수축, 호출 민감도, kinetics 문맥을 함께 움직여 "클론성"·"돌연변이 스펙트럼" 신호를 만들어냈다. 원인 특정(증폭 여부, 투입량, DNA 품질)은 실험팀 정보 대기.

---

## Figure 5. DNA methylation (5mC) shows no treatment effect

**a. Global CpG methylation.** Mean 5mC across 18.8 M autosomal CpGs covered ≥10× in all samples; 0.760–0.774, no condition effect. **[ready — 번들 fig4b 준비됨]**

**b. Region-level differential methylation, permutation-calibrated.** limma-trend on region log-odds (condition + run) for promoters (49,515), gene bodies (20,335), CpG islands (15,390), CGI shores (30,777) and 1-kb tiles (1.66 M). For each contrast, number of regions with p < 10⁻³ (point) against its within-run permutation null (1,000 permutations; 200 for tiles; violin). No region reaches FDR < 0.1 except one tile (D − A) and one gene body (drug). Treated-vs-treated contrasts (drug C − B, chronic D − C) are at the null in every region set (permutation p 0.21–0.97). **[ready — `$FW/dmr/limma/{set}.summary.tsv`, `{set}.regions.tsv.gz`; 방법 `fig_atlas_export/dmr_region_counts.py`, `dmr_limma.R`]**

**c. Untreated-versus-treated contrasts show a modest excess that tracks the library axis.** Contrasts involving A show ~1.5-fold more regions at p < 10⁻³ than the null (tiles: B − A 2,036 vs null median 1,280, p = 0.025; D − A 2,242 vs 1,227, p = 0.010; gene bodies B − A p = 0.018; shores p = 0.033; promoters p = 0.055; CGIs p = 0.29). Adding the STR-contraction index as a covariate returns tile counts to the null (B − A 1,345; D − A 1,285). The strongest D − A tiles follow the library order (e.g. chr13:50.13 Mb: A 0.19 < C 0.48 < B 0.60 < D 0.76; |r| with the index 0.80–0.94 for the top 10). Because A is both the only untreated condition and the condition with the distinct library (and SMRT cell slot 1), an untreated-versus-treated methylation difference cannot be separated from the technical difference in this design. **[ready — 같은 파일; 대표 영역 표는 `tile1k.regions.tsv.gz`의 p_DvsA 상위]**
  - 그림 형태: 대비 × 영역 집합 격자(관측/ null 비율, 공변량 보정 전·후), 오른쪽에 대표 타일 1개를 run 연결선 점 그림으로.

*Supplementary:* PCA(이전 4b, A 반복 샘플이 흩어짐)와 반복 샘플 간 메틸화 거리(이전 4d) — 같은 기술 축으로 설명될 수 있음을 명시하고 Supplementary로만. 이전 4c(CpG 단위 condition 효과), 4f(메틸화–클론성 상관), ASM 개수 — **[withdrawn]**.

---

## Figure 6. 5hmC and 6mA from the same reads

**a. Global 5hmC at CpGs per sample.** First sample (2-A01): 21.4 M CpGs, mean coverage 21.7×, global 5mC 70.7 %, 5hmC 1.1 %. **[pending — jasmine 재호출 진행 중, 전체 완료 약 10월 하순]**

**b. 5hmC by genomic element** (promoter, gene body of expressed vs silent genes, CpG island), with validation that 5hmC is enriched in gene bodies of active genes. Biological rationale: TET-mediated 5mC oxidation is sensitive to oxidative stress. **[pending]**
  - 해석 기준(미리 정함): 차이가 A < C < B < D(라이브러리 순서)를 따르면 artifact로 간주; 라이브러리 지표를 공변량으로 넣은 모형 병행; 처치군끼리의 대비(C − B, D − C)를 주 결과로.

**c. 5mC–5hmC relationship** at gene bodies per condition. **[pending]** **[decide]**

**d. Technical characterization of 6mA calls (Supplementary).** Per-sample 6mA call rate (~0.6 % of adenines in nuclear reads, far above credible levels of native mammalian 6mA), its 5-mer context dependence, chrM versus nuclear rate, and absence of mycoplasma (0 reads of 7 species in 1.2 M reads). Gene- and promoter-level 6mA show no condition effect beyond permutation; context-profile PC1 tracks read quality (Fig. 4e). **[ready — `$BM/6mA/stats/`, `$BM/figure_data/mycoplasma_summary.tsv`]**

> Fig 6 구성은 5hmC 결과에 따라 결정: 처치군 간 신호가 있으면 Fig 6 = 5hmC 단독, 없으면 5hmC·6mA 모두 Supplementary로 내리고 Fig 7을 Fig 6으로.

---

## Figure 7. Summary and design lessons

**Summary schematic.** Oxidative stress (acute or chronic H2O2, with or without drug) did not leave a detectable mark on the genome of cultured astrocytes — no increase in point mutations, indels or SVs (Fig. 2), no change in clonal composition (Fig. 3), and no global or focal 5mC change between treated conditions (Fig. 5). The between-condition gradients first observed were produced by a library-level difference concentrated in the untreated samples (Fig. 4). Design lessons for follow-up: randomize conditions across SMRT cell slots, prepare libraries in one batch with matched input and DNA quality, and record library QC (DIN/GQN, amplification) per sample. **[모식도]**

---

## 추가 분석·작업 목록

| 항목 | 작업 | 상태 |
|---|---|---|
| 전 게놈 확장 | MGP·재현성·read 근거 검사를 원래 전 게놈 호출에 적용(chr1–3 결과 확인용) | 미시작 — MGP 전 게놈 다운로드 약 10시간 |
| 6a–c | 5hmC 12샘플 완료 후 분석 | pending |
| 실험팀 정보 | 라이브러리 제작·DNA QC·마우스 계통·웰 준비 방식 | 사용자 문의 중 |
| fig-atlas 요청서 | 철회 패널 반영해 Delivery Tracker note 갱신(전달 전 상태 변경 없음) | 사용자 확인 후 |
