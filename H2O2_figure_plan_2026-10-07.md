# H2O2 astrocyte HiFi WGS — Figure 구성안 (캡션 수준)

**작성일:** 2026-10-07
**용도:** publication 수준 figure 제작 파이프라인에 넘길 패널별 캡션·데이터 출처·상태 정의
**스토리라인:** 돌연변이 증가 없음(Fig 2) → 대신 세포 집단 구조 변화(클론 균질화, Fig 3) → 후성유전체도 같은 방향(Fig 4–5) → 종합 모델(Fig 6)
**경로 기준:** `TA=/data_4tb/shared/wgs-tertiary-analysis/h2o2_analysis_results`, `BM=/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample`

패널 상태 표기: **[ready]** 데이터 있음, 그림만 다시 그리면 됨 · **[analysis]** 추가 분석 필요(짧음) · **[pending]** 실행 중인 계산 결과 대기 · **[decide]** 결과 보고 넣을지 결정

---

## 0. 확정된 설계 정보 (2026-10-07 사용자 확인)

- **그룹:** A = 무처치(control), B = 급성 H2O2, C = 급성 H2O2 + 약물, D = 만성 H2O2 + 약물. 2×2 요인 설계가 아님.
- **비교 대비(contrast):** H2O2 효과 = B − A, 급성 H2O2에서 약물 효과 = C − B, 약물 존재 하 만성 vs 급성 = D − C. 통계 모델은 `condition + run`.
- **배양:** 여러 개체가 섞인 primary culture. 상염색체 이형접합 SNV가 샘플당 9만–38만 개인 것은 개체 간 다형성 때문이며, 이 이형접합 자리들이 Fig 3 클론성 지표의 재료가 된다.
- **남은 확인 사항:** 네 웰이 하나의 혼합 세포 준비물을 나눈 것인지, 웰마다 다른 개체 pool인지. A/D에만 수컷 세포가 있고 B/C에는 없다는 점(Fig 1d)을 어떻게 서술할지가 이 답에 달려 있다.
  - 같은 준비물이면: B/C(급성 H2O2)에서 수컷 유래 세포가 사라졌다는 것 자체가 강한 선택의 증거가 된다.
  - 다른 pool이면: 출발 집단 차이로 서술하고, Fig 3은 "웰 간 출발 구성 차이를 감안해도 처치 순서와 일치"로 표현 수위를 낮춘다.

**스토리라인에서 달라지는 점:** 균질화 순서 A < C < B < D는 이제 다음처럼 읽힌다. 급성 H2O2가 집단을 균질화하고(B > A), 약물이 이를 부분적으로 막으며(C < B), 약물이 있어도 만성 노출은 가장 강한 균질화를 남긴다(D > B). 즉 "약물의 부분적 rescue"가 클론 구조 수준에서 보인다.

---

## 공통 스타일 규칙

- **그룹 색:** 4색 고정 팔레트(A는 회색 계열 control, B–D는 처치 강도 순으로 진해지는 단일 hue 계열 또는 colorblind-safe 4색). 모든 figure에서 동일.
- **run(배치) 표시:** 점 모양(●▲■)으로 run 1/2/3 구분. 같은 run의 샘플은 얇은 선으로 연결해 run 내 순위 재현성을 보여줌.
- **n이 작으므로:** 모든 boxplot은 개별 점(n=3)을 함께 표시. 막대그래프만 단독 사용 금지.
- **통계 표기:** 검정 이름·p·n을 캡션에 명시. permutation 결과는 null 분포 histogram + 관측값 세로선으로.
- **염색체:** 군간 비교는 상염색체만(성 교란 때문). 성염색체는 Fig 1d/Supp에서만.
- **출력:** 벡터(SVG/PDF), 패널 문자 a, b, c… 소문자 bold, 폰트 Arial/Helvetica 6–8 pt.

---

## Figure 1. Experimental design and a single long-read dataset that reads genome and epigenome together

**a. Study design (schematic).** Primary astrocytes pooled from multiple mice were cultured under four conditions — A, untreated; B, acute H2O2; C, acute H2O2 + drug; D, chronic H2O2 + drug —, each in three independent biological replicates processed in separate sequencing runs (run 1–3; n = 12). Genomic DNA from each well was sequenced on PacBio Revio (HiFi). **[ready — 모식도]**

**b. One dataset, multiple layers (schematic).** From the same HiFi reads: small variants (DeepVariant), structural variants and copy number (Sawfish), haplotype phasing (HiPhase), and base modifications from polymerase kinetics (5mC, 5hmC, 6mA; jasmine/pb-CpG-tools/modkit). Arrows indicate which layer feeds which figure. **[ready — 모식도]**

**c. Sequencing quality per sample.** HiFi read count (5.6–9.3 M), read-length N50 (14.3–16.9 kb), median read quality (Q29–Q36), mapping rate (≥ 99.9 %), mean autosomal depth (25–42×), and fraction of heterozygous variants phased (0.61–0.81) for the 12 samples (colour, condition; shape, run). **[ready — `/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/figure_data/sequencing_qc.tsv`]**
  - 주의(캡션 또는 Methods): 무처치(A)가 run 1·2에서 depth가 가장 낮고(25.4×, 29.8×) 세 run 모두에서 read quality가 가장 낮다(Q29–32 vs Q33–36). run 3의 A는 depth 38×로 다른 조건과 같은데도 이형접합 SNV 수와 메틸화 거리가 여전히 가장 커서 depth로 설명되지 않지만, 리뷰어 질문에 대비해 이 패널에서 보여준다.
  - phasing 비율은 이형접합 밀도를 따라 A에서 가장 높다(0.79 vs 0.62–0.72). stats_file의 phase_block_ng50 값은 대부분 0으로 기록돼 쓰지 않음.

**d. Composition of the starting cell populations.** Per-sample chrX/autosome and chrY/autosome depth ratios, and the number of autosomal heterozygous SNVs. A and D cultures contain male-derived cells (chrY/autosome ≈ 0.13, i.e. a male–female mixture), whereas B and C are female-only; heterozygous SNV counts differ up to four-fold between conditions. *(캡션 결론 문장은 §0-2 답에 따라 작성.)* **[ready — mosdepth summary, `TA/cohort/sample_metrics_combined.tsv`]**

---

## Figure 2. Treatment does not increase detectable mutational burden

**a. Variant counts per sample.** Autosomal PASS SNVs, indels, and SVs per sample. Differences in SNV counts track heterozygous-site density (Fig. 3) rather than new mutations. **[ready — `sample_metrics_combined.tsv`, 상염색체 재집계 확인]**

**b. Gene-level burden of protein-altering small variants.** Rank plot (or volcano) of HIGH/MODERATE-impact SNV/indel burden across conditions for 80 genes; no gene reaches FDR < 0.05 (top hit nominal p = 0.04, FDR = 1). **[ready — `TA/cohort/variant_burden_stats.csv`]**

**c. Gene-level burden of gene-disrupting SVs.** Same as b for SVs from single-sample calls (278 genes, 0 at FDR < 0.05). **[ready — `TA/cohort/sv_burden_stats.csv`]**

**d. Pathway-level burden.** Forest plot of 20 curated gene sets (oxidative stress, BER/NER/MMR/HR/NHEJ, apoptosis, ferroptosis, mitochondrion, …) for small variants and SVs; minimum FDR 0.33–0.39. **[ready — `geneset_burden_{variant,sv}_curated.csv`]**

**e. Mutational spectrum of sample-private SNVs.** Six-class substitution spectrum of SNVs called in only one of 12 samples with no alt-read support in the other 11 (alt reads ≥ 3, depth 15–100; autosomes). The C>A fraction follows treatment: untreated 7.5 %, acute H2O2 10.3 %, acute H2O2 + drug 8.6 %, chronic H2O2 + drug 11.0 % (quasi-binomial GLM condition + run, F = 15.3, p = 0.003; exact within-run permutation p = 0.0017, the minimum attainable; B vs A OR 1.38, C vs B OR 0.82, D vs C OR 1.31). Per 100,000 heterozygous sites, C>A rises (39.5 → 46.9 in acute H2O2, back to 40.0 with drug, 55.2 under chronic exposure) while C>T and T>C do not. Consistent with oxidative (8-oxoG-type) mutagenesis that the drug attenuates; the absolute excess is small (~10–16 C>A per sample). **[ready — `base_mod_per_sample/mutation_spectrum.{sh,py}`; 출력 `/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/figure_data/mutspec.*`]**
  - 주의(캡션 또는 Methods에 명시): private SNV 수가 샘플의 이형접합 SNV 수에 거의 비례(log-log 기울기 1.05, R² = 0.93)하므로 private 집합의 상당 부분은 체세포 변이가 아니라 웰마다 다른 개체 구성에서 온 생식세포 변이다. 그래서 절대 수가 아니라 클래스 비율과 이형접합 수로 나눈 값을 제시한다.
  - **SBS18 비교(완료):** H2O2 처리 조건(B+D 합산)에서 무처치 대비 늘어난 C>A의 3염기 문맥은 COSMIC SBS18(활성산소 8-oxoG 손상)과 코사인 0.76(C>A 16채널; 채널 섞기 기준 p = 0.04), SBS36(MUTYH 결손, 8-oxoG 수선 실패)과 0.80이다. 라이브러리 제작 중 산화 아티팩트 signature인 SBS45와는 0.58로 더 낮다. 약물 병용 조건(C)의 초과분은 SBS18과 닮지 않았다(0.41). 초과분이 이형접합 10만 개당 6–19개로 작아서 "SBS18과 부합" 수준으로만 쓴다. 패널 형태: 조건별 C>A 16채널 막대(이형접합 10만 개당) + SBS18 기준 막대. **[ready — `/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/figure_data/mutspec.CA16_per100k_het.tsv`, `/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/figure_data/sbs_signature_check.txt`]**
  - **아티팩트 필터 후에도 유지(완료):** depth > 샘플 평균 1.5배(다중 사본 정렬 의심) 또는 같은 샘플 후보가 10 kb 안에 있는 자리(군집)를 빼면 C>A 비율이 A 6.3 % → B 9.4 % → C 7.6 % → D 12.3 %로 조건 간 차이가 오히려 커진다. 정제 집합 통계(완료, `figure_data/mutspec.CA_stats.tsv`): 효과 크기는 커지지만(D vs A OR 2.08, GLM p = 0.007; B vs A OR 1.55, p = 0.049; D vs C OR 1.71, p = 0.036) 자리 수가 6,718개로 줄어 4조건 전체 검정은 경계 수준이다(F p = 0.034, 정확 permutation p = 0.064). VAF < 0.3에서도 같다(D vs A OR 2.08, p = 0.012; permutation p = 0.075). **캡션 권장 서술:** 주 결과는 정제 전 집합(permutation p = 0.0017)으로 제시하고, 정제 집합에서 방향과 효과 크기가 유지됨(효과 크기 증가, 전체 검정 p ≈ 0.06)을 함께 적는다. 약물 효과(C vs B)는 정제 집합에서 유의하지 않아(p = 0.26) "경향"으로만 쓴다. 표시 파일: `/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/figure_data/mutspec.private_sites.flagged.tsv.gz`.
  - IGV 목록: `figure_data/igv_CA_candidates.tsv`(변이 read 많은 순 120개, 판독 확인용 — 이 중 80개가 아티팩트 의심 표시), `figure_data/igv_CA_candidates_lowVAF_SBS18.tsv`(VAF 0.1–0.3, SBS18 주요 문맥, 표시된 자리 제외 125개 — 산화 손상 후보). 각각 `.batch.txt`로 IGV 스냅샷 자동화.

**f. Detection power.** With the observed between-replicate variation (residual SD of log private-SNV burden 0.18; 0.14 after normalizing to heterozygous-site count; df = 6), a 3-vs-3 run-paired comparison has 80 % power at α = 0.05 to detect a 1.65-fold (1.46-fold normalized) change in total burden. No condition differs at this scale (normalized burden 3.8–6.5 per 1,000 heterozygous sites, no treatment trend). **[ready]**

**g. Transcriptional strand asymmetry (Supplementary 후보).** Sample-private SNVs (artefact-filtered) in protein-coding/lncRNA genes split by transcribed (T) vs untranscribed (U) strand of the pyrimidine. Germline-dominated classes show the expected untranscribed-strand excess (T>C 57 % U, p < 10⁻⁴; C>T 54 %, p = 0.009), validating the annotation, whereas C>A is symmetric in every condition (123 U vs 123 T overall; treated vs untreated no difference) and does not depend on promoter-methylation-inferred activity. The treatment-associated C>A excess is therefore not transcription-coupled; power is limited (~50 C>A per condition). **[ready — `base_mod_per_sample/strand_asymmetry.py`, `figure_data/strand_asymmetry.*`]**
  - RNA-seq이 서버에 없어 발현은 프로모터 CpG 메틸화(12샘플 평균; <0.2 active, ≥0.6 silent)로 근사. 단일 분자 변이(요청서 ① 대체안)로 C>A 수가 늘면 재검정.

*Supplementary:* joint-called SV burden giving 112 "significant" genes that are all multi-copy gene families (alignment artefact) — 왜 single-sample을 채택했는지.

---

## Figure 3. Treatment is associated with progressive clonal homogenization of the cell population (core figure)

**a. Principle (schematic).** In a polyclonal bulk population heterozygous variants sit near VAF 0.5; when one lineage carrying loss of heterozygosity (LOH) expands, bulk VAF at those sites drifts toward 0 or 1. The fraction of heterozygous sites with VAF > 0.9 therefore measures dominance of a lineage. **[ready — 모식도]**

**b. VAF distributions.** Density of VAF at heterozygous SNVs for each sample, overlaid by condition (A → D). **[ready — `loh_evidence_svg/06_vaf_density_by_group.svg` 원자료]**

**c. Homozygous-peak fraction.** het_hom_peak_frac per sample (points, run-connected lines): A < C < B < D, reproduced in each of three runs; Kruskal–Wallis H = 9.5, p = 0.023. Acute H2O2 increases homogenization (B > A), the drug attenuates it (C < B), and chronic exposure produces the strongest effect despite the drug (D). **[ready — `sample_metrics_combined.tsv`]**

**d. Genotype discordance among replicates.** Fraction of sites with discordant genotypes across the three replicates of each condition from joint genotyping: A 0.53 > C 0.41 > B 0.40 > D 0.33. **[ready — `discordance_combined.tsv`]**

**e. Haplotype-directional allelic imbalance (long-read specific, independent evidence).** Fraction of phase blocks in which heterozygous SNVs are consistently skewed toward the same haplotype: A 0.12 < C 0.35 < B 0.63 < D 0.71. Unlike c–d, this metric depends on read-level phasing and not on per-site genotype counts. **[ready — `TA/cohort/hap_imbalance/*.tsv`]**

**f. Copy-number instability at multi-copy loci.** Range of integer copy number across replicates at seven multi-copy loci; mean instability rank A 2.93 > C 2.64 > B 2.36 > D 2.07 (descriptive; loci partly non-independent). **[ready — `copynum_instability_{rank,detail}.csv`]** *→ Supplementary 후보*

**g. Concordance across methods.** Heatmap of condition ranks (columns A–D) for each independent metric (rows c–f); all rows give the same order (Kendall's W = 1.0). Row annotation marks which metrics are statistically independent. **[ready — `heterogeneity_consensus_svg`]**

**h. Genome-wide map of condition-specific LOH.** Ideogram of LOH hotspots private to A (heterozygosity retained, n = 150) and to D (additional LOH, n = 65), with candidate genes labelled (Hsf2bp, Tusc3 shared by B+D; Piwil1 D-specific). Exploratory. **[ready — `loh_hotspot_cross_validation*.tsv`]**

*Supplementary:* HP-unassigned read fraction vs heterozygous-SNV count (ρ = −0.993) — 독립 증거가 아님을 명시; median het VAF의 depth 양자화 아티팩트.

---

## Figure 4. DNA methylation (5mC): globally stable, with focal differences that mirror population homogenization

**a. Global CpG methylation.** Mean 5mC across 18.8 M autosomal CpGs covered ≥10× in all samples; 0.760–0.774, no condition effect. **[ready — `cpg_group_signal_check.py` 출력]**

**b. Principal component analysis.** PCA of the 200,000 most variable CpGs. B and D replicates cluster tightly across runs; A replicates (and one C) are dispersed and each defines its own component (PC1 15.2 %, PC2 12.1 %, PC3 10.2 %). **[ready]**

**c. A real but modest condition effect.** Number of CpGs with condition effect p < 10⁻³ (model condition + run; F-test over B−A, C−B, D−C) observed (vertical line) versus 100 within-run label permutations (histogram); observed exceeds all permutations (empirical p = 0.02), with ~15 % excess over the null median. **[ready → 정식 버전은 exact permutation(576 조합)으로 교체 권장 — analysis]**

**d. Replicate dispersion of methylation by condition.** Depth-corrected mean squared CpG methylation difference between replicates of the same condition (18.8 M CpGs, binomial sampling variance subtracted): untreated 8.08 × 10⁻³ (all three pairs 7.4–8.7) versus 6.29–6.71 × 10⁻³ for treated conditions; unchanged at ≥ 20× coverage (9.2 M CpGs). Untreated replicates differ from each other more than replicates of any treated condition do, and more than samples from different conditions on average (7.13 × 10⁻³), whereas treated replicates are closer to each other than to other conditions. Mirrors the genetic heterogeneity in Fig. 3. **[ready — `base_mod_per_sample/cpg_replicate_dispersion.py`; 출력 `/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/figure_data/cpg_dispersion.*`]**
  - 그림 형태: 12×12 거리 heatmap(조건 순 정렬) + 조건별 within-pair 점 그림. 처치군 간 순서(B 6.71 > C 6.42 > D 6.29)는 차이가 작아 해석하지 않음.

**e. Region-level differential methylation.** Manhattan-style plot of DMRs called with the full design (condition + run), significance calibrated by within-run permutation of region statistics; top autosomal candidate (chr6:28.93 Mb, shared by B vs D and C vs D) shown as a methylation track with haplotype-resolved reads. **[analysis — 기존 pairwise DSS는 p-value 없음, 재분석 필요]**

**f. Methylation heterogeneity tracks genetic heterogeneity.** Per sample, mean depth-corrected methylation distance to its two same-condition replicates versus heterozygous-SNV count (Spearman ρ = 0.71, p = 0.010; run-centred r = 0.87, p = 0.0002) and versus homozygous-peak fraction (ρ = −0.62, p = 0.031). **[ready — `/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/figure_data/methylation_vs_clonality.tsv`]**
  - 해석 주의: 12개 점이 조건별로 묶여 있어 상관은 주로 조건 간 차이(A 대 나머지)에서 나온다. 또 메틸화 차이 중 일부는 유전형에 따라 메틸화가 달라지는 자리(meQTL, allele-specific methylation) 때문일 수 있어서, "세포 구성이 다양할수록 메틸화도 다양하다"까지만 쓰고 원인 방향은 주장하지 않는다.
  - 선택: 4d와 합쳐 한 패널로 쓰거나 Supplementary.

*Supplementary:* allele-specific methylation (ASM) region counts track heterozygous-site density (A 87–108 k → D 19–22 k) — 독립 신호가 아님을 명시; 성염색체 DMR이 상위 후보를 차지했던 것(성 교란) 제시.

---

## Figure 5. 5hmC and 6mA from the same reads

**a. Global 5hmC at CpGs per sample.** **[pending — jasmine 26.1.3 재호출 진행 중, run 1은 pileup만으로 약 2일 내, 전체 약 18일]**

**b. 5hmC by genomic element** (promoter, exon, intron, intergenic, CpG island) per condition. Biological rationale: TET-mediated 5mC oxidation is sensitive to oxidative stress. **[pending]**

**c. 5mC–5hmC relationship** at gene bodies per condition. **[pending]** **[decide]**

**d. Technical characterization of 6mA calls.** Per-sample 6mA call rate (~1 % of adenines at threshold 0.7 in the first chromosome examined — far above credible levels of native mammalian 6mA), its dependence on sequence context (5-mer; e.g. GGACT/CGACC enriched), mitochondrial versus nuclear rate, and fraction of reads from mycoplasma. Establishes that most calls reflect kinetic background rather than biology. **[pending — 6mA pileup·집계 진행 중(약 1–2일); 마이코플라스마 점검은 사용자 확인 대기]**
  - **마이코플라스마 점검(완료):** 12샘플 × read 10만 개(총 120만)에서 마이코플라스마 7종(*M. orale, hyorhinis, arginini, fermentans, hominis, salivarium, A. laidlawii*)에 정렬된 read 0개. 마우스에 정렬되지 않는 read는 0.006–0.04 %이고 GC가 낮은 미정렬 read는 샘플당 0–8개. 오염은 6mA 신호의 출처가 아님(검출 한계 약 0.003 % read). read당 6mA 밀도는 마우스 핵 read에서 12샘플 모두 0.56–0.65 %로 균일, chrM read에서는 더 낮음. `figure_data/mycoplasma_summary.tsv`

**e. 6mA enrichment by repeat class and condition.** Relative 6mA rate (feature / sample global) for RepeatMasker classes (e.g. satellite ~2.6× enriched in first data) and gene-level aggregation, model condition + run, permutation-calibrated. **[pending]** **[decide — 신호가 없거나 기술적이면 Fig 5에서 빼고 Supplementary로]**

> Fig 5 구성은 결과에 따라 결정: 5hmC에 그룹 신호가 있으면 Fig 5 = 5hmC 단독, 6mA는 Supplementary(기술적 한계 + DNA 손상 관련 kinetics 이상 가능성 언급).

---

## Figure 6. Model

**Summary schematic.** H2O2 does not measurably increase point mutations or SVs (Fig. 2) but is associated with reduced population diversity — fewer lineages dominating the culture (Fig. 3); the drug partially prevents this under acute exposure but not under chronic exposure. The loss of diversity is mirrored by convergent methylation states in acute H2O2 and chronic H2O2 + drug cultures (Fig. 4). Bulk epigenomic differences between conditions should therefore be interpreted in light of population composition. *(문구는 §0-2 답에 따라 "selection" 또는 "different starting populations" 중 선택.)* **[모식도]**

---

## 추가 분석 목록 (그림을 그리기 전에 필요)

| 패널 | 작업 | 규모 |
|---|---|---|
| 1c | ~~QC 지표~~ 완료 | — |
| 2e | ~~변이 스펙트럼, SBS18~~ 완료; IGV 확인 남음 | 짧음 |
| 2f | ~~검출력~~ 완료 | — |
| 4c | exact within-run permutation (576 조합) | 서버 여유 시 수 시간 |
| 4d | ~~반복 샘플 간 메틸화 거리~~ 완료 | — |
| 4e | condition + run 모델 DMR + permutation 보정 | 1–2일 |
| 5a–e | 5hmC/6mA 계산 완료 후 | 대기 |
