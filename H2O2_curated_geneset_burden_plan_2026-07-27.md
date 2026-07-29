# H2O2 WGS — 큐레이션 GO 세트 burden test 확대 계획

**작성일:** 2026-07-27
**배경:** 기존 GO burden은 산화스트레스·DNA repair **2개 세트만** 검정했고 유의한 것이 없었다([`H2O2_tertiary_analysis_results`](./H2O2_tertiary_analysis_results_2026-07-16.md) §3, Slide 2-5b). 개별 유전자 burden도 FDR 보정 후 0개였다. 본 작업은 **"가장 그럴듯한 2개만이 아니라 넓게 훑어도 유의한 경로가 없다"**는 완결성(negative-result rigor)을 확보하기 위해 큐레이션된 수십 개 GO 세트로 확대한다.

**중요 프레이밍**: 이건 **발견(discovery)이 아니라 완결성 확인**이다. n=3 + per-gene 노이즈 수준 결과를 고려하면 진짜 히트가 나올 가능성은 낮고, 나오면 오히려 위양성/아티팩트를 먼저 의심한다. 결과는 "넓게 봐도 없다"를 문서화하는 데 의의가 있다.

---

## 1. 방법 요약 (기존 파이프라인 재사용, 구조 변경 최소)

- `fetch_go_genesets.py`의 GO ID 매핑만 확대 → `resources/genesets_curated.json` 생성
- `geneset_burden.R`·Snakefile rule은 **세트 수에 무관하게 동작** → 코드 변경 불필요
- 기존 2-세트 출력(`geneset_burden_*.csv`)은 **건드리지 않고**, 큐레이션 결과는 `*_curated.csv`로 별도 저장(비파괴)
- 검정 방식은 기존과 동일: mode=variant는 `HIGH|MODERATE` 영향 SNV/indel, mode=sv는 `sv:cds/sv:utr`; 세트별 presence(Fisher/카이제곱) + count 기반 검정, FDR 보정, 효과크기(Cramér's V, max_prop_diff) 병기

## 2. 큐레이션 GO 세트 (H2O2 성상세포 생물학 겨냥, 다카피 아티팩트 회피)

산화스트레스·게놈안정성·세포사멸·신경/성상세포 테마 위주. **후각수용체·면역글로불린·vomeronasal 등 다카피 아티팩트 유발 GO는 의도적으로 제외**.

| 세트명 | 대표 GO ID | 테마 |
|---|---|---|
| oxidative_stress_response | GO:0006979, 0034599, 0000302, 0034614 | 산화스트레스/ROS 반응 (기존) |
| glutathione_redox | GO:0006749, 0045454, 0016209 | 글루타티온/redox 항상성 |
| dna_repair_core | GO:0006281, 0006974 | DNA 수선/손상반응 (기존) |
| dna_repair_ber | GO:0006284 | 염기절제수선(BER) |
| dna_repair_ner | GO:0006289 | 뉴클레오타이드절제수선(NER) |
| dna_repair_mmr | GO:0006298 | 불일치수선(MMR) |
| dna_repair_dsb_hr | GO:0000724 | 이중가닥절단-상동재조합(HR) |
| dna_repair_dsb_nhej | GO:0006303 | 이중가닥절단-NHEJ |
| dna_damage_checkpoint | GO:0031570, 0000077 | DNA 손상 체크포인트 |
| apoptosis | GO:0006915, 0043065 | 세포자멸 |
| intrinsic_apoptosis | GO:0097193, 0008630 | 내인성 세포자멸(DNA손상 유발) |
| ferroptosis | GO:0097707 | 페롭토시스(산화적 세포사) |
| autophagy | GO:0006914, 0016236 | 자가포식 |
| cellular_senescence | GO:0090398, 0007050 | 세포노화/세포주기정지 |
| cell_cycle | GO:0007049 | 세포주기 |
| er_stress_upr | GO:0034976, 0030968, 0006986 | ER 스트레스/UPR |
| mitochondrion | GO:0006119, 0007005 | 미토콘드리아/산화적인산화 |
| inflammatory_response | GO:0006954 | 염증반응(성상세포 반응성) |
| gliogenesis_astrocyte | GO:0042063, 0014002, 0048708 | 신경교/성상세포 분화 |
| hypoxia_response | GO:0001666, 0036293 | 저산소 반응 |
| calcium_signaling | GO:0006816, 0070588 | 칼슘 신호(성상세포) |

→ 약 21개 세트. (GAF에 직접 주석된 유전자만 모으므로, GAF에 없는 GO ID는 자동으로 빈 세트가 되어 결과에서 NA로 처리됨.)

## 3. 안전장치 / 정직성

1. **다카피 아티팩트 회피**: 후각·면역글로불린·vomeronasal GO 미포함(구성 단계에서 배제).
2. **NA 예상 명시**: variant 모드는 `HIGH|MODERATE` 필터라 많은 세트가 presence=0(전군 동일)→검정불가(NA)로 나올 것. 이는 오류가 아니라 "그 세트에 군 가를 고영향 변이 없음".
3. **효과크기 병기**: p값만이 아니라 Cramér's V·max_prop_diff를 함께 보고(n=3에서 p 과신 방지).
4. **FDR 전체 세트에 적용**: 세트 수 증가분까지 포함해 BH-FDR 재계산.
5. **GO closure 미전개 한계**: fetch 스크립트는 상/하위 계층을 펼치지 않고 직접 주석만 모음(기존과 동일 한계) — 부모 term은 일부 하위 유전자 누락 가능. 완결성 negative-control 목적엔 허용, 문서에 명시.

## 4. 실행 단계

1. `fetch_go_genesets.py`에 `--preset {original|curated}` 추가(기본 original → 기존 동작 불변), curated 컬렉션 정의.
2. `resources/mgi.gaf.gz`로 `resources/genesets_curated.json` 생성, 세트별 유전자 수 확인.
3. `geneset_burden.R`를 variant/sv 두 모드로 실행 → `geneset_burden_variant_curated.csv`, `geneset_burden_sv_curated.csv`(+plots).
4. 결과 요약: 세트별 p/FDR/효과크기 표, 유의(FDR<0.05) 세트 유무, NA 세트 수.

## 5. 예상 결과 & 해석 가이드

- **예상**: 유의(FDR<0.05) 세트 0개 또는 극소수. 만약 나오면 (1) 다카피 잔존 오염, (2) 위양성, (3) 성염색체 잔존을 먼저 점검.
- **결론 문구(예상)**: "산화스트레스·게놈안정성·세포사멸·신경 관련 ~21개 GO 세트로 확대해도 군을 가르는 유의한 경로 burden은 없었다 → 개별 유전자·경로 어느 층위에서도 범인이 없고, 신호는 게놈 전역 판세(Part 3)에 있다"는 기존 결론을 강화.

---

## 6. 결과 (완료, 2026-07-27)

`resources/genesets_curated.json`(20개 세트; cell_cycle는 GAF 직접주석 없어 자동 제외) 생성 후 두 모드 실행:
- 출력: `geneset_burden_variant_curated.csv`, `geneset_burden_sv_curated.csv` (+plots). 기존 2세트 출력은 그대로 보존.

| 모드 | 세트 수 | 검정가능 / 전부NA | **FDR<0.05 유의** | 최소 FDR | 최소 p(보정 전) |
|---|---|---|---|---|---|
| VARIANT (SNV+indel, HIGH\|MODERATE) | 20 | 3 / 17 | **0개** | 0.392 | apoptosis p_count=0.306 (Cramér's V=0.56) |
| SV (sv:cds/sv:utr) | 20 | 15 / 5 | **0개** | 0.325 | **mitochondrion p_count=0.036** → FDR=0.325 |

**해석:**
1. **어느 세트도 FDR<0.05를 통과하지 못함** — 두 모드 모두 유의 경로 0개. 기존 2세트 결과와 완전히 일치.
2. **SV에서 mitochondrion만 명목상 p=0.036**이었으나 20세트 다중검정 보정 후 FDR=0.325로 소멸 — 20개 검정 중 하나쯤 p<0.05가 나오는 우연 기대치 안(예상된 위양성). 효과크기도 특기할 것 없음.
3. **VARIANT는 17/20이 NA** — `HIGH|MODERATE` 필터 후 그 세트에 군 가를 고영향 SNV/indel이 없어 presence가 전군 동일(분산 없음)→검정불가. "신호 없음"의 또 다른 표현.
4. ferroptosis(유전자 1개), cell_cycle(직접주석 0) 등은 애초 검정력 없음 — 완결성 목록에 포함했으나 결과엔 기여 없음.

**결론(확정):** 가장 그럴듯한 2개가 아니라 **산화스트레스·글루타티온·DNA수선 하위경로(BER/NER/MMR/HR/NHEJ)·손상체크포인트·세포자멸·페롭토시스·자가포식·노화·ER스트레스·미토콘드리아·염증·신경교분화·저산소·칼슘신호 등 20개 경로로 넓혀도 군을 가르는 유의한 경로 burden은 없다.** → 개별 유전자든 기능 경로든 어느 층위에서도 "범인"이 없고, 신호는 게놈 전역 판세(Part 3 이질성 그라디언트)에 있다는 결론이 **완결성 있게 강화됨**. (단, GO closure 미전개·n=3 검정력 한계는 여전하며, 이는 "특정 경로 배제"가 아니라 "표적·확대 검정 모두 음성"으로 해석.)
