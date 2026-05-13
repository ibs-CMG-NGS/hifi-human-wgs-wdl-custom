# De Novo Assembly 분석 가이드

## 개요

`run_assembly=true` 플래그로 기존 reference-based 파이프라인과 동일한 HiFi reads로 de novo assembly 및 assembly QC를 수행한다. Reference-based 분석(SNV, SV, methylation 등)과 독립적으로 진행되며, 동일한 `miniwdl run` 호출 내에서 함께 실행된다.

### 분석 단계

```
HiFi reads (uBAM)
  │
  ├─ hifiasm ──────────────────────────► GFA (primary / hap1 / hap2)
  │     └─ gfa_to_fasta ──────────────► FASTA + 기본 통계 (size, contig count)
  │           │
  │           ├─ BUSCO (×3) ──────────► conserved gene 기반 completeness
  │           │
  │           ├─ meryl_count ──────────► k-mer database (tar.gz)
  │           │     └─ Merqury ───────► QV score + k-mer completeness
  │           │
  │           └─ minimap2_asm_to_ref   ► PAF (hap1, hap2 → reference)
  │                 └─ SyRI (×2) ─────► inversion / translocation / duplication
  │
  └─ QUAST ────────────────────────────► N50, NG50, misassemblies, genome fraction
```

---

## 활성화 방법

`inputs.json`에 플래그 추가:

```json
{
  "humanwgs_singleton.run_assembly": true,
  "humanwgs_singleton.busco_lineage": "mammalia_odb10"
}
```

또는 커맨드라인 오버라이드:

```bash
miniwdl run workflows/singleton.wdl \
  --input TG6283.inputs.json \
  --dir /path/to/output \
  humanwgs_singleton.run_assembly=true
```

기본값은 `run_assembly=false`이므로 기존 분석은 영향 없음.

---

## 파라미터

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `run_assembly` | `false` | de novo assembly + QC 실행 여부 |
| `busco_lineage` | `mammalia_odb10` | BUSCO lineage dataset |

### 종별 BUSCO lineage 권장값

| 종 | 권장 lineage |
|----|-------------|
| 마우스 (Mus musculus) | `mammalia_odb10` 또는 `murinae_odb10` |
| 인간 (Homo sapiens) | `primates_odb10` 또는 `vertebrata_odb10` |
| 기타 포유류 | `mammalia_odb10` |
| 척추동물 | `vertebrata_odb10` |
| 진핵생물 | `eukaryota_odb10` |

---

## 출력 파일

모든 출력은 `out/` 하위 디렉터리에 저장되며, 파일명 형식은 `{sample_id}.hifiasm.bp.{hap}.p_ctg.{analysis}.{ext}`.

### Assembly 파일

| 출력 변수 | 파일명 패턴 | 설명 |
|-----------|------------|------|
| `assembly_primary_gfa` | `{id}.hifiasm.bp.p_ctg.gfa` | Primary assembly GFA |
| `assembly_hap1_gfa` | `{id}.hifiasm.bp.hap1.p_ctg.gfa` | Haplotype 1 GFA |
| `assembly_hap2_gfa` | `{id}.hifiasm.bp.hap2.p_ctg.gfa` | Haplotype 2 GFA |
| `assembly_primary_fasta` | `{id}.hifiasm.bp.p_ctg.fa` | Primary assembly FASTA |
| `assembly_hap1_fasta` | `{id}.hifiasm.bp.hap1.p_ctg.fa` | Haplotype 1 FASTA |
| `assembly_hap2_fasta` | `{id}.hifiasm.bp.hap2.p_ctg.fa` | Haplotype 2 FASTA |

### BUSCO QC

| 출력 변수 | 파일명 패턴 | 설명 |
|-----------|------------|------|
| `assembly_busco_primary` | `{id}.hifiasm.bp.p_ctg.busco.short_summary.txt` | Primary BUSCO 결과 |
| `assembly_busco_hap1` | `{id}.hifiasm.bp.hap1.p_ctg.busco.short_summary.txt` | Hap1 BUSCO 결과 |
| `assembly_busco_hap2` | `{id}.hifiasm.bp.hap2.p_ctg.busco.short_summary.txt` | Hap2 BUSCO 결과 |

### Merqury QC

| 출력 변수 | 파일명 패턴 | 설명 |
|-----------|------------|------|
| `assembly_merqury_qv` | `{id}.hifiasm.merqury.qv` | QV scores |
| `assembly_merqury_stats` | `{id}.hifiasm.merqury.completeness.stats` | K-mer completeness |

### Assembly-to-Reference (minimap2 + SyRI + QUAST)

| 출력 변수 | 파일명 패턴 | 설명 |
|-----------|------------|------|
| `assembly_minimap2_hap1_paf` | `{id}.hifiasm.bp.hap1.p_ctg.asm_to_ref.paf` | Hap1 → ref 정렬 PAF |
| `assembly_minimap2_hap2_paf` | `{id}.hifiasm.bp.hap2.p_ctg.asm_to_ref.paf` | Hap2 → ref 정렬 PAF |
| `assembly_syri_hap1_out` | `{id}.hifiasm.bp.hap1.p_ctg.syri.out` | Hap1 구조변이 전체 목록 |
| `assembly_syri_hap1_summary` | `{id}.hifiasm.bp.hap1.p_ctg.syri.summary` | Hap1 구조변이 요약 |
| `assembly_syri_hap2_out` | `{id}.hifiasm.bp.hap2.p_ctg.syri.out` | Hap2 구조변이 전체 목록 |
| `assembly_syri_hap2_summary` | `{id}.hifiasm.bp.hap2.p_ctg.syri.summary` | Hap2 구조변이 요약 |
| `assembly_quast_report` | `{id}.hifiasm.bp.p_ctg.quast.report.tsv` | QUAST 품질 리포트 (TSV) |
| `assembly_quast_html` | `{id}.hifiasm.bp.p_ctg.quast.report.html` | QUAST 인터랙티브 리포트 |

---

## Stats TSV 추가 항목

`run_assembly=false`일 때는 모두 `"0"`.

### Assembly 기본 통계

| 컬럼 | 설명 |
|------|------|
| `assembly_primary_size` | Primary 어셈블리 총 크기 (bp) |
| `assembly_primary_ctgs` | Primary contig 수 |
| `assembly_hap1_size` | Hap1 어셈블리 총 크기 (bp) |
| `assembly_hap1_ctgs` | Hap1 contig 수 |
| `assembly_hap2_size` | Hap2 어셈블리 총 크기 (bp) |
| `assembly_hap2_ctgs` | Hap2 contig 수 |

### BUSCO (primary 기준)

| 컬럼 | 설명 |
|------|------|
| `busco_primary_complete` | Complete BUSCO % |
| `busco_primary_single` | Complete & single-copy % |
| `busco_primary_missing` | Missing BUSCO % |

### Merqury

| 컬럼 | 설명 |
|------|------|
| `merqury_qv` | Assembly QV score |
| `merqury_completeness` | K-mer completeness (%) |

### QUAST

| 컬럼 | 설명 |
|------|------|
| `quast_n50` | Assembly N50 (bp) |
| `quast_ng50` | Assembly NG50 대비 reference (bp) |
| `quast_misassemblies` | Misassembly 수 |
| `quast_genome_fraction` | Genome fraction covered (%) |

### SyRI (hap1, hap2 각각)

| 컬럼 | 설명 |
|------|------|
| `syri_hap1_inv` | Hap1 inversion 수 |
| `syri_hap1_trans` | Hap1 translocation 수 |
| `syri_hap1_dup` | Hap1 duplication 수 |
| `syri_hap2_inv` | Hap2 inversion 수 |
| `syri_hap2_trans` | Hap2 translocation 수 |
| `syri_hap2_dup` | Hap2 duplication 수 |

---

## 도구 및 컨테이너

| 도구 | 버전 | 컨테이너 | 용도 |
|------|------|----------|------|
| hifiasm | 0.19.5 | `quay.io/biocontainers/hifiasm:0.19.5--h5b5514e_1` | De novo assembly |
| BUSCO | 5.7.1 | `quay.io/biocontainers/busco:5.7.1--pyhdfd78af_1` | Completeness QC |
| Merqury/meryl | 1.3 | `quay.io/biocontainers/merqury:1.3--hdfd78af_4` | K-mer QC |
| minimap2 | 2.30 | `quay.io/biocontainers/minimap2:2.30--h577a1d6_0` | Assembly → ref 정렬 |
| SyRI | 1.7.1 | `quay.io/biocontainers/syri:1.7.1--py311hcf77733_1` | 구조변이/재배열 탐지 |
| QUAST | 5.3.0 | `quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2` | Assembly 품질 평가 |

---

## 리소스 요구사항 (마우스 기준 순차 실행)

| Task | CPU | RAM | 예상 소요 시간 |
|------|-----|-----|---------------|
| hifiasm | 32 | 64 GB | ~3h |
| gfa_to_fasta (×3) | 2 | 4 GB | ~5min |
| BUSCO (×3) | 8 | 16 GB | ~1h each |
| meryl_count | 16 | 32 GB | ~30min |
| Merqury | 8 | 16 GB | ~30min |
| minimap2 (×2) | 8 | 16 GB | ~30min each |
| SyRI (×2) | 4 | 8 GB | ~15min each |
| QUAST | 8 | 16 GB | ~20min |
| **전체 추가 시간** | — | — | **~8–10h** |

> `task_concurrency=1` 환경(이 서버)에서는 모든 task가 순차 실행됩니다.

---

## 설계 노트

- **SyRI는 hap1/hap2만 실행**: primary assembly는 collapsed diploid로 per-haplotype SV 해석이 덜 명확하여 phased haplotype에만 적용
- **minimap2 preset `asm5`**: 동일 종(intra-species) 비교 표준값; 이종 비교 시 `asm20`으로 변경
- **QUAST `--no-icarus`**: 대용량 Icarus HTML 뷰어 비활성화로 불필요한 파일 생성 방지
- **BUSCO offline fallback**: HPC 오프라인 환경에서 `--offline` 먼저 시도 후 실패 시 온라인 다운로드
- **run_assembly=false 기본값**: 기존 reference-based 파이프라인 동작에 영향 없음
