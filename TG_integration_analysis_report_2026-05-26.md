# Transgene Integration Site Analysis Report

**Date:** 2026-05-26  
**Analyst:** ygkim  
**Pipeline:** PacBio HiFi WGS — hifiasm de novo assembly + hybrid reference  
**Reference:** GRCm39 (mouse) + Gencode vM36 annotation

---

## Summary

| Sample | Transgene | Breakpoint (GRCm39) | Gene | Copy # | Zygosity | Status |
|--------|-----------|----------------------|------|--------|----------|--------|
| TG-6283 B01 | gfa2+CreERT2+SV40pA | chr2:5,559,779 | Camk1d intron | 2 (tandem) | Possibly homozygous | Complete |
| TG-6102 A01 | Aldh1l1-CreERT2 | chr6:90,534,206 | Aldh1l1 exon (+22 bp) | 1 | Heterozygous (hap2) | Complete |
| TG-6425 C01 | Aldh1l1-EGFP | — | — | — | — | Pipeline running |

---

## 1. TG-6283 B01

### 1.1 Analysis Overview

- **HiFi reads:** 36x depth, hifiasm v0.25.0
- **Method:** De novo assembly → chimeric contig identification → hybrid reference with TG FASTA → chimeric HiFi read detection

### 1.2 Breakpoint

| Parameter | Value |
|-----------|-------|
| Chromosome | chr2 |
| Left breakpoint | 5,559,779 |
| Right breakpoint | 5,559,783 |
| Reference deletion | 4 bp |
| Support | **20 chimeric reads** + 44 TG-mapped reads |

### 1.3 Insert Structure

```
chr2 host ─────────[5.559 Mb]──────────────────────────────
                         ↓ insertion
Contig: ──host──|← TG copy 2 (RC) ←|← 403 bp spacer ←|← TG copy 1 (RC) ←|──host──
                  5,315 bp                               5,315 bp
                 Total insert: ~10,681 bp (+ 403 bp spacer = ~11,084 bp)
```

- **Copy number:** 2 tandem copies, both in reverse complement orientation
- **Copy coordinates (hap1 contig h1tg000090l):**
  - Copy 1: contig[2,531,427–2,536,742] → chr2 ~5,552,416 (interpolated)
  - Copy 2: contig[2,537,145–2,542,461] → chr2 ~5,558,076 (interpolated)
- **Spacer between copies:** 403 bp
- **TG FASTA coverage:** 100% for both copies (correct FASTA confirmed)

### 1.4 Genomic Context

- **Gene:** Camk1d (Calcium/calmodulin-dependent protein kinase ID)
- **Feature:** Intronic insertion (large ~120 kb intron)
- **Expression:** Neuronal — predominantly expressed in brain and nervous system
- **Insertion is intergenic to adjacent exons** — no direct exon disruption

### 1.5 Haplotype Structure

| | hap1 (h1tg000090l) | hap2 (h2tg000084l) |
|--|--|--|
| Contig length | 3,048,605 bp | 3,042,992 bp |
| Ref alignment blocks | 1 (single, mapq=60) | 3 (fragmented) |
| TG copies | 2 | 2 |
| Interpretation | Clean insertion | Additional structural variant possible |

The fragmented alignment of hap2 (3 blocks vs hap1's single block) suggests a possible additional structural variant near the insertion site on hap2, independent of the TG insertion.

---

## 2. TG-6102 A01

### 2.1 Analysis Overview

- **HiFi reads:** hifiasm de novo assembly
- **Initial analysis:** Failed — 0 chimeric reads, 37.5% TG coverage (only 1,967/5,239 bp matched)
- **Root cause:** TG FASTA (`Aldh1l1_CreERT2_TG-6102.fa`) mismatched actual insert by 57%
- **Solution:** Extracted actual TG insert (2,563 bp) from chimeric contig → rebuilt hybrid reference

### 2.2 TG FASTA Mismatch Analysis

| | Value |
|--|--|
| FASTA header | `260427_Aldh1l1CreERT2_full_reconstructed_5UTR+CreERT2+SV40pA+Aldh1l1_CDS_3UTR` |
| FASTA length | 5,239 bp |
| Assembly-matched bp | 2,253 bp (43.0%) |
| Unmatched bp | 2,986 bp (**57.0%**) |

**Coverage map (FASTA positions):**
```
1                                                                        5,239
|                                                                            |
···██████████████████████████████░·········································░████
   200        2,166          4,954                                       5,239

Matched:   200–2,166   (1,967 bp, 99.34% identity — 5'UTR + CreERT2 N-terminal)
           4,954–5,239 (286 bp, 100.0% — Aldh1l1 3'UTR end)
Missing:   2,167–4,953 (2,787 bp — CreERT2 C-terminal/hinge + SV40pA estimated)
```

**Interpretation:** The FASTA is labeled `full_reconstructed`, suggesting it was built in silico from reference sequences rather than sequenced directly from the actual plasmid. The ~2,787 bp gap (CreERT2 C-terminal to SV40pA region) does not match the inserted sequence.

### 2.3 Extracted Insert Sequence

Since the original FASTA was unreliable, the actual TG insert was extracted from chimeric contig h2tg000140l:

| Parameter | Value |
|-----------|-------|
| Contig | h2tg000140l (hap2, 189,979 bp) |
| Insert position in contig | 82,886–85,449 |
| Extracted length | **2,563 bp** |
| Host flanking (left) | chr6:90,534,279–90,542,165 (contig region: 85,449–93,000) |
| Host flanking (right) | chr6:90,526,421–90,534,206 (contig region: 51,339–82,886) |
| Output FASTA | `h2tg000140l_82886_85449_TG_insert.fa` |

Notable sequence features in extracted insert:
- NheI restriction site (GCTAGC) at position 1
- NotI restriction site (GCGGCCGGCC) near 3' end
- These are characteristic cloning sites for CreERT2 expression vectors

**Caveat — circular reasoning risk:**  
The extracted insert was derived from the same HiFi reads that produced the assembly. Therefore:
- Breakpoint coordinates are reliable (host-flanking sequence from reference)
- Insert sequence identity requires independent validation (Sanger/short-read sequencing of actual plasmid)

### 2.4 Breakpoint (from extracted FASTA hybrid reference)

| Parameter | Value |
|-----------|-------|
| Chromosome | chr6 |
| Breakpoint | **90,534,206** (±4 bp) |
| Reference gap | chr6:90,534,206 / 90,534,279 **(73 bp deletion)** |
| Range | chr6:90,519,758–90,549,328 |
| Support | **17 chimeric reads** |

The 73 bp reference deletion at the insertion site is characteristic of CRISPR-mediated targeted insertion (indel at cut site), consistent with the Aldh1l1 locus-directed knock-in design.

### 2.5 Genomic Context

- **Gene:** Aldh1l1 (Aldehyde dehydrogenase 1 family member L1)
- **Feature:** **Exon** — insertion is 22 bp from exon boundary
- **Design intent:** Knock-in at Aldh1l1 exon (intended to place CreERT2 under Aldh1l1 promoter)
- **Zygosity:** Heterozygous — hap1 has no TG signal; hap2 only

### 2.6 Pending Validation

- [ ] Obtain original Aldh1l1-CreERT2 construct sequence (Addgene accession, GenBank, or plasmid map)
- [ ] BLAST extracted 2,563 bp sequence to identify which CreERT2 variant was used
- [ ] Confirm 73 bp deletion by PCR across insertion junction
- [ ] Re-run full integration pipeline with verified TG FASTA

---

## 3. Pipeline Bug Fix

### `extract_integration_coords` — incorrect breakpoint coordinate

**File:** `workflows/wdl-common/wdl/tasks/transgene_integration.wdl`  
**Commit:** `bed8cd4`

#### Problem

The original awk logic took the **midpoint** of the first PAF alignment block with nmatch ≥ min_match_bp. For a chimeric contig spanning 9+ Mb on chr6, this produced a wildly incorrect position:

```
chr6:81,268,183 – 90,534,206  →  midpoint = chr6:85,901,194  ✗
```

This caused `extract_region_reads` to extract reads from chr6:85.6–86.2 Mb (wrong), returning 0 chimeric reads.

#### Fix

Replace midpoint logic with gap-finding logic: sort PAF by contig query position, find consecutive alignment blocks where the **gap between them** (500–50,000 bp) corresponds to the TG insert, use the reference coordinate at the gap boundary.

```bash
# Before (wrong):
awk '$10 >= min { print $6 > "integration_chr.txt"; \
                  print int(($8+$9)/2) > "integration_pos.txt"; exit }'

# After (correct):
sort -k1,1 -k3,3n "$PAF" | awk -v min_gap=MIN '
  prev_qname != "" && $1 == prev_qname {
    gap = $3 - prev_qend
    if (gap >= min_gap && gap <= 50000) {
      chr = prev_tname; pos = prev_tend  # right boundary of left block
      print chr > "integration_chr.txt"
      print pos > "integration_pos.txt"
      exit
    }
  }
  { prev_qname=$1; prev_qend=$4; prev_tname=$6; prev_tstart=$8; prev_tend=$9; prev_strand=$5 }
' -
```

Result: correctly returns chr6:90,534,206 instead of chr6:85,901,194.

---

## 4. Biological Side-Effect Assessment

### TG-6283 B01 — Camk1d Intronic Insertion

| Risk factor | Assessment |
|-------------|------------|
| Insertion type | Intron (large, ~120 kb) |
| Gene function | Neuronal kinase (Camk1d) |
| Exon disruption | None directly |
| Splice site risk | Low (far from exon boundaries) |
| Zygosity | Possibly homozygous (both haplotypes show insertion) |
| TG expression | gfa2+CreERT2+SV40pA — inducible system |

**Conclusion: Moderate concern**  
Intronic insertion with intact exons is generally low-risk. However:
- Possibly homozygous — both copies of Camk1d may carry the insertion
- Camk1d is expressed in neurons; disruption of intronic regulatory elements is possible
- Recommend: RNA-seq of brain tissue to confirm Camk1d expression is unaffected

### TG-6102 A01 — Aldh1l1 Exonic Knock-in

| Risk factor | Assessment |
|-------------|------------|
| Insertion type | Exon (22 bp from boundary) |
| Gene function | Aldehyde dehydrogenase (astrocyte marker) |
| Exon disruption | **1 allele only (heterozygous)** |
| Design intent | Knock-in — CreERT2 driven by endogenous Aldh1l1 promoter |
| Zygosity | Heterozygous by design |
| Wild-type allele | Intact (hap1) |

**Conclusion: Low concern**  
This is an intended knock-in design. Heterozygosity means one functional Aldh1l1 allele remains. The intact hap1 allele preserves normal Aldh1l1 expression. The 73 bp reference deletion at the cut site is consistent with CRISPR editing.  
Recommend: Junction PCR to confirm integration, Aldh1l1 expression check in astrocytes.

---

## 5. TG-6425 C01 — Status

- **Transgene:** Aldh1l1-EGFP
- **HiFi reads:** `/mnt/Ext_1tb_B/r84285_20260219_052427/1_C01/hifi_reads/m84285_260219_093939_s3.hifi_reads.bc2026.bam`
- **Pipeline:** `humanwgs_singleton` (PID 1216218, started 2026-05-26)
- **Output directory:** `/mnt/Ext_1tb_A/tg-integration-denovo-assembly/TG6425_C01`
- **Integration inputs:** `TG6425.tg_integration.inputs.json` (ready)
- **Post-completion command:**
  ```bash
  bash run_tg_integration.sh TG-6425 TG6425.tg_integration.inputs.json \
    /mnt/Ext_1tb_A/tg-integration-denovo-assembly/TG6425_C01_integration_analysis
  ```

---

## 6. Pending Items

| Priority | Task | Notes |
|----------|------|-------|
| High | TG-6102 original FASTA validation | Contact researcher — Addgene/GenBank/plasmid map for Aldh1l1-CreERT2 |
| High | TG-6425 C01 integration analysis | After humanwgs_singleton completes |
| Medium | TG-6283 RNA-seq follow-up | Check Camk1d expression in brain tissue |
| Medium | TG-6102 junction PCR | Confirm 73 bp deletion at chr6:90,534,206 |
| Low | TG-6903 D01 preparation | BAM at `/mnt/Ext_1tb_B/r84285_20260219_052427/1_D01/hifi_reads/m84285_260219_114241_s4.hifi_reads.bc2027.bam` (49 GB) |

---

## 7. Key Files

| File | Description |
|------|-------------|
| `workflows/wdl-common/wdl/tasks/transgene_integration.wdl` | Integration pipeline (bug fixed, commit bed8cd4) |
| `TG6425.tg_integration.inputs.json` | TG-6425 integration analysis inputs |
| `TG6425.inputs.json` | TG-6425 humanwgs_singleton inputs |
| `run_TG6425_C01.sh` | TG-6425 pipeline run script |
| `TG6102.tg_integration.extracted.inputs.json` | TG-6102 integration inputs (extracted FASTA) |
| `/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6102_A01_integration_analysis/tg_sequence_analysis/h2tg000140l_82886_85449_TG_insert.fa` | Extracted 2,563 bp actual TG insert (TG-6102) |
| `/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6102_A01_integration_analysis/tg_sequence_analysis/TG6102_sequence_comparison_report.txt` | TG FASTA mismatch analysis report |
| `/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6283_B01_integration_analysis/TG-6283_integration_report.html` | TG-6283 final HTML report |
| `/mnt/JJ_dis_8tb/tg-integration-denovo-assembly/TG6102_A01_integration_analysis/20260526_095733_transgene_integration/call-integration_report/out/report_html/TG-6102.Aldh1l1_CreERT2.integration_report.html` | TG-6102 final HTML report (extracted FASTA) |
