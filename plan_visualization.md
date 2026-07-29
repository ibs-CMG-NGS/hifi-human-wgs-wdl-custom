# 3차 파이프라인: Transgene Integration 시각화 계획

**작성일:** 2026-05-27  
**대상 샘플:** TG-6283 B01 (chr2:5,559,779), TG-6102 A01 (chr6:90,534,206)  
**위치:** 2차 파이프라인(`transgene_integration.wdl`) 출력을 입력으로 사용하는 독립 파이프라인

---

## 1. 파이프라인 분리 근거

| 항목 | 2차 파이프라인 | 3차 파이프라인 |
|------|--------------|--------------|
| 목적 | 삽입 위치 탐색/검출 | 결과 시각화/해석 |
| 입력 | BAM, assembly FASTA | breakpoint 좌표, PAF, chimeric BAM |
| 출력 | 삽입 좌표, chimeric reads, HTML 보고서 | 논문 수준 그림 (PNG/SVG) |
| 실행 방식 | 자동화 (WDL) | 연구자가 반복 조정 |
| 도구 | bioinformatics CLI | Python / R 시각화 라이브러리 |

→ **WDL 대신 Python 스크립트 모음 + Snakemake**로 관리 (그림은 색상/범위/레이블 수작업 조정이 필수적)

---

## 2. 목표 그림 유형 (참고 논문 기준)

| Panel | 내용 | 참고 |
|-------|------|------|
| A | 트랜스진 구조 모식도 (도메인 블록) | α-Cre paper Fig.1a |
| B | 게놈 전체 삽입 위치 (karyogram) | α-Cre paper Fig.1b |
| C | 삽입 부위 시퀀스 커버리지 (hap1/hap2) | α-Cre paper Fig.1c |
| D | 컨티그 내부 구조 (host/TG 블록) | Casr paper Fig. top |
| E | 탠덤/내부 concatemer reads 시각화 | Casr paper Fig. C |
| F | **Nucleotide dot plot** (read vs 레퍼런스) | Ccr2/Ccr5 paper Fig. B, D, E |

---

## 3. Panel별 상세 계획

### Panel A — 트랜스진 구조 모식도

**목표:** 트랜스진 도메인(promoter, CDS, polyA 등)을 컬러 블록 화살표로 표현

**필요 데이터:**
- 트랜스진 도메인 좌표 (연구팀 construct map 또는 Addgene에서 확인)
  - TG-6283: gfa2 promoter / CreERT2 / SV40pA (총 4,309 bp)
  - TG-6102: Aldh1l1 5'UTR / CreERT2 / SV40pA (총 2,563 bp extracted)

**도구:**
```bash
pip install dna_features_viewer matplotlib
```

**스크립트:** `scripts/draw_construct.py`
```python
from dna_features_viewer import GraphicFeature, GraphicRecord

features = [
    GraphicFeature(start=0,    end=800,  strand=+1, color="#f9c74f", label="gfa2 promoter"),
    GraphicFeature(start=800,  end=2630, strand=+1, color="#90be6d", label="CreERT2"),
    GraphicFeature(start=2630, end=2880, strand=+1, color="#f94144", label="SV40pA"),
]
record = GraphicRecord(sequence_length=4309, features=features)
ax, _ = record.plot(figure_width=10)
```

**미결 사항:** 연구팀으로부터 도메인 정확한 bp 좌표 수령 필요

---

### Panel B — 게놈 전체 삽입 위치 (Karyogram)

**목표:** GRCm39 마우스 염색체 이디오그램에 삽입 위치 마커 표시

**필요 데이터:**
- TG-6283: chr2:5,559,779
- TG-6102: chr6:90,534,206
- GRCm39 chromosome size (UCSC 또는 Ensembl)

**도구 (R):**
```r
# BiocManager::install("karyoploteR")
library(karyoploteR)

kp <- plotKaryotype(genome="mm39")
kpPoints(kp, chr="chr2", x=5559779,  y=0.5, col="red", cex=2, pch=25)
kpPoints(kp, chr="chr6", x=90534206, y=0.5, col="red", cex=2, pch=25)
kpAddBaseNumbers(kp)
```

**스크립트:** `scripts/draw_karyogram.R`

---

### Panel C — 삽입 부위 시퀀스 커버리지

**목표:** 삽입 부위 ±300 kb 구간의 hap1/hap2 read coverage + 유전자 구조 트랙

**필요 데이터 생성:**

```bash
# 1. haplotagged BAM에서 hap별 분리
samtools view -b -d HP:1 TG-6283.GRCm39.haplotagged.bam \
  chr2:5259779-5859779 > hap1_region.bam && samtools index hap1_region.bam

samtools view -b -d HP:2 TG-6283.GRCm39.haplotagged.bam \
  chr2:5259779-5859779 > hap2_region.bam && samtools index hap2_region.bam

# 2. BigWig 변환
bamCoverage -b hap1_region.bam -o hap1.bw --binSize 100 --normalizeUsing RPKM
bamCoverage -b hap2_region.bam -o hap2.bw --binSize 100 --normalizeUsing RPKM

# 3. 유전자 exon BED (Camk1d / Aldh1l1)
grep "Camk1d" gencode.vM36.annotation.sorted.gtf | awk '$3=="exon"' | \
  awk '{print $1"\t"$4"\t"$5"\t.\t0\t"$7}' > Camk1d_exons.bed
```

**도구:**
```bash
pip install pyGenomeTracks deeptools
```

**pyGenomeTracks 설정 (`tracks_TG6283.ini`):**
```ini
[hap1]
file = hap1.bw
title = Haplotype 1
color = #4472C4
height = 3
min_value = 0

[hap2]
file = hap2.bw
title = Haplotype 2
color = #ED7D31
height = 3
min_value = 0

[TG insertion]
file = tg_insertion.bed
title = TG 삽입 위치
color = red
height = 0.5

[genes]
file = Camk1d_exons.bed
title = Camk1d
height = 2

[x-axis]
```

```bash
pyGenomeTracks --tracks tracks_TG6283.ini \
  --region chr2:5259779-5859779 \
  -o TG6283_coverage.png --dpi 300
```

**스크립트:** `scripts/generate_coverage_tracks.sh`, `scripts/tracks_TG6283.ini`, `scripts/tracks_TG6102.ini`

---

### Panel D — 컨티그 내부 구조 블록도

**목표:** chimeric contig 내 host/TG 도메인 구간을 블록으로 표시 + 아래에 HiFi reads pileup

**필요 데이터 (이미 보유):**
- `chimeric_hap1_ref_alignment_fresh.paf` → host 블록 좌표
- `tg_vs_hap1.paf` → TG 블록 좌표
- HiFi reads → chimeric contig 재정렬 (신규 생성 필요)

```bash
# HiFi reads를 chimeric contig에 재정렬
minimap2 -ax map-hifi --secondary=no chimeric_contigs_hap1.fa \
  TG-6283.hifi_reads.fastq.gz | samtools sort -o reads_vs_contig.bam
samtools index reads_vs_contig.bam
```

**도구:** `matplotlib`, `pysam`

**스크립트:** `scripts/draw_contig_structure.py`
```python
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

fig, (ax_contig, ax_reads) = plt.subplots(2, 1, figsize=(14, 6),
                                            gridspec_kw={'height_ratios':[1, 3]})

# 컨티그 블록 (PAF에서 파싱)
ax_contig.barh(0, 2400000, height=0.4, color='#b0bec5', left=0)   # host (left)
ax_contig.barh(0, 5315,    height=0.7, color='#4CAF50', left=2531427, label='TG copy 1')
ax_contig.barh(0, 403,     height=0.4, color='#9e9e9e', left=2536742)  # spacer
ax_contig.barh(0, 5316,    height=0.7, color='#2196F3', left=2537145, label='TG copy 2')
ax_contig.barh(0, 500000,  height=0.4, color='#b0bec5', left=2542461)  # host (right)

ax_contig.set_xlim(2400000, 2700000)
ax_contig.set_yticks([])
ax_contig.legend(loc='upper right')
ax_contig.set_title('h1tg000090l — chr2:5.56 Mb insertion (TG-6283)')

# reads pileup은 pysam으로 추가
plt.tight_layout()
plt.savefig('TG6283_contig_structure.png', dpi=300, bbox_inches='tight')
```

---

### Panel E — 탠덤 구조 상세 (선택)

**목표:** TG-6283의 2카피 탠덤 구조를 reads level에서 보여주는 그림 (Casr 논문 panel C 스타일)

**필요 데이터:**
- chimeric contig vs TG FASTA alignment (이미 보유: `tg_vs_hap1.paf`)
- 개별 HiFi reads 내 TG 반복 구조 (reads_vs_contig.bam)

**도구:** matplotlib custom (reads를 수평 바로 그리고, 도메인 컬러 매핑)

---

### Panel F — Nucleotide Dot Plot (read vs 레퍼런스)

**목표:** 개별 chimeric HiFi read를 x축, 연결된 레퍼런스 서열(트랜스진 + 숙주 플랭킹)을 y축으로 놓아 각 read 내부의 alignment 대각선을 그림으로 표현  
→ 탠덤 카피 수, 방향(forward/reverse), host-TG 경계, E.coli 오염 여부를 한 눈에 확인

**참고 논문 패널:** Ccr2/Ccr5 transgene paper Fig. B, D, E  
- x축: 개별 read 위치 (kb)  
- y축: 레퍼런스 스택 (E.coli backbone → TG copy 1 → TG copy 2 → host 5' flank → host 3' flank)  
- 파란 대각선 = forward 매칭, 빨간 대각선 = reverse 매칭  

**필요 데이터:**

```bash
# 1. 레퍼런스 FASTA 연결 (y축 구성)
#    TG-6283: transgene + chr2 플랭킹 ±5 kb
cat Aldh1l1-EGFP_TG-6283.fa \
    <(samtools faidx mouse_GRCm39.fasta chr2:5554779-5564779) \
    > TG6283_dotplot_ref.fa

# TG-6102: transgene + chr6 플랭킹 ±5 kb
cat Aldh1l1-CreERT2_TG-6102.fa \
    <(samtools faidx mouse_GRCm39.fasta chr6:90529206-90539206) \
    > TG6102_dotplot_ref.fa

# 2. chimeric reads FASTA 추출 (hybrid_ref_analysis에서 확보)
#    TG-6283: chimeric_hap1 contig에 정렬된 reads
samtools fasta reads_vs_contig.bam > chimeric_reads_TG6283.fa

# 3. minimap2로 PAF 생성 (read × ref dot plot용)
minimap2 -x map-hifi --secondary=no -c \
    TG6283_dotplot_ref.fa chimeric_reads_TG6283.fa \
    > TG6283_dotplot.paf
```

**도구:**
```bash
pip install matplotlib numpy pandas
# 추가 옵션 (GUI dot plot): gepard (Java), MUMmer/mummerplot
```

**스크립트:** `scripts/draw_dotplot.py`
```python
import matplotlib.pyplot as plt
import matplotlib.collections as mc
import numpy as np
import pandas as pd

def parse_paf(paf_path, min_mapq=10, min_block=200):
    cols = ["qname","qlen","qstart","qend","strand",
            "tname","tlen","tstart","tend","nmatch","alen","mapq"]
    df = pd.read_csv(paf_path, sep="\t", header=None, usecols=range(12), names=cols)
    return df[(df.mapq >= min_mapq) & (df.alen >= min_block)]

def build_ref_offsets(ref_fasta):
    """각 레퍼런스 contig의 y축 누적 시작 좌표 계산"""
    from Bio import SeqIO
    offsets, pos = {}, 0
    for rec in SeqIO.parse(ref_fasta, "fasta"):
        offsets[rec.id] = pos
        pos += len(rec.seq)
    return offsets, pos

def draw_dotplot(paf_path, ref_fasta, out_png, top_reads=30):
    df = parse_paf(paf_path)
    ref_offsets, total_ref = build_ref_offsets(ref_fasta)

    # 가장 긴 chimeric reads 선택
    top = (df.groupby("qname")["alen"].sum()
             .nlargest(top_reads).index)
    df = df[df.qname.isin(top)]

    # reads를 x축 기준 정렬 → 행 배치
    read_order = (df.groupby("qname")["qlen"].first()
                    .sort_values(ascending=False).index)
    read_ypos = {name: i for i, name in enumerate(read_order)}
    read_ymax = df.groupby("qname")["qlen"].first()

    fig, ax = plt.subplots(figsize=(16, 10))
    fwd_segs, rev_segs = [], []

    for _, row in df.iterrows():
        rx = read_ypos[row.qname]          # read 행 인덱스 → y panel
        ref_y_off = ref_offsets.get(row.tname, 0)

        # 각 read panel의 x 범위: [rx*scale, (rx+1)*scale]
        # 실제로는 각 read를 별도 subplot으로 그리는 방식
        seg_x = [row.qstart, row.qend]
        seg_y = [ref_y_off + row.tstart, ref_y_off + row.tend]

        if row.strand == "+":
            fwd_segs.append(list(zip(seg_x, seg_y)))
        else:
            rev_segs.append(list(zip([row.qstart, row.qend],
                                     [ref_y_off + row.tend, ref_y_off + row.tstart])))

    ax.add_collection(mc.LineCollection(fwd_segs, colors="#4472C4", linewidths=0.8, alpha=0.8))
    ax.add_collection(mc.LineCollection(rev_segs, colors="#C0392B", linewidths=0.8, alpha=0.8))

    # y축 레퍼런스 구간 구분선 + 레이블
    for name, off in ref_offsets.items():
        ax.axhline(off, color="gray", lw=0.5, ls="--")
        ax.text(-500, off + 200, name, fontsize=7, va="bottom")

    ax.set_xlabel("Read position (bp)")
    ax.set_ylabel("Reference position (bp)")
    ax.set_title(f"Nucleotide dot plot — chimeric reads vs reference")
    ax.autoscale()
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"Saved: {out_png}")

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--paf",  required=True)
    p.add_argument("--ref",  required=True, help="연결된 레퍼런스 FASTA")
    p.add_argument("--out",  required=True)
    p.add_argument("--top",  type=int, default=30)
    args = p.parse_args()
    draw_dotplot(args.paf, args.ref, args.out, args.top)
```

**실행 예시:**
```bash
python scripts/draw_dotplot.py \
    --paf  output/TG6283/TG6283_dotplot.paf \
    --ref  output/TG6283/TG6283_dotplot_ref.fa \
    --out  output/TG6283/TG6283_dotplot.png \
    --top  50
```

**TG-6283 기대 결과:** 각 chimeric read에서 TG copy 1 → TG copy 2 대각선이 두 번 연속 나타남 (탠덤 2카피 확인)  
**TG-6102 기대 결과:** TG 대각선이 1회만 나타남 (단일 카피)  

**미결 사항:**
- `reads_vs_contig.bam` 생성 필요 (Panel D 전처리와 공유)
- 레퍼런스 FASTA 경계를 논문 수준으로 정교하게 조정 (E.coli backbone 포함 여부 검토)

---

## 4. 파일 구조

```
workflows/
└── visualization/                  ← 3차 파이프라인 디렉토리
    ├── README.md
    ├── Snakefile                   ← (선택) 전체 워크플로우 자동화
    ├── scripts/
    │   ├── draw_construct.py       ← Panel A
    │   ├── draw_karyogram.R        ← Panel B
    │   ├── generate_coverage_tracks.sh  ← Panel C 전처리
    │   ├── tracks_TG6283.ini       ← Panel C pyGenomeTracks 설정
    │   ├── tracks_TG6102.ini
    │   ├── draw_contig_structure.py ← Panel D
    │   ├── draw_dotplot.py         ← Panel F
    │   └── parse_paf.py            ← PAF 파싱 공통 유틸
    └── output/
        ├── TG6283/
        └── TG6102/
```

---

## 5. 환경 설정

```bash
# Python 환경 (hifi-human-wgs conda에 추가 또는 별도 환경)
pip install dna_features_viewer pyGenomeTracks deeptools pysam matplotlib biopython

# R 패키지
Rscript -e 'BiocManager::install("karyoploteR")'
```

---

## 6. 실행 우선순위

| 순서 | 작업 | 소요 시간 | 선행 조건 |
|------|------|----------|----------|
| 1 | IGV 스크린샷 (빠른 확인) | 30분 | haplotagged BAM |
| 2 | Panel C — pyGenomeTracks 커버리지 | 2~3시간 | deeptools 설치 |
| 3 | Panel D — 컨티그 구조 블록도 | 3~4시간 | PAF 파일 (보유) |
| 4 | Panel B — karyoploteR | 1시간 | breakpoint 좌표 (보유) |
| 5 | Panel A — 트랜스진 모식도 | 1시간 | 연구팀 도메인 좌표 수령 후 |
| 6 | Panel E — 탠덤 reads 시각화 | 4~6시간 | Panel D 완료 후 |
| 7 | Panel F — Nucleotide dot plot | 2~3시간 | reads_vs_contig.bam (Panel D와 공유) |

---

## 7. 미결 사항

- [ ] 연구팀으로부터 TG-6283 / TG-6102 construct 도메인 정확한 bp 좌표 수령 (Panel A)
- [ ] TG-6425 C01 파이프라인 완료 후 동일 시각화 적용
- [ ] TG-6103 D01 분석 완료 후 4샘플 비교 그림 추가 검토
