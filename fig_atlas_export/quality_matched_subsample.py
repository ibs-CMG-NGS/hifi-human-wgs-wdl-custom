"""품질 맞춘 BAM: read 품질(rq → Q) 구간별 유지 확률로 read를 추려, 12샘플이 같은 품질 분포·같은 depth가 되게 한다.

유지 확률표(keep_prob_bins.tsv)는 chr1:50–70 Mb의 primary read 길이 가중 Q 분포와 샘플 평균 depth로
구간별 목표 depth = 12샘플 최솟값을 정해 만든 것(목표 합계 22.8x).
선택은 read 이름의 CRC32로 결정해, 같은 read의 primary/supplementary가 함께 남거나 빠진다.
stdin: samtools view -h SAM, stdout: SAM
사용: samtools view -h BAM chr1 chr2 chr3 | python quality_matched_subsample.py KEEP_TSV SAMPLE | samtools view -b ...
"""
import math
import re
import sys
import zlib

import pandas as pd

keep = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
row = keep.loc[sys.argv[2]]
bins = []
for col, p in row.items():
    lo, hi = map(float, re.findall(r"[\d.]+", col))
    bins.append((lo, hi, float(p)))

n_in = n_out = 0
out = sys.stdout
for line in sys.stdin:
    if line.startswith("@"):
        out.write(line)
        continue
    n_in += 1
    f = line.split("\t", 11)
    m = re.search(r"\trq:f:([0-9.eE+-]+)", line)
    q = -10 * math.log10(max(1 - float(m.group(1)), 1e-6)) if m else 20.0
    p = 0.0
    for lo, hi, pp in bins:
        if lo <= q < hi:
            p = pp
            break
    if zlib.crc32(f[0].encode()) / 2**32 < p:
        out.write(line)
        n_out += 1
print(f"{sys.argv[2]}: kept {n_out}/{n_in} records", file=sys.stderr)
