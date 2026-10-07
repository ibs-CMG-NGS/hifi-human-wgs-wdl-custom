"""Fig 2e: 샘플 고유(private) SNV의 치환 스펙트럼 — 산화 손상(8-oxoG → C>A) 흔적 확인.

입력: 12샘플 DeepVariant VCF를 bcftools merge한 표(stdin):
      CHROM POS REF ALT [GT AD DP] × 12   (PASS + RefCall, 상염색체, biallelic SNV)
정의:
  called(s)   : GT에 alt 포함 (PASS 호출)
  evidence(s) : AD의 alt read ≥ 1 (RefCall 포함). 레코드가 없으면 DeepVariant 후보 조건
                (alt read ≥2, ≥12%) 미달 → alt 증거 없음으로 간주
  private(s)  : s만 called, 다른 11개 샘플 모두 evidence 없음, s의 alt read ≥ 3, DP 15–100
                → 여러 개체 혼합 배양에서 독립 웰 간 공유될 수 없는 체세포 변이 후보
  shared_all  : 12샘플 모두 called → 생식세포 다형성 기준 스펙트럼
출력: OUT.private_sites.tsv.gz (사이트별), OUT.spectrum.tsv (샘플×클래스×VAF 구간 카운트)
사용: bcftools ... | python mutation_spectrum.py samples.tsv OUT_PREFIX
"""
import gzip
import sys

import numpy as np
import pandas as pd

FASTA = "/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/mouse_GRCm39.fasta"
meta = pd.read_csv(sys.argv[1], sep="\t")
OUT = sys.argv[2]
SIDS = meta["sample_id"].tolist()
NS = len(SIDS)
COMP = str.maketrans("ACGT", "TGCA")

fai = pd.read_csv(f"{FASTA}.fai", sep="\t", header=None, usecols=[0, 1, 2, 3, 4],
                  names=["chrom", "len", "off", "lb", "lw"]).set_index("chrom")


def load_seq(chrom):
    r = fai.loc[chrom]
    nlines = (r.len + r.lb - 1) // r.lb
    with open(FASTA, "rb") as f:
        f.seek(r.off)
        raw = f.read(int(r.len + nlines * (r.lw - r.lb)))
    return raw.replace(b"\n", b"").replace(b"\r", b"").upper().decode()


def classify(ref, alt, ctx):
    """피리미딘 기준 6-class와 96-context (5'-X[REF>ALT]Y-3')."""
    if ref in "GA":
        ref, alt, ctx = ref.translate(COMP), alt.translate(COMP), ctx.translate(COMP)[::-1]
    return f"{ref}>{alt}", f"{ctx[0]}[{ref}>{alt}]{ctx[2]}"


def vaf_bin(v):
    return "VAF<0.2" if v < 0.2 else ("VAF0.2-0.35" if v < 0.35 else "VAF>=0.35")


counts = {}  # (sample, set, vafbin, class6) -> n
ctx96 = {}   # (sample, set, ctx96) -> n
cur_chrom, seq = None, None
fo = gzip.open(f"{OUT}.private_sites.tsv.gz", "wt")
fo.write("chrom\tpos\tref\talt\tsample\tgroup\tclass6\tctx96\tvaf\tdp\talt_reads\n")
grp = dict(zip(SIDS, meta["group"]))

for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    chrom, pos, ref, alt = f[0], int(f[1]), f[2], f[3]
    if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
        continue
    called, evid, info = [], [], {}
    for i in range(NS):
        gt, ad, dp = f[4 + 3 * i: 7 + 3 * i]
        if gt in (".", "./.", ".|."):
            continue
        ads = ad.split(",")
        a = int(ads[1]) if len(ads) > 1 and ads[1] != "." else 0
        r = int(ads[0]) if ads[0] != "." else 0
        if a >= 1:
            evid.append(i)
        if "1" in gt:
            called.append(i)
            info[i] = (a, r, int(dp) if dp != "." else a + r)
    if not called:
        continue
    if chrom != cur_chrom:
        seq, cur_chrom = load_seq(chrom), chrom
    ctx = seq[pos - 2: pos + 1]
    if len(ctx) != 3 or "N" in ctx or ctx[1] != ref:
        continue
    c6, c96 = classify(ref, alt, ctx)
    if len(called) == NS:
        for i in range(NS):
            a, r, dp = info[i]
            k = (SIDS[i], "shared_all", vaf_bin(a / max(a + r, 1)), c6)
            counts[k] = counts.get(k, 0) + 1
            ctx96[(SIDS[i], "shared_all", c96)] = ctx96.get((SIDS[i], "shared_all", c96), 0) + 1
    elif len(called) == 1 and evid == called:
        i = called[0]
        a, r, dp = info[i]
        if a >= 3 and 15 <= dp <= 100:
            v = a / max(a + r, 1)
            k = (SIDS[i], "private", vaf_bin(v), c6)
            counts[k] = counts.get(k, 0) + 1
            ctx96[(SIDS[i], "private", c96)] = ctx96.get((SIDS[i], "private", c96), 0) + 1
            fo.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{SIDS[i]}\t{grp[SIDS[i]]}\t{c6}\t{c96}\t{v:.3f}\t{dp}\t{a}\n")
fo.close()

rows = [dict(sample=s, set=st, vaf_bin=vb, class6=c, n=n) for (s, st, vb, c), n in counts.items()]
spec = pd.DataFrame(rows)
spec["group"] = spec["sample"].map(grp)
spec.to_csv(f"{OUT}.spectrum.tsv", sep="\t", index=False)
pd.DataFrame([dict(sample=s, set=st, ctx96=c, n=n) for (s, st, c), n in ctx96.items()]).to_csv(
    f"{OUT}.ctx96.tsv", sep="\t", index=False)

# 요약 출력
pd.set_option("display.width", 250)
for st in ("private", "shared_all"):
    d = spec[spec["set"] == st].groupby(["sample", "class6"])["n"].sum().unstack(fill_value=0)
    d = d.loc[[s for s in SIDS if s in d.index]]
    tot = d.sum(axis=1)
    frac = d.div(tot, axis=0)
    frac.insert(0, "n_total", tot)
    frac.insert(0, "group", [grp[s] for s in frac.index])
    print(f"\n## {st} SNV spectrum (fraction)\n" + frac.to_string(float_format=lambda x: f"{x:.3f}"))
d = spec[spec["set"] == "private"].groupby(["sample", "vaf_bin"])["n"].sum().unstack(fill_value=0)
print("\n## private SNV by VAF bin\n" + d.loc[[s for s in SIDS if s in d.index]].to_string())
