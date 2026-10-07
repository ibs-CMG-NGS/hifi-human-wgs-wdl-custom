"""6mA per-sample pileup(run_6mA_pileup.sh 출력)을 유전체 특징별로 집계.

입력 행: chrom, pos(0-based), out_strand, N_valid, N_mod — A 또는 T 참조 위치마다 +/- 두 행.
실제 6mA 호출은 (참조 A, +)와 (참조 T, -) 행에만 있다(정방향/역방향 가닥의 A).
(A, -)와 (T, +) 행은 커버리지만 있고 수정 호출이 사실상 0이라(1e-7 수준) 포함하면
비율이 절반으로 희석되므로 제외하고, qc_strand 층에만 기록한다.

출력(tsv.gz): layer, feature, n_sites, n_hit_sites(N_mod>0), N_valid, N_mod
layer:
  chrom          염색체 전체 (chrM 포함, 미배치 contig 제외)
  qc_strand      참조 염기(A/T) × 출력 strand(+/-) 균형 점검
  element        promoter(TSS±1kb) > exon > intron > intergenic (GENCODE vM36)
  cgi            CpG island 안/밖
  rep_class / rep_family / rep_name   RepeatMasker (UCSC mm39), 비반복은 none
  gene_body / gene_promoter           protein_coding + lncRNA 유전자 (gene_id|name|type)
  bin100k        100 kb 고정 구간
  ctx5           A 중심 5-mer (A가 있는 가닥 기준)
사용: python aggregate_6mA.py SAMPLE.6mA.tsv.gz OUT.tsv.gz [max_rows]
"""
import gzip
import sys
import time

import numpy as np
import pandas as pd

RES = "/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39"
FASTA = f"{RES}/mouse_GRCm39.fasta"
GTF = f"{RES}/gencode.vM36.annotation.gtf.gz"
RMSK = f"{RES}/rmsk.mm39.txt.gz"
CGI = f"{RES}/cpgIslandExt.txt.gz"
CHROMS = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY", "chrM"]
PROM = 1000
BIN = 100_000
GENE_TYPES = {"protein_coding", "lncRNA"}
CHUNK = 20_000_000

src, out = sys.argv[1], sys.argv[2]
max_rows = int(sys.argv[3]) if len(sys.argv) > 3 else None
t0 = time.time()


def log(msg):
    print(f"[{time.time() - t0:7.0f}s] {msg}", flush=True)


# ---------- 참조 서열 ----------
fai = pd.read_csv(f"{FASTA}.fai", sep="\t", header=None, usecols=[0, 1, 2, 3, 4],
                  names=["chrom", "len", "off", "lb", "lw"]).set_index("chrom")
CHROM_LEN = {c: int(fai.loc[c, "len"]) for c in CHROMS}
LUT = np.full(256, 4, dtype=np.uint8)  # A0 C1 G2 T3 N4
for i, b in enumerate(b"ACGT"):
    LUT[b] = i
    LUT[b + 32] = i


def load_seq(chrom):
    r = fai.loc[chrom]
    nlines = (r.len + r.lb - 1) // r.lb
    with open(FASTA, "rb") as f:
        f.seek(r.off)
        raw = f.read(int(r.len + nlines * (r.lw - r.lb)))
    raw = np.frombuffer(raw, dtype=np.uint8)
    return LUT[raw[(raw != 10) & (raw != 13)]][: r.len]


# ---------- 주석 ----------
def intervals_by_chrom(df):
    return {c: g for c, g in df.groupby("chrom")}


log("loading GTF")
gtf = pd.read_csv(GTF, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                  names=["chrom", "type", "start", "end", "strand", "attr"])
gtf = gtf[gtf["type"].isin(["gene", "exon"]) & gtf["chrom"].isin(CHROMS)]
gtf["start"] -= 1  # 0-based half-open
genes = gtf[gtf["type"] == "gene"].copy()
genes["gene_id"] = genes["attr"].str.extract(r'gene_id "([^"]+)"', expand=False)
genes["gene_name"] = genes["attr"].str.extract(r'gene_name "([^"]+)"', expand=False)
genes["gene_type"] = genes["attr"].str.extract(r'gene_type "([^"]+)"', expand=False)
exons = gtf[gtf["type"] == "exon"][["chrom", "start", "end"]]
genes_all = genes[["chrom", "start", "end", "strand"]]
genes = genes[genes["gene_type"].isin(GENE_TYPES)].reset_index(drop=True)
genes["label"] = genes["gene_id"] + "|" + genes["gene_name"] + "|" + genes["gene_type"]
GENE_LABELS = ["none"] + genes["label"].tolist()
genes["idx"] = np.arange(1, len(genes) + 1)
assert len(GENE_LABELS) < 65535
tss_all = genes_all.copy()
tss_all["tss"] = np.where(tss_all["strand"] == "+", tss_all["start"], tss_all["end"] - 1)

log("loading RepeatMasker")
rmsk = pd.read_csv(RMSK, sep="\t", header=None, usecols=[5, 6, 7, 10, 11, 12],
                   names=["chrom", "start", "end", "name", "cls", "fam"])
rmsk = rmsk[rmsk["chrom"].isin(CHROMS)]
REP = {}
for col in ("cls", "fam", "name"):
    labels = ["none"] + sorted(rmsk[col].unique())
    REP[col] = labels
    rmsk[col + "_i"] = pd.Categorical(rmsk[col], categories=labels).codes
assert len(REP["name"]) < 65535

cgi = pd.read_csv(CGI, sep="\t", header=None, usecols=[1, 2, 3], names=["chrom", "start", "end"])
cgi = cgi[cgi["chrom"].isin(CHROMS)]

ELEMENTS = ["intergenic", "intron", "exon", "promoter"]
G_BY = intervals_by_chrom(genes)
GA_BY = intervals_by_chrom(genes_all)
TSS_BY = intervals_by_chrom(tss_all)
EX_BY = intervals_by_chrom(exons)
R_BY = intervals_by_chrom(rmsk)
C_BY = intervals_by_chrom(cgi)

BIN_OFF, BIN_LABELS = {}, []
for c in CHROMS:
    BIN_OFF[c] = len(BIN_LABELS)
    BIN_LABELS += [f"{c}:{b * BIN}" for b in range((CHROM_LEN[c] + BIN - 1) // BIN)]
CTX_LABELS = ["".join("ACGT"[(k >> (2 * (3 - j))) & 3] for j in range(4)) for k in range(256)]
CTX_LABELS = [f"{s[:2]}A{s[2:]}" for s in CTX_LABELS] + ["N"]


def paint(arr, df, val=None, col=None):
    if df is None:
        return
    s, e = df["start"].to_numpy(), df["end"].to_numpy()
    v = df[col].to_numpy() if col else None
    for i in range(len(s)):
        arr[s[i]:e[i]] = val if col is None else v[i]


def build_chrom(c):
    n = CHROM_LEN[c]
    seq = load_seq(c)
    elem = np.zeros(n, np.uint8)
    paint(elem, GA_BY.get(c), 1)
    paint(elem, EX_BY.get(c), 2)
    t = TSS_BY.get(c)
    if t is not None:
        prom = pd.DataFrame({"start": np.clip(t["tss"] - PROM, 0, n), "end": np.clip(t["tss"] + PROM + 1, 0, n)})
        paint(elem, prom, 3)
    cg = np.zeros(n, np.uint8)
    paint(cg, C_BY.get(c), 1)
    lab = {"element": elem, "cgi": cg}
    r = R_BY.get(c)
    for col, dt in (("cls", np.uint8), ("fam", np.uint8), ("name", np.uint16)):
        a = np.zeros(n, dt)
        paint(a, r, col=col + "_i")
        lab["rep_" + {"cls": "class", "fam": "family", "name": "name"}[col]] = a
    gb = np.zeros(n, np.uint16)
    gp = np.zeros(n, np.uint16)
    g = G_BY.get(c)
    if g is not None:
        g = g.assign(len=g["end"] - g["start"]).sort_values("len", ascending=False)  # 겹치면 짧은 유전자 우선
        paint(gb, g, col="idx")
        tss = np.where(g["strand"] == "+", g["start"], g["end"] - 1)
        gpdf = pd.DataFrame({"start": np.clip(tss - PROM, 0, n), "end": np.clip(tss + PROM + 1, 0, n),
                             "idx": g["idx"].to_numpy()})
        paint(gp, gpdf, col="idx")
    lab["gene_body"], lab["gene_promoter"] = gb, gp
    return seq, lab


LAYERS = {
    "element": ELEMENTS, "cgi": ["nonCGI", "CGI"],
    "rep_class": REP["cls"], "rep_family": REP["fam"], "rep_name": REP["name"],
    "gene_body": GENE_LABELS, "gene_promoter": GENE_LABELS,
    "bin100k": BIN_LABELS, "ctx5": CTX_LABELS, "chrom": CHROMS,
    "qc_strand": ["A+", "A-", "T+", "T-"],
}
acc = {L: np.zeros((len(v), 4), np.float64) for L, v in LAYERS.items()}  # n_sites, n_hit, N_valid, N_mod


def add(layer, idx, nv, nm):
    K = len(LAYERS[layer])
    acc[layer][:, 0] += np.bincount(idx, minlength=K)
    acc[layer][:, 1] += np.bincount(idx, weights=(nm > 0).astype(float), minlength=K)
    acc[layer][:, 2] += np.bincount(idx, weights=nv, minlength=K)
    acc[layer][:, 3] += np.bincount(idx, weights=nm, minlength=K)


cur, seq, lab = None, None, None
carry = None
nrows = 0


def process(df):
    """같은 염색체의 행들. 위치별로 +/- 행을 합산한 뒤 집계."""
    global cur, seq, lab
    c = df["chrom"].iat[0]
    if c not in CHROM_LEN:
        return
    if c != cur:
        log(f"building annotation arrays for {c}")
        seq, lab = build_chrom(c)
        cur = c
    pos = df["pos"].to_numpy(np.int64)
    st = df["strand"].to_numpy()
    nv = df["nv"].to_numpy(np.float64)
    nm = df["nm"].to_numpy(np.float64)
    base = seq[pos]
    # qc: 출력 strand별 균형
    q = np.where(base == 0, 0, 2) + (st == "-").astype(np.int64)
    ok = (base == 0) | (base == 3)
    add("qc_strand", q[ok], nv[ok], nm[ok])
    # A가 있는 가닥의 행만 사용 (위치당 1행)
    keep = ((base == 0) & (st == "+")) | ((base == 3) & (st == "-"))
    upos, nv, nm, base = pos[keep], nv[keep], nm[keep], base[keep]
    add("chrom", np.full(len(upos), CHROMS.index(c)), nv, nm)
    for L in ("element", "cgi", "rep_class", "rep_family", "rep_name", "gene_body", "gene_promoter"):
        add(L, lab[L][upos].astype(np.int64), nv, nm)
    add("bin100k", BIN_OFF[c] + upos // BIN, nv, nm)
    # 5-mer 문맥 (A 가닥 기준)
    n = len(seq)
    offs = np.arange(-2, 3)
    idx = np.clip(upos[:, None] + offs[None, :], 0, n - 1)
    k = seq[idx].astype(np.int64)  # (m,5)
    rev = base == 3
    k[rev] = (3 - k[rev][:, ::-1]) % 8  # 역상보 (N=4 → 3-4=-1 → 7, N 판정 유지)
    bad = (k > 3).any(axis=1)
    flank = k[:, [0, 1, 3, 4]]
    code = (flank[:, 0] << 6) | (flank[:, 1] << 4) | (flank[:, 2] << 2) | flank[:, 3]
    code[bad] = 256
    add("ctx5", code, nv, nm)


reader = pd.read_csv(src, sep="\t", header=None, names=["chrom", "pos", "strand", "nv", "nm"],
                     dtype={"chrom": str, "pos": np.int64, "strand": str, "nv": np.int32, "nm": np.int32},
                     compression="gzip", chunksize=CHUNK, nrows=max_rows)
for chunk in reader:
    nrows += len(chunk)
    if carry is not None:
        chunk = pd.concat([carry, chunk], ignore_index=True)
    # 청크 경계에서 마지막 위치는 다음 청크와 합쳐야 하므로 보류
    last_c, last_p = chunk["chrom"].iat[-1], chunk["pos"].iat[-1]
    tail = (chunk["chrom"] == last_c) & (chunk["pos"] == last_p)
    carry = chunk[tail]
    chunk = chunk[~tail]
    for c, g in chunk.groupby("chrom", sort=False):
        process(g)
    log(f"{nrows:,} rows read (at {last_c}:{last_p:,})")
if carry is not None and len(carry):
    process(carry)

rows = []
for L, labels in LAYERS.items():
    a = acc[L]
    for i in np.nonzero(a[:, 0])[0]:
        rows.append((L, labels[i], int(a[i, 0]), int(a[i, 1]), int(a[i, 2]), int(a[i, 3])))
pd.DataFrame(rows, columns=["layer", "feature", "n_sites", "n_hit_sites", "N_valid", "N_mod"]).to_csv(
    out, sep="\t", index=False, compression="gzip")
log(f"wrote {out} ({len(rows):,} rows) from {nrows:,} input rows")
