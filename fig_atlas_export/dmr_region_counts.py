"""Fig 4e DMR 입력: 상염색체 영역별·샘플별 5mC 집계(메틸화 read 수, 전체 read 수, CpG 수).

입력(MOD 인자)
  5mC  (기본): 원 파이프라인 pb-CpG-tools combined bed(모델 기반 5mC %, coverage). 메틸화 read 수 = mod% × cov / 100.
  5hmC       : base_modification_per_sample.wdl의 modkit CpG pileup TSV(chr, 0-based pos, N_valid, N_mod; jasmine 26.1.3)
영역 집합 (상염색체만; 성염색체는 성 구성 차이로 제외)
  promoter : GENCODE vM36 protein_coding·lncRNA 유전자 TSS -1000 ~ +500
  genebody : protein_coding 유전자 전체
  cgi      : UCSC cpgIslandExt
  shore    : CGI 양쪽 2 kb (좌우 별도 영역)
  tile1k   : 1 kb 고정 타일
영역 합은 샘플마다 자기 CpG의 누적합을 searchsorted로 구하므로 영역이 겹쳐도 된다.
영역 유지: 모든 샘플에서 CpG ≥ 5개, 영역 coverage 합 ≥ 50
출력: OUT_DIR/{set}.meth.tsv.gz, {set}.cov.tsv.gz (행=영역, 열=샘플), {set}.regions.tsv.gz
사용: python dmr_region_counts.py samples.tsv OUT_DIR [5mC|5hmC]
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd

RES = "/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39"
meta = pd.read_csv(sys.argv[1], sep="\t")
OUT = sys.argv[2]
MOD = sys.argv[3] if len(sys.argv) > 3 else "5mC"
HMC_DIR = "/mnt/OIM_HDD2/base_mod_discovery_pilot_output/h2o2_per_sample/5hmC"
os.makedirs(OUT, exist_ok=True)
AUTO = [f"chr{i}" for i in range(1, 20)]
MIN_CPG, MIN_COV = 5, 50


def cpg_path(wgs_bam):
    run_dir = wgs_bam.split("/out/")[0]
    sid = wgs_bam.split("/")[-1].replace(".GRCm39.haplotagged.bam", "")
    return f"{run_dir}/out/cpg_combined_bed/{sid}.GRCm39.cpg_pileup.combined.bed.gz"


def load_sample(r):
    """상염색체 CpG별 (chrom, begin, meth read 수, cov)"""
    if MOD == "5hmC":
        df = pd.read_csv(f"{HMC_DIR}/{r.sample_id}/_LAST/call-modbam_cpg_pileup/work/{r.sample_id}.5hmC.tsv.gz",
                         sep="\t", header=None, names=["chrom", "begin", "cov", "meth"], dtype={"chrom": str})
        df = df[df.chrom.isin(AUTO)]
        return df.assign(cov=df["cov"].astype(float), meth=df["meth"].astype(float))
    df = pd.read_csv(cpg_path(r.wgs_haplotagged_bam), sep="\t", comment="#", usecols=[0, 1, 3, 5],
                     names=["chrom", "begin", "mod", "cov"], header=None, dtype={"chrom": str},
                     compression="gzip", low_memory=False)
    df = df[df.chrom.isin(AUTO)]
    df["cov"] = pd.to_numeric(df["cov"], errors="coerce")
    df = df.dropna()
    return df.assign(meth=df["mod"].astype(float) * df["cov"] / 100.0)


# 영역 정의
genes = []
with gzip.open(f"{RES}/gencode.vM36.annotation.gtf.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        x = line.split("\t", 9)
        if x[2] != "gene" or x[0] not in AUTO:
            continue
        a = x[8]
        gtype = a.split('gene_type "')[1].split('"')[0]
        if gtype not in ("protein_coding", "lncRNA"):
            continue
        name = a.split('gene_name "')[1].split('"')[0]
        genes.append((x[0], int(x[3]) - 1, int(x[4]), x[6], name, gtype))
g = pd.DataFrame(genes, columns=["chrom", "start", "end", "strand", "name", "gtype"])
tss = np.where(g.strand == "+", g.start, g.end - 1)
prom = pd.DataFrame(dict(chrom=g.chrom, start=np.where(g.strand == "+", tss - 1000, tss - 500),
                         end=np.where(g.strand == "+", tss + 500, tss + 1000), name=g.name, gtype=g.gtype))
body = g[g.gtype == "protein_coding"][["chrom", "start", "end", "name", "gtype"]]
cgi = pd.read_csv(f"{RES}/cpgIslandExt.txt.gz", sep="\t", header=None, usecols=[1, 2, 3])
cgi.columns = ["chrom", "start", "end"]
cgi = cgi[cgi.chrom.isin(AUTO)].reset_index(drop=True)
cgi["name"] = cgi.chrom + ":" + cgi.start.astype(str)
cgi["gtype"] = "cgi"
shore = pd.concat([cgi.assign(end=cgi.start, start=cgi.start - 2000, name=cgi.name + "_L", gtype="shore"),
                   cgi.assign(start=cgi.end, end=cgi.end + 2000, name=cgi.name + "_R", gtype="shore")])
fai = pd.read_csv(f"{RES}/mouse_GRCm39.fasta.fai", sep="\t", header=None, usecols=[0, 1], names=["c", "n"]).set_index("c")
tiles = pd.concat([pd.DataFrame(dict(chrom=c, start=np.arange(0, int(fai.loc[c, "n"]), 1000)))
                   for c in AUTO]).assign(end=lambda d: d.start + 1000)
tiles["name"] = tiles.chrom + ":" + tiles.start.astype(str)
tiles["gtype"] = "tile"
SETS = {"promoter": prom, "genebody": body, "cgi": cgi, "shore": shore, "tile1k": tiles}
for k in SETS:
    SETS[k] = SETS[k].assign(start=lambda d: d.start.clip(lower=0)).reset_index(drop=True)

res = {k: {"meth": {}, "cov": {}, "ncpg": {}} for k in SETS}
for _, r in meta.iterrows():
    df = load_sample(r)
    for k, R in SETS.items():
        m_out = np.zeros(len(R)); c_out = np.zeros(len(R)); n_out = np.zeros(len(R))
        for c, sub in df.groupby("chrom"):
            o = np.argsort(sub.begin.to_numpy())
            pos = sub.begin.to_numpy()[o]
            cov = sub["cov"].to_numpy(float)[o]
            meth = sub["meth"].to_numpy(float)[o]
            cm, cc = np.concatenate([[0], np.cumsum(meth)]), np.concatenate([[0], np.cumsum(cov)])
            idx = np.flatnonzero(R.chrom.to_numpy() == c)
            lo = np.searchsorted(pos, R.start.to_numpy()[idx])
            hi = np.searchsorted(pos, R.end.to_numpy()[idx])
            m_out[idx], c_out[idx], n_out[idx] = cm[hi] - cm[lo], cc[hi] - cc[lo], hi - lo
        res[k]["meth"][r.sample_id], res[k]["cov"][r.sample_id], res[k]["ncpg"][r.sample_id] = m_out, c_out, n_out
    print(f"loaded {r.sample_id}", flush=True)

for k, R in SETS.items():
    M, C, N = (pd.DataFrame(res[k][x]) for x in ("meth", "cov", "ncpg"))
    keep = (N.min(axis=1) >= MIN_CPG) & (C.min(axis=1) >= MIN_COV)
    R2 = R[keep.to_numpy()].copy()
    R2["ncpg"] = N[keep].median(axis=1).to_numpy()
    R2["region_id"] = [f"{k}_{i}" for i in range(len(R2))]
    R2.to_csv(f"{OUT}/{k}.regions.tsv.gz", sep="\t", index=False)
    M[keep].round(2).set_axis(R2.region_id).to_csv(f"{OUT}/{k}.meth.tsv.gz", sep="\t")
    C[keep].set_axis(R2.region_id).to_csv(f"{OUT}/{k}.cov.tsv.gz", sep="\t")
    print(f"{k}: {keep.sum()} / {len(R)} regions kept", flush=True)
