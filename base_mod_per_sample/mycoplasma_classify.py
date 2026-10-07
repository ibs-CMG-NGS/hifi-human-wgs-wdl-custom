"""mycoplasma_check.sh 보조: minimap2 SAM(stdin)에서 read별 분류 + 6mA 밀도 기록, 또는 요약.

read 행: read, length, target(mouse_nuclear/mouse_chrM/종 이름/unmapped), mapq, gc, n_A_sites, n_6mA
- 6mA: MM의 A+a / T-a 호출 중 ML>=179(0.7) 개수. implicit 모드라 목록에 없는 A는 비수정.
- 분모: A+a는 read(원래 방향)의 A 수, T-a는 T 수.
사용: ... | python mycoplasma_classify.py OUT.reads.tsv.gz
      python mycoplasma_classify.py --summary OUTDIR
"""
import glob
import gzip
import os
import sys

ML_MIN = 179
COMP = str.maketrans("ACGTN", "TGCAN")


def species_of(rname, desc):
    if rname in desc:
        return desc[rname]
    return "mouse_chrM" if rname == "chrM" else "mouse_nuclear"


def classify(out):
    desc = {}
    for line in open("/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/mycoplasma_panel.fa"):
        if line.startswith(">"):
            acc, _, rest = line[1:].strip().partition(" ")
            desc[acc] = "_".join(rest.split()[:2])
    with gzip.open(out + ".tmp", "wt") as fo:
        fo.write("read\tlength\ttarget\tmapq\tgc\tn_A_sites\tn_6mA\n")
        for line in sys.stdin:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            flag = int(f[1])
            if flag & 0x900:
                continue
            seq = f[9]
            if flag & 16:
                seq = seq.translate(COMP)[::-1]
            target = "unmapped" if flag & 4 else species_of(f[2], desc)
            tags = {t[:2]: t[5:] for t in f[11:] if t[:2] in ("MM", "ML")}
            n_sites = n_mod = 0
            if "MM" in tags and "ML" in tags:
                ml = [int(x) for x in tags["ML"].split(",")[1:]] if tags["ML"] else []
                i = 0
                for grp in tags["MM"].rstrip(";").split(";"):
                    parts = grp.split(",")
                    code = parts[0].rstrip(".?")
                    k = (len(parts) - 1) * (len(code) - 2)
                    if code in ("A+a", "T-a"):
                        n_sites += seq.count(code[0])
                        n_mod += sum(1 for x in ml[i:i + k] if x >= ML_MIN)
                    i += k
            gc = (seq.count("G") + seq.count("C")) / max(len(seq), 1)
            fo.write(f"{f[0]}\t{len(seq)}\t{target}\t{f[4]}\t{gc:.3f}\t{n_sites}\t{n_mod}\n")
    os.replace(out + ".tmp", out)


def summary(outdir):
    import pandas as pd
    rows = []
    for p in sorted(glob.glob(f"{outdir}/*.reads.tsv.gz")):
        sid = os.path.basename(p).replace(".reads.tsv.gz", "")
        d = pd.read_csv(p, sep="\t")
        tot_bp = d["length"].sum()
        for t, g in d.groupby("target"):
            rows.append({"sample": sid, "target": t, "reads": len(g), "frac_reads": len(g) / len(d),
                         "frac_bp": g["length"].sum() / tot_bp, "median_gc": g["gc"].median(),
                         "6mA_per_A": g["n_6mA"].sum() / max(g["n_A_sites"].sum(), 1),
                         "frac_reads_6mA>1%": (g["n_6mA"] / g["n_A_sites"].clip(lower=1) > 0.01).mean()})
        u = d[d["target"] == "unmapped"]
        rows.append({"sample": sid, "target": "unmapped_lowGC(<0.35)", "reads": int((u["gc"] < 0.35).sum()),
                     "frac_reads": (u["gc"] < 0.35).sum() / len(d)})
    s = pd.DataFrame(rows)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 500)
    print(s.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    s.to_csv(f"{outdir}/mycoplasma_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    if sys.argv[1] == "--summary":
        summary(sys.argv[2])
    else:
        classify(sys.argv[1])
