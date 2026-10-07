"""품질 맞춘 재호출(chr1–3) vs 원래 호출(같은 chr1–3)의 클론성·C>A 지표 비교.

dataset = original (원 파이프라인 phased VCF / joint VCF) | matched (품질·depth 맞춘 BAM의 DeepVariant / GLnexus)
샘플 지표 (PASS, chr1–3):
  het_snv, het_del, het_ins, het_hom_peak_frac_snv (SNV VAF > 0.9), het_snv_vaf_frac_lt0.2
  private SNV (12샘플 merge, 다른 11개에 alt read 없음, alt ≥ 3, DP 15–100; mutation_spectrum.py와 같은 규칙):
    n_private, n_private_CA, CA_frac, CA_per_1M_het_snv, nonCA_per_1M_het_snv
조건 지표: joint VCF SNV 3-way/pairwise 불일치율 (joint_discordance_pairs.py와 같은 정의)
검정: 각 샘플 지표에 대해 condition + run 선형모형 F (로그 변환 없이), run 안 정확 permutation(13,824)
출력: OUT.samples.tsv, OUT.pairs.tsv, OUT.summary.txt
사용: python qm_metrics.py QM_DIR OUT_PREFIX
"""
import itertools
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import stats

BCF = "/home/ygkim/program/anaconda3/envs/wgs-tertiary-pipeline/bin/bcftools"
REPO = "/data_4tb/hifi-human-wgs-wdl-custom"
CHR = "chr1,chr2,chr3"
QM, OUT = sys.argv[1], sys.argv[2]
meta = pd.read_csv(f"{REPO}/base_mod_per_sample/samples.tsv", sep="\t")
SIDS = meta.sample_id.tolist()
GRP = dict(zip(meta.sample_id, meta.group))
RUN = dict(zip(meta.sample_id, meta.run))
JOINT = {"A": "/mnt/hdd1_1tb/h2o2-wgs/joint-H2O2-A01/20260714_054055_joint",
         "B": "/mnt/hdd1_1tb/h2o2-wgs/joint-H2O2-B01/20260714_062818_joint",
         "C": "/mnt/hdd2_1tb/h2o2-wgs/joint-H2O2-C01/20260714_071318_joint",
         "D": "/mnt/hdd2_1tb/h2o2-wgs/joint-H2O2-D01/20260714_075510_joint"}


def vcf_of(ds, sid):
    if ds == "matched":
        return f"{QM}/dv/{sid}/{sid}.vcf.gz"
    o = meta.set_index("sample_id").loc[sid, "wgs_haplotagged_bam"].split("/out/")[0] + "/out"
    return f"{o}/phased_small_variant_vcf/{sid}.GRCm39.small_variants.phased.vcf.gz"


def sample_metrics(ds, sid):
    cmd = [BCF, "query", "-t", CHR, "-i", 'FILTER="PASS"', "-f", "%REF\t%ALT[\t%GT\t%VAF]\n", vcf_of(ds, sid)]
    snv_v, ndel, nins = [], 0, 0
    for line in subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.splitlines():
        ref, alt, gt, vaf = line.split("\t")
        al = gt.replace("|", "/").split("/")
        if len(al) != 2 or al[0] == al[1] or "." in al:
            continue
        alts = alt.split(",")
        if len(ref) == 1 and all(len(a) == 1 for a in alts):
            try:
                snv_v.append(max(float(x) for x in vaf.split(",") if x != "."))
            except ValueError:
                pass
        elif len(ref) > max(len(a) for a in alts):
            ndel += 1
        else:
            nins += 1
    v = np.array(snv_v)
    return dict(het_snv=len(v), het_del=ndel, het_ins=nins, het_hom_peak_frac_snv=float((v > 0.9).mean()),
                het_snv_vaf_frac_lt0_2=float((v < 0.2).mean()))


def private(ds):
    files = [vcf_of(ds, s) for s in SIDS]
    view = " ".join(f"<({BCF} view -t {CHR} -v snps -m2 -M2 -Ou {f})" for f in files)
    cmd = f"{BCF} merge -m none --force-samples {view} -Ou | {BCF} query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT\\t%AD\\t%DP]\\n'"
    cnt = {s: [0, 0] for s in SIDS}
    p = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        ref, alt = f[2], f[3]
        called, evid, info = [], [], {}
        for i in range(len(SIDS)):
            gt, ad, dp = f[4 + 3 * i: 7 + 3 * i]
            if gt in (".", "./.", ".|."):
                continue
            ads = ad.split(",")
            a = int(ads[1]) if len(ads) > 1 and ads[1] != "." else 0
            if a >= 1:
                evid.append(i)
            if "1" in gt:
                called.append(i)
                info[i] = (a, int(dp) if dp != "." else a)
        if len(called) == 1 and evid == called:
            a, dp = info[called[0]]
            if a >= 3 and 15 <= dp <= 100:
                ca = (ref, alt) in (("C", "A"), ("G", "T"))
                cnt[SIDS[called[0]]][0] += 1
                cnt[SIDS[called[0]]][1] += ca
    p.wait()
    return cnt


def discordance(ds):
    rows = []
    for g, sub in meta.groupby("group"):
        if ds == "matched":
            cmd = f"{BCF} view -v snps -m2 -M2 -Ou {QM}/joint/{g}.vcf.gz | {BCF} query -f '[%GT\\t]\\n'"
            names = subprocess.run([BCF, "query", "-l", f"{QM}/joint/{g}.vcf.gz"], capture_output=True, text=True).stdout.split()
        else:
            import glob
            vcfs = sorted(glob.glob(f"{JOINT[g]}/out/split_joint_small_variant_vcfs/*/*.vcf.gz"))
            names = [subprocess.run([BCF, "query", "-l", v], capture_output=True, text=True).stdout.strip() for v in vcfs]
            cmd = f"{BCF} merge -r {CHR} {' '.join(vcfs)} -Ou | {BCF} view -v snps -m2 -M2 -Ou | {BCF} query -f '[%GT\\t]\\n'"
        pw = {pr: [0, 0] for pr in itertools.combinations(range(3), 2)}
        for line in subprocess.run(["bash", "-c", cmd], capture_output=True, text=True).stdout.splitlines():
            gts = [x for x in line.split("\t") if x]
            if len(gts) != 3:
                continue
            n = []
            for x in gts:
                a = x.replace("|", "/").split("/")
                n.append(None if len(a) != 2 or "." in a else tuple(sorted(a)))
            for (i, j), c in pw.items():
                if n[i] is None or n[j] is None or (n[i] == ("0", "0") and n[j] == ("0", "0")):
                    continue
                c[0] += 1
                c[1] += n[i] != n[j]
        for (i, j), (nn, k) in pw.items():
            rows.append(dict(dataset=ds, condition=g, sample_a=names[i], sample_b=names[j], n_sites=nn,
                             n_discordant=k, discordant_frac=k / nn if nn else np.nan))
    return rows


def _ssr(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    return float(r @ r)


def exact_perm_F(y):
    """condition + run 선형모형의 condition F와 run 안 정확 permutation p (13,824 배치)."""
    runs_lab = np.array([RUN[s] for s in SIDS])
    R = np.column_stack([np.ones(len(y))] + [(runs_lab == r).astype(float) for r in (2, 3)])
    ssr_red = _ssr(R, y)
    def F(gs):
        G = np.column_stack([(np.array(gs) == g).astype(float) for g in "BCD"])
        ssr = _ssr(np.column_stack([R, G]), y)
        return ((ssr_red - ssr) / 3) / (ssr / (len(y) - 6))
    g0 = [GRP[s] for s in SIDS]
    obs = F(g0)
    runs = [[i for i, s in enumerate(SIDS) if RUN[s] == r] for r in (1, 2, 3)]
    ge = tot = 0
    for a in itertools.permutations(runs[0]):
        for b in itertools.permutations(runs[1]):
            for c in itertools.permutations(runs[2]):
                gs = list(g0)
                for src, dst in ((runs[0], a), (runs[1], b), (runs[2], c)):
                    for i, j in zip(src, dst):
                        gs[i] = g0[j]
                tot += 1
                ge += F(gs) >= obs - 1e-9
    return obs, ge / tot


rows, prows = [], []
for ds in ("original", "matched"):
    pv = private(ds)
    for s in SIDS:
        m = sample_metrics(ds, s)
        n, ca = pv[s]
        m.update(dataset=ds, sample_id=s, condition=GRP[s], run=RUN[s], n_private=n, n_private_CA=ca,
                 CA_frac=ca / n if n else np.nan, CA_per_1M_het_snv=ca / m["het_snv"] * 1e6,
                 nonCA_per_1M_het_snv=(n - ca) / m["het_snv"] * 1e6)
        rows.append(m)
        print(ds, s, m, flush=True)
    prows += discordance(ds)
d = pd.DataFrame(rows)
d.to_csv(f"{OUT}.samples.tsv", sep="\t", index=False)
pdf = pd.DataFrame(prows)
pdf.to_csv(f"{OUT}.pairs.tsv", sep="\t", index=False)

lines = []
for ds in ("original", "matched"):
    sub = d[d.dataset == ds].set_index("sample_id").loc[SIDS]
    lines.append(f"## {ds} (chr1-3)")
    for m in ("het_snv", "het_del", "het_ins", "het_hom_peak_frac_snv", "het_snv_vaf_frac_lt0_2", "n_private",
              "CA_frac", "CA_per_1M_het_snv", "nonCA_per_1M_het_snv"):
        means = sub.groupby("condition")[m].mean()
        order = {r: "".join(g.sort_values(m).condition) for r, g in sub.groupby("run")}
        Fv, pp = exact_perm_F(sub[m].to_numpy(float))
        lines.append(f"  {m}: means " + " ".join(f"{g}={means[g]:.5g}" for g in "ABCD")
                     + f" | order low->high per run {order} | F={Fv:.2f} exact perm p={pp:.4f}")
    pp_ = pdf[pdf.dataset == ds].groupby("condition").discordant_frac.agg(["mean", "min", "max"])
    lines.append("  SNV pairwise discordance: " + " ".join(f"{g}={r['mean']:.4f}[{r['min']:.4f}-{r['max']:.4f}]" for g, r in pp_.iterrows()))
open(f"{OUT}.summary.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
