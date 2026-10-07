"""Fig 3d 내보내기: 조건(웰)별 joint genotype 불일치율, 3-way와 반복 쌍(pairwise).

정의는 wgs-tertiary-analysis scripts/joint_genotype_discordance.py와 같다:
  위상 무시 GT 비교, 비교 대상 = 해당 샘플 모두 유효 GT이고 적어도 하나가 non-ref인 사이트
  3-way 불일치 = 세 GT가 모두 같지 않음; pairwise 불일치 = 두 GT가 다름
차이: 상염색체(chr1–19)만 사용(요청서 규칙). 비교를 위해 전체 염색체 3-way 값도 함께 기록.
출력(tsv): OUT.wells.tsv (condition, scope, n_sites, n_discordant, discordant_frac),
           OUT.pairs.tsv (condition, sample_a, sample_b, n_sites, n_discordant, discordant_frac)
SNV_ONLY=1 환경변수를 주면 REF·ALT가 모두 1염기인 사이트(SNV)만 비교한다.
사용: [SNV_ONLY=1] python joint_discordance_pairs.py OUT_PREFIX
"""
import glob
import itertools
import os
import subprocess
import sys

BCF = "/home/ygkim/program/anaconda3/envs/wgs-tertiary-pipeline/bin/bcftools"
JOINT = {"A": "/mnt/hdd1_1tb/h2o2-wgs/joint-H2O2-A01/20260714_054055_joint",
         "B": "/mnt/hdd1_1tb/h2o2-wgs/joint-H2O2-B01/20260714_062818_joint",
         "C": "/mnt/hdd2_1tb/h2o2-wgs/joint-H2O2-C01/20260714_071318_joint",
         "D": "/mnt/hdd2_1tb/h2o2-wgs/joint-H2O2-D01/20260714_075510_joint"}
AUTOS = {f"chr{i}" for i in range(1, 20)}
out = sys.argv[1]
SNV_ONLY = os.environ.get("SNV_ONLY") == "1"


def norm(gt):
    a = gt.replace("|", "/").split("/")
    if len(a) != 2 or "." in a:
        return None
    return tuple(sorted(a))


wells, pairs = [], []
for cond, d in JOINT.items():
    vcfs = sorted(glob.glob(f"{d}/out/split_joint_small_variant_vcfs/*/*.vcf.gz"))
    names = subprocess.run(f"for v in {' '.join(vcfs)}; do {BCF} query -l $v; done", shell=True,
                           capture_output=True, text=True, check=True).stdout.split()
    sel = "| {BCF} view -v snps -m2 -M2 -Ou ".format(BCF=BCF) if SNV_ONLY else ""
    cmd = f"{BCF} merge {' '.join(vcfs)} -Ou {sel}| {BCF} query -f '%CHROM[\\t%GT]\\n'"
    tri = {"autosomes": [0, 0], "all_chromosomes": [0, 0]}
    pw = {pr: [0, 0] for pr in itertools.combinations(range(3), 2)}
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True) as p:
        for line in p.stdout:
            f = line.rstrip("\n").split("\t")
            g = [norm(x) for x in f[1:4]]
            auto = f[0] in AUTOS
            if all(x is not None for x in g) and not all(x == ("0", "0") for x in g):
                disc = len(set(g)) > 1
                tri["all_chromosomes"][0] += 1
                tri["all_chromosomes"][1] += disc
                if auto:
                    tri["autosomes"][0] += 1
                    tri["autosomes"][1] += disc
            if not auto:
                continue
            for (i, j), c in pw.items():
                if g[i] is None or g[j] is None or (g[i] == ("0", "0") and g[j] == ("0", "0")):
                    continue
                c[0] += 1
                c[1] += g[i] != g[j]
    for scope, (n, k) in tri.items():
        wells.append((cond, scope, n, k, k / n))
    for (i, j), (n, k) in pw.items():
        pairs.append((cond, names[i], names[j], n, k, k / n))
    print(cond, tri, flush=True)

with open(f"{out}.wells.tsv", "w") as fo:
    fo.write("condition\tscope\tn_sites\tn_discordant\tdiscordant_frac\n")
    for r in wells:
        fo.write("\t".join(map(str, r[:4])) + f"\t{r[4]!r}\n")
with open(f"{out}.pairs.tsv", "w") as fo:
    fo.write("condition\tsample_a\tsample_b\tn_sites\tn_discordant\tdiscordant_frac\n")
    for r in pairs:
        fo.write("\t".join(map(str, r[:5])) + f"\t{r[5]!r}\n")
