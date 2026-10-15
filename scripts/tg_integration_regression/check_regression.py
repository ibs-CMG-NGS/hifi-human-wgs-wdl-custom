#!/usr/bin/env python3
"""Compare transgene_integration outputs against truth.tsv.

Usage: check_regression.py <sample> <miniwdl run dir (the one holding outputs.json)> [tolerance_bp=20]
A truth junction passes when a read junction cluster (read_junctions_tsv, or the breakpoint_range of older
outputs) lies within tolerance on the same side and chromosome. Exit 1 if any truth junction is missing.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sample, run = sys.argv[1], sys.argv[2]
tol = int(sys.argv[3]) if len(sys.argv) > 3 else 20
truth = [l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, "truth.tsv")) if l.strip() and not l.startswith("#")]
truth = [t for t in truth if t[0] == sample]
out = json.load(open(os.path.join(run, "outputs.json")))
get = lambda k: out.get(f"transgene_integration.{k}")

found = []
jt = get("read_junctions_tsv")
if jt and os.path.exists(jt):
    for l in list(open(jt))[1:]:
        side, chrom, pos, n = l.split("\t")
        found.append((side, chrom, int(pos), int(n)))
if get("stat_breakpoint_range"):   # strongest left/right pair, even when one side has a single read
    import re
    m = re.match(r"(\S+?):([\d,]+)-([\d,]+)", get("stat_breakpoint_range"))
    if m:
        found += [("left", m.group(1), int(m.group(2).replace(",", "")), 0),
                  ("right", m.group(1), int(m.group(3).replace(",", "")), 0)]
print(f"{sample}: primary {get('stat_breakpoint_chr')}:{get('stat_breakpoint_pos')}  range {get('stat_breakpoint_range')}")
print(f"  de novo primary {get('stat_integration_chr')}:{get('stat_integration_pos')}")
print(f"  QC warnings: {get('stat_qc_warnings')}")
ok = True
for _, _, _, side, chrom, pos, *note in truth:
    hit = [f for f in found if f[0] == side and f[1] == chrom and abs(f[2] - int(pos)) <= tol]
    print(f"  {'PASS' if hit else 'FAIL'} {side:5s} {chrom}:{pos}" + (f"  ({hit[0][3]} reads at {hit[0][2]})" if hit else ""))
    ok &= bool(hit)
sys.exit(0 if ok else 1)
