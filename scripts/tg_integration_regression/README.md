# transgene_integration regression

Truth set: `truth.tsv` (junctions verified at read level, 2026-10-14/15).

```bash
# full run (TG-6102: construct_type=BAC, bac_region=chr6:90478421-90706261)
miniwdl run workflows/transgene_integration.wdl --input <sample>.inputs.json --dir <out>
python3 scripts/tg_integration_regression/check_regression.py <sample> <out>/_LAST

# quick check of one task's python block on saved PAFs
python3 scripts/tg_integration_regression/run_task_python.py extract_integration_coords <workdir> \
  hap1_tg_paf=... hap2_tg_paf=... hap1_ref_paf=... hap2_ref_paf=... min_match_bp=500 min_flank_bp=2000 \
  construct_type=plasmid bac_region= transgene_fasta=<construct.fa>
```

2026-10-15 run (`/data_3tb/shared/tg-integration-denovo-assembly/wdl_regression_20261015/run2`): all 6 truth junctions PASS.
TG-6102 used a mini BAM (saved chr6 region reads re-aligned) because the full BAM is not available.
