#!/bin/bash
# 00_prepare_superpop_labels.sh -- one-off preparation of superpopulation
# labels for the 452 samples. RUN THIS ONCE INTERACTIVELY on a Gadi login
# node (login nodes have internet; compute nodes don't).
#
# What it does:
#   1) Downloads the 1000 Genomes integrated call-set sample panel file
#      (~30 KB), which maps sample -> population -> superpopulation.
#   2) Joins that against the 452-sample metadata to produce
#      config/metadata_with_pop.tsv, adding a `superpopulation` column.
#   3) Reports the per-superpopulation counts so you can sanity-check
#      against the expected AFR 122 / SAS 101 / EAS 93 / EUR 71 / AMR 65.
#
# No SU cost -- runs on the login node, seconds to complete.
#
# Usage:
#   bash 00_prepare_superpop_labels.sh
#   # verify counts match the expected superpop breakdown
#   # then submit ancestry_eur.pbs

set -euo pipefail

CONFIG=/scratch/cy94/rs4477/1kgp-cis-mqtl-pipeline/config
IN_META=$CONFIG/metadata_452.tsv       # sample sex platform basecaller (no pop)
OUT_META=$CONFIG/metadata_with_pop.tsv
PANEL_URL=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
PANEL_LOCAL=$CONFIG/1000g_sample_panel.tsv

# 1) fetch the 1000G panel if not already present
if [[ ! -s "$PANEL_LOCAL" ]]; then
    echo "[$(date +%H:%M:%S)] downloading 1000G panel file"
    curl -fSL "$PANEL_URL" -o "$PANEL_LOCAL"
fi
echo "[$(date +%H:%M:%S)] panel file: $(wc -l < $PANEL_LOCAL) lines"
head -3 "$PANEL_LOCAL"
echo

# 2) join by sample_id; the panel's columns are: sample pop super_pop gender
#    the metadata's first column is `sample`
python3 - <<PY
import pandas as pd, sys
meta = pd.read_csv("$IN_META", sep="\t")
if "sample" not in meta.columns:
    sys.exit(f"expected column 'sample' in metadata; have {list(meta.columns)}")

panel = pd.read_csv("$PANEL_LOCAL", sep="\t")
# panel columns: sample, pop, super_pop, gender  (super_pop = superpopulation)
panel = panel.rename(columns={"sample": "sample", "super_pop": "superpopulation", "pop": "population"})

merged = meta.merge(panel[["sample", "population", "superpopulation"]], on="sample", how="left")
missing = merged["superpopulation"].isna().sum()
if missing:
    print(f"WARNING: {missing} samples had no superpopulation match", file=sys.stderr)
    print(merged.loc[merged["superpopulation"].isna(), "sample"].head(10).to_list(), file=sys.stderr)

merged.to_csv("$OUT_META", sep="\t", index=False)
print(f"wrote $OUT_META  ({len(merged)} rows)")
print()
print("counts by superpopulation:")
print(merged["superpopulation"].value_counts().to_string())
PY

echo
echo "[$(date +%H:%M:%S)] DONE"
echo "Expected: AFR 122  SAS 101  EAS 93  EUR 71  AMR 65"
