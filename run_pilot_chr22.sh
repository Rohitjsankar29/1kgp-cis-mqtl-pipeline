#!/usr/bin/env bash
set -euo pipefail

# Tiny pilot helper. On Gadi, use PBS for methylation extraction.
# This script documents the run order and runs only non-PBS stages.

source config/gadi.env
mkdir -p logs "$GDATA_ROOT/methylation/bedmethyl" "$GDATA_ROOT/methylation/matrix" "$GDATA_ROOT/covariates" "$GDATA_ROOT/tensorqtl" "$GDATA_ROOT/qc"

echo "Step 1: build BAM manifest"
python scripts/manifest/00_build_manifest.py

echo "Step 2: submit methylation extraction manually, starting with one sample"
echo "qsub -J 1-1 scripts/methylation/02_stream_bam_extract_methylation.sh"
echo "After it succeeds, try: qsub -J 1-2 scripts/methylation/02_stream_bam_extract_methylation.sh"

echo "Step 3: after bedMethyl files exist, aggregate matrix"
echo "python scripts/matrix/01_aggregate_bedmethyl.py \\"
echo "  --input-dir $GDATA_ROOT/methylation/bedmethyl \\"
echo "  --pattern '*_combined.bedmethyl.gz' \\"
echo "  --output-dir $GDATA_ROOT/methylation/matrix \\"
echo "  --min-coverage 5 --min-samples 0.8"

echo "Step 4: build covariates after methylation/genotype PCs exist"
echo "python scripts/covariates/02_prepare_covariates.py --help"

echo "Step 5: run TensorQTL"
echo "python scripts/tensorqtl/05_run_tensorqtl_pilot.py --help"
