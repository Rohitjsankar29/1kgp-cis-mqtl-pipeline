#!/usr/bin/env bash
set -euo pipefail
source config/gadi.env

# 1) manifest: default 2 samples
python scripts/manifest/00_build_manifest.py --out config/bam_manifest.tsv --n-samples 2

# 2) methylation extraction: on PBS use jobs/01_extract_chr_modkit.pbs instead
# qsub -J 1-2 -v CHROM=chr22,PIPELINE_DIR=$PWD jobs/01_extract_chr_modkit.pbs

# 3) after modkit jobs complete, build methylation matrix
python scripts/matrix/02_build_methylation_matrix.py \
  --manifest config/bam_manifest.tsv \
  --bedmethyl-root "$GDATA_ROOT/methylation/bedmethyl_chr" \
  --chrom chr22 \
  --out-bed "$GDATA_ROOT/methylation/matrix/pilot.chr22.methylation.Mval.bed.gz" \
  --out-beta "$GDATA_ROOT/methylation/matrix/pilot.chr22.methylation.beta.bed.gz" \
  --out-pca "$GDATA_ROOT/methylation/matrix/pilot.chr22.meth_pca.tsv" \
  --out-qc-dir "$GDATA_ROOT/methylation/qc/pilot.chr22"

# 4) covariates; add --geno-pca when genotype PCA exists
python scripts/covariates/03_build_covariates.py \
  --manifest config/bam_manifest.tsv \
  --meth-pca "$GDATA_ROOT/methylation/matrix/pilot.chr22.meth_pca.tsv" \
  --out "$GDATA_ROOT/covariates/pilot.chr22.covariates.tsv"
