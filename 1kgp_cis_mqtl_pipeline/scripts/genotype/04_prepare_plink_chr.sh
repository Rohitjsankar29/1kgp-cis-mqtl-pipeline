#!/usr/bin/env bash
# Merge per-sample phased VCFs into PLINK files for a pilot chromosome.
# Edit VCF_DIR to where small VCFs were downloaded.
set -euo pipefail
MANIFEST=${1:?manifest.tsv}
VCF_DIR=${2:?per-sample phased VCF dir}
CHROM=${3:?chr22}
OUT_PREFIX=${4:?output plink prefix}
THREADS=${5:-4}
mkdir -p "$(dirname "$OUT_PREFIX")"
SAMPLES=$(mktemp)
awk 'NR>1 {print $1}' "$MANIFEST" > "$SAMPLES"
# This assumes VCFs are already bgzipped/indexed and contain the selected samples.
# For full production, use bcftools merge by chromosome before plink2.
bcftools merge --threads "$THREADS" -r "$CHROM" -S "$SAMPLES" "$VCF_DIR"/*.vcf.gz -Oz -o "${OUT_PREFIX}.${CHROM}.vcf.gz"
tabix -p vcf "${OUT_PREFIX}.${CHROM}.vcf.gz"
plink2 --vcf "${OUT_PREFIX}.${CHROM}.vcf.gz" --set-all-var-ids '@:#:$r:$a' --make-bed --out "$OUT_PREFIX"
plink2 --bfile "$OUT_PREFIX" --pca approx 10 --out "${OUT_PREFIX}.pca"
