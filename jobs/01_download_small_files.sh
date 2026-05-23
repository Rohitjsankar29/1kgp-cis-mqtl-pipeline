#!/bin/bash
# =============================================================================
# Phase 1 Download: Small files (~75 GB)
# Metadata, QC stats, per-sample VCFs (SNV + SV), GVCFs
#
# Run from: gadi-dm.nci.org.au (data-mover node — external network access)
# Submit as: qsub -q copyq 01_download_small_files.sh
#
# Set GDATA_ROOT before submitting:
#   export GDATA_ROOT=/g/data/<project>/1kgp-mqtl
# =============================================================================

#PBS -N 1kgp_download_small
#PBS -P cy94
#PBS -q copyq
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l walltime=12:00:00
#PBS -l jobfs=10GB
#PBS -l storage=gdata/PROJ                # Replace PROJ with project code
#PBS -l wd
#PBS -o logs/download_small.o
#PBS -e logs/download_small.e

set -euo pipefail

S3_BASE="https://s3.amazonaws.com/1000g-ont/PROCESSED_DATA/ALIGNED_TO_HG38"
GDATA_ROOT="${GDATA_ROOT:-/g/data/PROJ/1kgp-mqtl}"  # Replace PROJ

echo "[$(date)] Starting Phase 1 download to ${GDATA_ROOT}"

# --- Create directory structure ---
mkdir -p \
    "${GDATA_ROOT}/metadata" \
    "${GDATA_ROOT}/qc/cramino" \
    "${GDATA_ROOT}/qc/hp_dp_stats" \
    "${GDATA_ROOT}/variants/snv/per_sample_phased" \
    "${GDATA_ROOT}/variants/snv/gvcf" \
    "${GDATA_ROOT}/variants/snv/joint" \
    "${GDATA_ROOT}/variants/sv/per_sample" \
    "${GDATA_ROOT}/variants/sv/joint" \
    "${GDATA_ROOT}/methylation/bedmethyl" \
    "${GDATA_ROOT}/methylation/haplotype" \
    "${GDATA_ROOT}/references/hg38"

# --- Metadata ---
echo "[$(date)] Downloading metadata..."
curl -sL "https://s3.amazonaws.com/1000g-ont/PROCESSED_DATA/1kGP_LRSC_500_ONT_Metadata.tsv" \
    -o "${GDATA_ROOT}/metadata/1kGP_LRSC_500_ONT_Metadata.tsv"
curl -sL "https://s3.amazonaws.com/1000g-ont/PROCESSED_DATA/READ_ME.txt" \
    -o "${GDATA_ROOT}/metadata/READ_ME.txt"

# --- QC: cramino alignment stats ---
echo "[$(date)] Downloading cramino QC stats..."
SAMPLE_IDS=$(awk 'NR>1 {print $1}' "${GDATA_ROOT}/metadata/1kGP_LRSC_500_ONT_Metadata.tsv")

for SAMPLE in ${SAMPLE_IDS}; do
    URL="${S3_BASE}/CRAMINO/${SAMPLE}_aligned_hg38_cramino.txt"
    OUT="${GDATA_ROOT}/qc/cramino/${SAMPLE}_cramino.txt"
    [ -f "${OUT}" ] && continue
    curl -sfL "${URL}" -o "${OUT}" 2>/dev/null || echo "WARN: cramino missing for ${SAMPLE}"
done
echo "[$(date)] cramino: done"

# --- SNV phased VCFs (Clair3) ---
echo "[$(date)] Downloading Clair3 phased VCFs..."
# Generate URL list from BAM filenames (basename determines VCF name)
aws s3 ls --no-sign-request s3://1000g-ont/PROCESSED_DATA/ALIGNED_TO_HG38/CLAIR3/PHASED_VCF/ \
    | awk '{print $4}' | grep '\.vcf\.gz$' > /tmp/phased_vcf_list.txt

while read -r FNAME; do
    OUT="${GDATA_ROOT}/variants/snv/per_sample_phased/${FNAME}"
    [ -f "${OUT}" ] && continue
    aws s3 cp --no-sign-request \
        "s3://1000g-ont/PROCESSED_DATA/ALIGNED_TO_HG38/CLAIR3/PHASED_VCF/${FNAME}" \
        "${OUT}"
done < /tmp/phased_vcf_list.txt
echo "[$(date)] Clair3 phased VCFs: done"

# --- SNV GVCFs (Clair3) ---
echo "[$(date)] Downloading Clair3 GVCFs..."
aws s3 ls --no-sign-request s3://1000g-ont/PROCESSED_DATA/ALIGNED_TO_HG38/CLAIR3/GVCF/ \
    | awk '{print $4}' | grep -v '^$' > /tmp/gvcf_list.txt

while read -r FNAME; do
    OUT="${GDATA_ROOT}/variants/snv/gvcf/${FNAME}"
    [ -f "${OUT}" ] && continue
    aws s3 cp --no-sign-request \
        "s3://1000g-ont/PROCESSED_DATA/ALIGNED_TO_HG38/CLAIR3/GVCF/${FNAME}" \
        "${OUT}"
done < /tmp/gvcf_list.txt
echo "[$(date)] Clair3 GVCFs: done"

# --- SV VCFs (Sniffles2) ---
echo "[$(date)] Downloading Sniffles2 SV VCFs..."
aws s3 ls --no-sign-request s3://1000g-ont/PROCESSED_DATA/ALIGNED_TO_HG38/SNIFFLES_v2.6.2/ \
    | awk '{print $4}' | grep '\.vcf\.gz$' > /tmp/sv_vcf_list.txt

while read -r FNAME; do
    OUT="${GDATA_ROOT}/variants/sv/per_sample/${FNAME}"
    [ -f "${OUT}" ] && continue
    aws s3 cp --no-sign-request \
        "s3://1000g-ont/PROCESSED_DATA/ALIGNED_TO_HG38/SNIFFLES_v2.6.2/${FNAME}" \
        "${OUT}"
done < /tmp/sv_vcf_list.txt
echo "[$(date)] Sniffles2 SV VCFs: done"

echo "[$(date)] Phase 1 download complete."
echo "Storage used: $(du -sh ${GDATA_ROOT} | cut -f1)"
