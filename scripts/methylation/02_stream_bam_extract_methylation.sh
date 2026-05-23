#!/bin/bash
# =============================================================================
# Phase 2: Per-sample BAM streaming + methylation extraction
#
# Strategy (Decision R2): Download BAM to $PBS_JOBFS (local node scratch),
# run modkit pileup, save bedMethyl to /g/data, delete BAM.
# Never store the full cohort (58.4 TB) simultaneously.
#
# Submitted by: deploy/submit_pipeline.sh (which runs git pull first)
# Usage: qsub -J 1-500 scripts/00_download/02_stream_bam_extract_methylation.sh
#
# Requires:
#   PIPELINE_DIR   path to repo clone on Gadi
#   COMMIT         git commit hash (set at submission time)
# =============================================================================

#PBS -N 1kgp_methyl
#PBS -P cy94
#PBS -q normal
#PBS -l ncpus=16
#PBS -l mem=64GB
#PBS -l walltime=10:00:00
#PBS -l jobfs=250GB
#PBS -l storage=scratch/cy94+gdata/cy94+gdata/de95
#PBS -l wd
#PBS -j oe
#PBS -o /scratch/cy94/sb8857/1kgp-mqtl/logs/

set -euo pipefail

# --- Source config ---
source "${PIPELINE_DIR}/config/gadi.env"
echo "=== 1KGP methylation extraction ==="
echo "Commit:      ${COMMIT:-unknown}"
echo "Node:        $(hostname)"
echo "Array index: ${PBS_ARRAY_INDEX}"
echo "Time:        $(date)"
echo ""

# --- Load modules ---
module load samtools/1.22
module load python3/3.12.1   # required for awscli (installed via pip --user with python3/3.12.1)

# --- Get sample manifest line ---
MANIFEST="${PIPELINE_DIR}/config/bam_manifest.tsv"
SAMPLE_LINE=$(sed -n "${PBS_ARRAY_INDEX}p" "${MANIFEST}")
SAMPLE_ID=$(echo "${SAMPLE_LINE}" | cut -f1)
BAM_URL=$(echo "${SAMPLE_LINE}" | cut -f2)
BAM_BASENAME=$(basename "${BAM_URL}")
BASECALLER=$(echo "${SAMPLE_LINE}" | cut -f3)   # for logging/covariate

echo "Sample:      ${SAMPLE_ID}"
echo "Basecaller:  ${BASECALLER}"
echo "BAM:         ${BAM_BASENAME}"
echo ""

# --- Skip if already done ---
OUT_DIR="${GDATA_ROOT}/methylation"
OUT_BED="${OUT_DIR}/bedmethyl/${SAMPLE_ID}_combined.bedmethyl.gz"
OUT_HP1="${OUT_DIR}/haplotype/${SAMPLE_ID}_hp1.bedmethyl.gz"
OUT_HP2="${OUT_DIR}/haplotype/${SAMPLE_ID}_hp2.bedmethyl.gz"

if [ -f "${OUT_BED}" ] && [ -f "${OUT_HP1}" ] && [ -f "${OUT_HP2}" ]; then
    echo "[$(date)] ${SAMPLE_ID}: outputs already exist — skipping."
    exit 0
fi

# --- Stage BAM to $PBS_JOBFS (local node SSD, fast I/O) ---
BAM_LOCAL="${PBS_JOBFS}/${BAM_BASENAME}"
BAI_LOCAL="${BAM_LOCAL}.bai"

echo "[$(date)] Staging BAM to local node scratch (${PBS_JOBFS})..."
${AWS} s3 cp --no-sign-request "${BAM_URL}" "${BAM_LOCAL}"
${AWS} s3 cp --no-sign-request "${BAM_URL}.bai" "${BAI_LOCAL}"
echo "[$(date)] BAM staged: $(du -sh ${BAM_LOCAL} | cut -f1)"

# --- Run modkit pileup (Decision R4, validated 2026-02-26) ---
# Flags:
#   --modified-bases 5mC    extract 5mC only (discards 5hmC if present)
#   --cpg                   CpG context only (harmonises guppy vs dorado)
#   --combine-strands       merge +/- strand counts per position (requires --ref)
#   --phased                split into HP1 + HP2 + combined outputs
#   --min-coverage 5        filter low-coverage sites
#
# Output files:
#   <prefix>_combined.bedmethyl
#   <prefix>_1.bedmethyl    (HP1)
#   <prefix>_2.bedmethyl    (HP2)
MODKIT_OUT="${PBS_JOBFS}/modkit_out"
mkdir -p "${MODKIT_OUT}"

echo "[$(date)] Running modkit pileup (16 threads)..."
${MODKIT} pileup \
    --threads ${PBS_NCPUS} \
    --ref "${REFERENCE}" \
    --modified-bases 5mC \
    --cpg \
    --combine-strands \
    --phased \
    --prefix "${SAMPLE_ID}" \
    "${BAM_LOCAL}" \
    "${MODKIT_OUT}/"

echo "[$(date)] modkit complete. Output files:"
ls -lh "${MODKIT_OUT}/"

# --- Compress and index outputs ---
echo "[$(date)] Compressing and indexing..."
# modkit with --prefix produces: <prefix>_combined.bedmethyl, <prefix>_hp1.bedmethyl, <prefix>_hp2.bedmethyl
for SUFFIX in combined hp1 hp2; do
    SRC="${MODKIT_OUT}/${SAMPLE_ID}_${SUFFIX}.bedmethyl"
    [ -f "${SRC}" ] || { echo "WARN: missing ${SRC}"; continue; }
    bgzip -@ 4 "${SRC}"
    tabix -p bed "${SRC}.gz"
done

# --- Move to persistent storage (gdata) ---
mkdir -p "${OUT_DIR}/bedmethyl" "${OUT_DIR}/haplotype"

mv "${MODKIT_OUT}/${SAMPLE_ID}_combined.bedmethyl.gz"     "${OUT_BED}"
mv "${MODKIT_OUT}/${SAMPLE_ID}_combined.bedmethyl.gz.tbi" "${OUT_BED}.tbi"
mv "${MODKIT_OUT}/${SAMPLE_ID}_hp1.bedmethyl.gz"          "${OUT_HP1}"
mv "${MODKIT_OUT}/${SAMPLE_ID}_hp1.bedmethyl.gz.tbi"      "${OUT_HP1}.tbi"
mv "${MODKIT_OUT}/${SAMPLE_ID}_hp2.bedmethyl.gz"          "${OUT_HP2}"
mv "${MODKIT_OUT}/${SAMPLE_ID}_hp2.bedmethyl.gz.tbi"      "${OUT_HP2}.tbi"

# --- Cleanup BAM from local scratch ---
echo "[$(date)] Cleaning up staged BAM..."
rm -f "${BAM_LOCAL}" "${BAI_LOCAL}"

echo ""
echo "[$(date)] ${SAMPLE_ID}: done."
echo "  Combined: ${OUT_BED}"
echo "  HP1:      ${OUT_HP1}"
echo "  HP2:      ${OUT_HP2}"
