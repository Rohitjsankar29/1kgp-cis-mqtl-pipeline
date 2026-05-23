#!/usr/bin/env bash
# run_modkit_pileup.sh
# =====================
# Run modkit pileup on a single modBAM file to produce per-CpG bedMethyl output.
# Usage: bash run_modkit_pileup.sh <sample_id> <modbam> <ref> <outdir> [threads]
#
# Dependencies: modkit >= 0.2.0, samtools, tabix
# Author: Buckberry Lab

set -euo pipefail

SAMPLE=$1
MODBAM=$2
REF=$3
OUTDIR=$4
THREADS=${5:-16}

mkdir -p "${OUTDIR}"

LOG="${OUTDIR}/${SAMPLE}.modkit.log"
OUTFILE="${OUTDIR}/${SAMPLE}.bedmethyl.gz"

echo "[$(date)] Starting modkit pileup for ${SAMPLE}" | tee -a "${LOG}"

modkit pileup \
    --cpg \
    --ref "${REF}" \
    --modified-bases 5mC \
    --combine-strands \
    --threads "${THREADS}" \
    --bgzf \
    --log "${LOG}" \
    "${MODBAM}" \
    "${OUTFILE}"

# Index the output
tabix -p bed "${OUTFILE}"

# Coverage filter: keep only CpGs with >= 10x coverage
FILTERED="${OUTDIR}/${SAMPLE}.bedmethyl.cov10.gz"
zcat "${OUTFILE}" | awk '$10 >= 10' | bgzip > "${FILTERED}"
tabix -p bed "${FILTERED}"

echo "[$(date)] Done: ${FILTERED}" | tee -a "${LOG}"

# Summary stats
echo "[$(date)] CpG counts:" | tee -a "${LOG}"
echo "  Total CpGs: $(zcat "${OUTFILE}" | wc -l)" | tee -a "${LOG}"
echo "  CpGs >=10x: $(zcat "${FILTERED}" | wc -l)" | tee -a "${LOG}"
