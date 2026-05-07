#!/usr/bin/env bash
# Extract CpG 5mC methylation for one sample and chromosome.
# This streams the chromosome from S3, so no full BAM is stored.
set -euo pipefail

MANIFEST=${1:?manifest.tsv}
ROW_NUM=${2:?row number, 1-based excluding header}
CHROM=${3:?chromosome, e.g. chr22}
REF=${4:?hg38.fa}
OUTDIR=${5:?output directory}
THREADS=${6:-8}
MIN_COV=${7:-5}

command -v samtools >/dev/null
command -v modkit >/dev/null
command -v bgzip >/dev/null
command -v tabix >/dev/null

LINE=$(awk -v n=$((ROW_NUM+1)) 'NR==n {print}' "$MANIFEST")
SAMPLE_ID=$(echo "$LINE" | cut -f1)
BAM_URL=$(echo "$LINE" | cut -f2)
MODS=$(echo "$LINE" | cut -f9)

[ -n "$SAMPLE_ID" ] || { echo "No sample for row $ROW_NUM" >&2; exit 1; }

SAMPLE_DIR="$OUTDIR/$CHROM/$SAMPLE_ID"
mkdir -p "$SAMPLE_DIR"
LOG="$SAMPLE_DIR/${SAMPLE_ID}.${CHROM}.modkit.log"
RAW="$SAMPLE_DIR/${SAMPLE_ID}.${CHROM}.bedmethyl.gz"
FILTERED="$SAMPLE_DIR/${SAMPLE_ID}.${CHROM}.cov${MIN_COV}.bedmethyl.gz"

if [ -s "$FILTERED" ] && [ -s "$FILTERED.tbi" ]; then
  echo "[$(date)] $SAMPLE_ID $CHROM exists; skipping"
  exit 0
fi

MOD_FLAG="5mC"
if echo "$MODS" | grep -q "5hmC"; then
  MOD_FLAG="5mCG_5hmCG"
fi

echo "[$(date)] Starting $SAMPLE_ID $CHROM" | tee "$LOG"
echo "BAM: $BAM_URL" | tee -a "$LOG"

samtools view -b -@ "$THREADS" "$BAM_URL" "$CHROM" | \
  modkit pileup \
    --cpg \
    --ref "$REF" \
    --modified-bases "$MOD_FLAG" \
    --combine-strands \
    --threads "$THREADS" \
    --bgzf \
    --log "$LOG" \
    - "$RAW"

tabix -p bed "$RAW"
bgzip -dc "$RAW" | awk -v c="$MIN_COV" '$10 >= c' | bgzip -@ 4 > "$FILTERED"
tabix -p bed "$FILTERED"

echo "Total CpGs: $(bgzip -dc "$RAW" | wc -l)" | tee -a "$LOG"
echo "CpGs >=${MIN_COV}x: $(bgzip -dc "$FILTERED" | wc -l)" | tee -a "$LOG"
echo "[$(date)] Done $SAMPLE_ID $CHROM" | tee -a "$LOG"
