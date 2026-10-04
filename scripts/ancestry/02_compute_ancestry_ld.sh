#!/bin/bash
set -uo pipefail
export MAMBA_EXE=/home/571/rs4477/.local/bin/micromamba
export MAMBA_ROOT_PREFIX=/scratch/cy94/rs4477/micromamba
eval "$($MAMBA_EXE shell hook --shell bash --root-prefix $MAMBA_ROOT_PREFIX)"
micromamba activate plink

G=/g/data/cy94/rs4477/downstream/genome_452
IN_DIR=/scratch/cy94/rs4477/coloc/ldref
KEEP_DIR=$G/ancestry/keep
OUT_DIR=$G/ancestry/plink
mkdir -p "$OUT_DIR"

POPS="${@:-EUR}"
for POP in $POPS; do
    KEEP="$KEEP_DIR/keep_${POP}.txt"
    [ -s "$KEEP" ] || { echo "MISSING: $KEEP"; continue; }
    N=$(wc -l < "$KEEP")
    mkdir -p "$OUT_DIR/$POP"
    for CHR in $(seq 1 22); do
        OUT_BFILE=$OUT_DIR/$POP/chr${CHR}
        [ -f "${OUT_BFILE}.bed" ] && { echo "skip $POP chr$CHR (exists)"; continue; }
        echo "[$(date +%H:%M:%S)] $POP  chr$CHR  (n=$N)"
        plink2 --bfile $IN_DIR/chr${CHR} --keep "$KEEP" \
               --maf 0.05 --geno 0.05 --make-bed --out "$OUT_BFILE" --silent
    done
done
echo "[$(date +%H:%M:%S)] DONE"
