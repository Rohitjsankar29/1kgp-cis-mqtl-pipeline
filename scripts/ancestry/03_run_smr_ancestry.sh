#!/bin/bash
# 03_run_smr_ancestry.sh -- ancestry-restricted SMR re-run.
# Matches the invocation pattern in scripts/coloc/run_smr_disease.sh
# (per-chromosome, --thread-num 2, --diff-freq-prop 0.99), pointing at the
# ancestry-restricted PLINK panels from step 02.
#
# Cost: ~1 s per probe x ~1000 probes per chr x 22 chr x 3 diseases = a few
# minutes wall-time per ancestry. ~2-3 SU per ancestry on normalq.
#
# Usage:
#   POPS="EUR"                  # primary sensitivity
#   POPS="EUR AFR SAS EAS AMR"  # stretch
#   ./03_run_smr_ancestry.sh $POPS

set -euo pipefail
G=/g/data/cy94/rs4477/downstream/genome_452
SMR=/scratch/cy94/rs4477/tools/squashfs-root/usr/bin/smr
export LD_LIBRARY_PATH=/scratch/cy94/rs4477/tools/squashfs-root/usr/lib:${LD_LIBRARY_PATH:-}
BESD_DIR=/scratch/cy94/rs4477/coloc/besd
GWAS_DIR=/scratch/cy94/rs4477/coloc/gwas

# Trait code -> GWAS .ma filename (matches your existing convention)
declare -A GWAS_MA=(
    [t2d]=t2d_hg38.ma
    [cvd]=cad_hg38.ma
    [eGFR]=eGFR_hg38.ma
)

POPS="${@:-EUR}"

for POP in $POPS; do
    PANEL_DIR=$G/ancestry/plink/$POP
    OUT_DIR=$G/ancestry/smr/$POP
    mkdir -p "$OUT_DIR"

    if [[ ! -f "$PANEL_DIR/chr22.bed" ]]; then
        echo "MISSING: $PANEL_DIR/chr22.bed  --  run 02_compute_ancestry_ld.sh first"
        continue
    fi

    for DISEASE in "${!GWAS_MA[@]}"; do
        MA=$GWAS_DIR/${GWAS_MA[$DISEASE]}
        [[ -s "$MA" ]] || { echo "MISSING: $MA"; continue; }

        echo "[$(date +%H:%M:%S)] $POP  ·  $DISEASE"
        cd "$OUT_DIR"

        for N in $(seq 1 22); do
            [[ -s smr_${DISEASE}_chr${N}.smr ]] && continue
            [[ -f $PANEL_DIR/chr${N}.bed ]] || continue

            $SMR --bfile          $PANEL_DIR/chr${N} \
                 --gwas-summary   $MA \
                 --beqtl-summary  $BESD_DIR/mqtl_chr${N} \
                 --diff-freq-prop 0.99 \
                 --thread-num     2 \
                 --out            smr_${DISEASE}_chr${N} \
                 2>&1 | grep -iE "saved|error" | tail -1
        done

        # concatenate per-chr outputs -> smr_<DISEASE>_all.smr
        if [[ -s smr_${DISEASE}_chr22.smr ]]; then
            head -1 smr_${DISEASE}_chr22.smr > smr_${DISEASE}_all.smr
            for N in $(seq 1 22); do
                tail -n +2 smr_${DISEASE}_chr${N}.smr 2>/dev/null >> smr_${DISEASE}_all.smr
            done
            TOT=$(($(wc -l < smr_${DISEASE}_all.smr)-1))
            SMRSIG=$(awk -v t=$TOT -F'\t' 'NR>1 && $19!="NA" && $19+0<0.05/t' \
                     smr_${DISEASE}_all.smr | wc -l)
            HEIDIPASS=$(awk -v t=$TOT -F'\t' 'NR>1 && $19!="NA" && $19+0<0.05/t \
                             && $20!="NA" && $20+0>0.05' \
                        smr_${DISEASE}_all.smr | wc -l)
            echo "$POP $DISEASE: $TOT tested, SMR-sig $SMRSIG, HEIDI-pass $HEIDIPASS"
        fi
    done
done

echo "[$(date +%H:%M:%S)] DONE -- per-ancestry SMR outputs under $G/ancestry/smr/"
