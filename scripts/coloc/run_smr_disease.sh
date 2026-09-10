#!/bin/bash
DISEASE=$1
MA=$2
SMR=/scratch/cy94/rs4477/tools/smr-1.4.2-linux-x86_64/squashfs-root/usr/bin/smr
cd /scratch/cy94/rs4477/coloc
for N in $(seq 1 22); do
  [ -s smr_${DISEASE}_chr$N.smr ] && continue
  $SMR --bfile ldref/chr$N --gwas-summary $MA --beqtl-summary besd/mqtl_chr$N --out smr_${DISEASE}_chr$N --thread-num 2 --diff-freq-prop 0.99 2>&1 | grep -iE "saved|error" | tail -1
done
head -1 smr_${DISEASE}_chr22.smr > smr_${DISEASE}_all.smr
for N in $(seq 1 22); do tail -n +2 smr_${DISEASE}_chr$N.smr 2>/dev/null >> smr_${DISEASE}_all.smr; done
TOT=$(($(wc -l < smr_${DISEASE}_all.smr)-1))
echo "$DISEASE: $TOT CpGs tested"
echo "SMR-sig: $(awk -v t=$TOT -F'\t' 'NR>1 && $19!="NA" && $19+0<0.05/t' smr_${DISEASE}_all.smr|wc -l)"
echo "SMR-sig+HEIDI-pass: $(awk -v t=$TOT -F'\t' 'NR>1 && $19!="NA" && $19+0<0.05/t && $20!="NA" && $20+0>0.05' smr_${DISEASE}_all.smr|wc -l)"
