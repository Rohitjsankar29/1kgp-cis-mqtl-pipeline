#!/bin/bash
# Snakemake cluster-generic SUBMIT wrapper for Gadi/PBS.
# Snakemake calls:  qsub.sh <queue> <ncpus> <mem_mb> <walltime> <gpu> <rule> <jobscript>
# It must print ONLY the PBS job id on stdout (qsub already does this).
# NOTE: mem is passed as an integer (MB) to avoid the space in "16 GB" splitting args.
set -euo pipefail
QUEUE=$1; NCPUS=$2; MEM_MB=$3; WALL=$4; GPU=$5; RULE=$6; JOBSCRIPT=$7

L="ncpus=${NCPUS},mem=${MEM_MB}MB,walltime=${WALL},storage=gdata/cy94+scratch/cy94"
# only request GPUs when the rule actually asks for them (normal queue rejects ngpus>0)
if [ "${GPU}" -gt 0 ]; then L="ngpus=${GPU},${L}"; fi

# Forward the env vars the in-job snakemake needs (PBS jobs start with a clean
# environment and do NOT inherit the login shell's exports). Without this, the
# in-job snakemake writes its source cache to ~/.cache and dies on home quota.
exec qsub \
  -P cy94 \
  -q "${QUEUE}" \
  -l "${L}" \
  -l wd \
  -j oe \
  -o /scratch/cy94/rs4477/logs/ \
  -N "smk_${RULE}" \
  -v XDG_CACHE_HOME=/g/data/cy94/rs4477/.cache,MAMBA_ROOT_PREFIX=/g/data/cy94/rs4477/micromamba \
  "${JOBSCRIPT}"
