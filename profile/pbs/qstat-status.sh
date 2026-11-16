#!/bin/bash
# Snakemake cluster-generic STATUS checker for Gadi/PBS.
# Snakemake calls:  qstat-status.sh <jobid>
# It must print exactly one of: running | success | failed
set -uo pipefail
JOBID=$1

# -x keeps finished jobs queryable for a while after they leave the queue
OUT=$(qstat -x -f "${JOBID}" 2>/dev/null) || { echo success; exit 0; }  # gone from history -> assume done

STATE=$(printf '%s\n' "$OUT" | awk -F'= ' '/job_state/ {print $2; exit}')
case "${STATE}" in
  R|Q|H|E|B|S|U|W|T|M) echo running ;;
  F|C)
      EC=$(printf '%s\n' "$OUT" | awk -F'= ' '/Exit_status/ {print $2; exit}')
      if [ "${EC:-1}" = "0" ]; then echo success; else echo failed; fi ;;
  *) echo running ;;   # unknown/transient -> keep polling rather than fail
esac
