#!/usr/bin/env python3
"""compute_lambda_gw.py - genomic inflation (lambda) for the genome-wide 452
run: REAL permutation scan vs the PERMUTED-GENOTYPE NULL, per chromosome and
pooled. Writes calibration_gw.txt (the genome-wide analogue of the chr22/100
pilot's lambda 1.446 vs lambda_null 1.104 check).

Locate the real permutation outputs first, e.g.:
  ls /scratch/cy94/rs4477/1kgp-cis-mqtl/tensorqtl/genome_452/
then pass the matching glob as --real-glob.
"""
import argparse, glob, os
import numpy as np
import pandas as pd
from scipy.stats import chi2

ap = argparse.ArgumentParser()
ap.add_argument("--real-glob", required=True,
                help="glob for the real per-chrom permutation result files")
ap.add_argument("--null-glob",
                default="/g/data/cy94/rs4477/downstream/genome_452/chr*/null/*.null.seed*.txt.gz")
ap.add_argument("--out",
                default="/g/data/cy94/rs4477/downstream/genome_452/calibration_gw.txt")
a = ap.parse_args()

MED = chi2.ppf(0.5, 1)  # 0.4549

def lam(p):
    p = np.asarray(p, float)
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    return float(np.median(chi2.isf(p, 1)) / MED) if len(p) else float("nan")

def collect(pattern, label):
    rows, pooled = [], []
    for f in sorted(glob.glob(pattern)):
        d = pd.read_csv(f, sep="\t")
        pcol = "pval_beta" if "pval_beta" in d.columns else \
               ("pval_perm" if "pval_perm" in d.columns else "pval_nominal")
        p = d[pcol].to_numpy(float)
        rows.append((os.path.basename(f), len(p), lam(p)))
        pooled.append(p)
    pooled = np.concatenate(pooled) if pooled else np.array([])
    return rows, lam(pooled), len(pooled), label

out = []
for pattern, label in [(a.real_glob, "REAL"), (a.null_glob, "NULL (shuffled genotypes)")]:
    rows, l_all, n_all, lab = collect(pattern, label)
    out.append(f"\n{lab}: pooled lambda = {l_all:.3f}  ({n_all:,} CpGs, {len(rows)} files)")
    for name, n, l in rows:
        out.append(f"  {name:45s} n={n:>8,}  lambda={l:.3f}")

txt = ("GENOME-WIDE CALIBRATION (n=452) - real scan vs permuted-genotype null\n"
       + "=" * 68 + "\n" + "\n".join(out) + "\n\n"
       "Interpretation: if REAL lambda >> 1 while NULL collapses toward 1, the\n"
       "inflation depends on the genuine genotype-phenotype correspondence -\n"
       "real cis-mQTL signal, not confounding/software artefact (same logic as\n"
       "the chr22/100 pilot: 1.446 real vs 1.104 null).\n")
with open(a.out, "w") as o:
    o.write(txt)
print(txt)
