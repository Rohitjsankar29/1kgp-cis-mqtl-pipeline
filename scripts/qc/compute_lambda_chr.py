#!/usr/bin/env python3
"""compute_lambda_chr.py — genomic inflation (lambda_GC) for a single chromosome.

Reads:
  --real  — an EXISTING tensorQTL permutation output from the genome-wide 452
            run (no re-computation needed).
  --null  — the permuted-genotype-null permutation outputs produced fresh by
            05d_permuted_null.py (in the calibration PBS below), one file per
            seed.

Computes lambda_GC on both sets of per-CpG empirical p-values and writes a
2-line calibration report.

Interpretation:
    lambda_real  >> 1.0  AND  lambda_null  ~ 1.0   ==>  genuine cis-mQTL signal
    lambda_real  >> 1.0  AND  lambda_null  >> 1.0  ==>  residual confounding
    lambda_real  ~ 1.0                             ==>  no discovery signal (unusual)
"""
import argparse, glob, os, sys
import numpy as np
import pandas as pd
from scipy.stats import chi2

MED = chi2.ppf(0.5, 1)   # 0.4549 -- expected median chi2 under the null


def lam(p):
    p = np.asarray(p, float)
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    return float(np.median(chi2.isf(p, 1)) / MED) if len(p) else float("nan")


def load(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f"no files match: {pattern}")
    frames = []
    for f in files:
        # tensorQTL permutation output is tab-separated with a header
        try:
            df = pd.read_csv(f, sep="\t")
        except Exception as e:
            sys.exit(f"failed to read {f}: {e}")
        # tensorQTL uses pval_beta (beta-approx) as the per-CpG empirical p
        col = "pval_beta" if "pval_beta" in df.columns else (
              "pval_perm" if "pval_perm" in df.columns else None)
        if col is None:
            sys.exit(f"{f}: neither pval_beta nor pval_perm present\n"
                     f"columns: {list(df.columns)[:10]}")
        frames.append(df[[col]].rename(columns={col: "p"}))
    return pd.concat(frames, ignore_index=True)["p"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", required=True,
                    help="glob for the EXISTING real-pass permutation output "
                         "(e.g. /scratch/.../tensorqtl/genome_452/combined/chr22.permutation.txt)")
    ap.add_argument("--null", required=True,
                    help="glob for the permuted-null-pass outputs")
    ap.add_argument("--chr", required=True, help="chromosome label (for the report)")
    ap.add_argument("--out",  required=True)
    a = ap.parse_args()

    p_real = load(a.real)
    p_null = load(a.null)

    lam_real = lam(p_real)
    lam_null = lam(p_null)

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write(f"# calibration report -- chr{a.chr} at n=452\n")
        fh.write(f"n_cpgs_real\t{len(p_real)}\n")
        fh.write(f"n_cpgs_null\t{len(p_null)}\n")
        fh.write(f"lambda_real\t{lam_real:.4f}\n")
        fh.write(f"lambda_null\t{lam_null:.4f}\n")
        fh.write(f"# interpretation: lambda_real >> 1 with lambda_null ~ 1  =>  genuine cis-mQTL signal\n")
        fh.write(f"# chr22-100 pilot reference: lambda_real=1.446, lambda_null=1.104\n")
    print(open(a.out).read())


if __name__ == "__main__":
    main()
