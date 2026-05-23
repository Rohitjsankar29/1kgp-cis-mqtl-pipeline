#!/usr/bin/env python3
"""
11_score_prioritise_variants.py

Compute a simple CpG-centric prioritisation score from the feature matrix.

This implements the dissertation's starting weighted scoring model:
    score_i = sum_j w_j * standardised(feature_ij)

It combines:
    - association strength
    - fine-mapping PIP
    - distance to CpG
    - SV evidence
    - functional annotation overlaps

Input:
    --features annotated feature matrix TSV

Output:
    ranked_prioritised_variants.tsv
"""

import argparse
import logging
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)


def zscore(s):
    s = pd.to_numeric(s, errors="coerce")
    if s.notna().sum() < 2 or s.std(skipna=True) == 0:
        return pd.Series(np.zeros(len(s)), index=s.index)
    return (s - s.mean(skipna=True)) / s.std(skipna=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--features", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--top-per-cpg", type=int, default=20)
    args = p.parse_args()

    df = pd.read_csv(args.features, sep="\t")

    score = pd.Series(0.0, index=df.index)

    # Association strength: smaller p-value is better.
    for pcol in ["pval_nominal", "pval_beta", "p_value", "pval"]:
        if pcol in df.columns:
            df["neglog10p"] = -np.log10(pd.to_numeric(df[pcol], errors="coerce").clip(lower=1e-300))
            score += 2.0 * zscore(df["neglog10p"])
            break

    # Effect size magnitude
    for bcol in ["slope", "beta", "effect_size", "b"]:
        if bcol in df.columns:
            score += 1.0 * zscore(pd.to_numeric(df[bcol], errors="coerce").abs())
            break

    # Fine-mapping evidence
    if "pip" in df.columns:
        score += 3.0 * zscore(df["pip"])

    # Distance: closer is better, so use negative distance
    for dcol in ["distance_to_cpg", "sv_distance_to_cpg", "distance"]:
        if dcol in df.columns:
            score += 1.0 * zscore(-pd.to_numeric(df[dcol], errors="coerce"))
            break

    # SV evidence
    if "is_sv" in df.columns:
        score += 0.75 * pd.to_numeric(df["is_sv"], errors="coerce").fillna(0)

    if "sv_overlaps_cpg" in df.columns:
        score += 1.0 * pd.to_numeric(df["sv_overlaps_cpg"], errors="coerce").fillna(0)

    # Functional annotations: any overlap feature contributes
    overlap_cols = [c for c in df.columns if c.startswith("variant_overlap_") or c.startswith("cpg_overlap_")]
    if overlap_cols:
        df["n_functional_overlaps"] = df[overlap_cols].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1)
        score += 1.5 * zscore(df["n_functional_overlaps"])

    df["prioritisation_score"] = score

    sort_cols = ["phenotype_id", "prioritisation_score"]
    if "phenotype_id" in df.columns:
        df = df.sort_values(sort_cols, ascending=[True, False])
        df["rank_within_cpg"] = df.groupby("phenotype_id")["prioritisation_score"].rank(ascending=False, method="first")
        df = df[df["rank_within_cpg"] <= args.top_per_cpg]
    else:
        df = df.sort_values("prioritisation_score", ascending=False)

    df.to_csv(args.out, sep="\t", index=False)
    log.info(f"Written: {args.out} ({len(df):,} rows)")


if __name__ == "__main__":
    main()
