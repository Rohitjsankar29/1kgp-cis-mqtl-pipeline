
#!/usr/bin/env python3

"""
Build final variant-CpG feature matrix for prioritisation.

Combines:
  - TensorQTL association results
  - SuSiE fine-mapping results
  - Functional annotation features
  - Distance from variant to CpG
  - SV indicator if available
"""

import argparse
from pathlib import Path
import logging
import pandas as pd
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)


def parse_variant_pos(variant_id: str) -> int:
    try:
        return int(str(variant_id).split(":")[1])
    except Exception:
        return np.nan


def parse_cpg_pos(cpg_id: str) -> int:
    try:
        return int(str(cpg_id).split(":")[1])
    except Exception:
        return np.nan


def main():
    parser = argparse.ArgumentParser(description="Build variant-CpG prioritisation feature matrix")
    parser.add_argument("--tensorqtl", required=True, type=Path)
    parser.add_argument("--finemap", required=True, type=Path)
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--sv-list", type=Path, help="Optional list of SV variant IDs")

    args = parser.parse_args()

    log.info("Loading TensorQTL results")
    tqtl = pd.read_csv(args.tensorqtl, sep="\t")

    log.info("Loading fine-mapping results")
    fine = pd.read_csv(args.finemap, sep="\t")

    log.info("Loading annotations")
    annot = pd.read_csv(args.annotation, sep="\t")

    key_cols = ["phenotype_id", "variant_id"]

    df = tqtl.merge(fine, on=key_cols, how="left")
    df = df.merge(annot, on=key_cols, how="left")

    df["pip"] = df["pip"].fillna(0)
    df["in_credible_set"] = df["credible_set"].notna().astype(int)

    df["variant_pos"] = df["variant_id"].apply(parse_variant_pos)
    df["cpg_pos"] = df["phenotype_id"].apply(parse_cpg_pos)
    df["distance_to_cpg"] = (df["variant_pos"] - df["cpg_pos"]).abs()

    if args.sv_list and args.sv_list.exists():
        sv_ids = set(pd.read_csv(args.sv_list, header=None)[0].astype(str))
        df["is_sv"] = df["variant_id"].astype(str).isin(sv_ids).astype(int)
    else:
        df["is_sv"] = 0

    numeric_cols = [
        "pval_nominal",
        "slope",
        "pip",
        "distance_to_cpg",
        "promoter_overlap",
        "enhancer_overlap",
        "atac_overlap",
        "h3k27ac_overlap",
        "cpg_island_overlap",
        "is_sv",
        "in_credible_set",
    ]

    for col in numeric_cols:
        if col not in df.columns:
            df[col] = 0

    df["neg_log10_p"] = -np.log10(df["pval_nominal"].clip(lower=1e-300))

    keep_cols = key_cols + [
        "slope",
        "pval_nominal",
        "neg_log10_p",
        "pip",
        "credible_set",
        "in_credible_set",
        "distance_to_cpg",
        "is_sv",
        "promoter_overlap",
        "enhancer_overlap",
        "atac_overlap",
        "h3k27ac_overlap",
        "cpg_island_overlap",
    ]

    df = df[keep_cols].copy()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    log.info(f"Feature matrix written: {args.out}")
    log.info(f"Rows: {len(df):,}")


if __name__ == "__main__":
    main()
  
