</> Python

#!/usr/bin/env python3

"""
Annotate variant-CpG pairs with regulatory genomic features.

Inputs:
  - Variant-CpG table from TensorQTL / fine-mapping
  - Variant BED file
  - Annotation BED files such as promoters, enhancers, ATAC peaks, histone marks

Output:
  - Annotated variant-CpG feature table
"""

import argparse
import logging
from pathlib import Path
import pandas as pd
import pyranges as pr

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)


def read_variant_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"variant_id", "phenotype_id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    return df


def variant_id_to_bed(df: pd.DataFrame) -> pd.DataFrame:
    """
    Expected variant_id format: chr:pos:ref:alt
    Converts variant position to BED-like interval.
    """
    parts = df["variant_id"].str.split(":", expand=True)
    df["Chromosome"] = parts[0]
    df["Start"] = parts[1].astype(int) - 1
    df["End"] = parts[1].astype(int)
    return df


def add_overlap_feature(df: pd.DataFrame, annotation_bed: Path, feature_name: str) -> pd.DataFrame:
    log.info(f"Adding annotation: {feature_name}")

    variants = pr.PyRanges(df[["Chromosome", "Start", "End", "variant_id"]])
    annot = pr.read_bed(str(annotation_bed))

    overlap = variants.join(annot)
    overlap_ids = set(overlap.df["variant_id"])

    df[feature_name] = df["variant_id"].isin(overlap_ids).astype(int)
    return df


def main():
    parser = argparse.ArgumentParser(description="Annotate variant-CpG pairs")
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)

    parser.add_argument("--promoters", type=Path)
    parser.add_argument("--enhancers", type=Path)
    parser.add_argument("--atac", type=Path)
    parser.add_argument("--h3k27ac", type=Path)
    parser.add_argument("--cpg-islands", type=Path)

    args = parser.parse_args()

    df = read_variant_table(args.input)
    df = variant_id_to_bed(df)

    annotations = {
        "promoter_overlap": args.promoters,
        "enhancer_overlap": args.enhancers,
        "atac_overlap": args.atac,
        "h3k27ac_overlap": args.h3k27ac,
        "cpg_island_overlap": args.cpg_islands,
    }

    for feature, bed in annotations.items():
        if bed and bed.exists():
            df = add_overlap_feature(df, bed, feature)
        else:
            log.warning(f"Skipping {feature}; BED file not provided")
            df[feature] = 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    log.info(f"Written: {args.out}")


if __name__ == "__main__":
    main()
