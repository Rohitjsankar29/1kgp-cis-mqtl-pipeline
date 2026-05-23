#!/usr/bin/env python3
"""
10_annotate_functional.py

Annotate variant-CpG feature matrix with functional genomic evidence.

This script uses BED files downloaded from ENCODE, Roadmap Epigenomics,
or other public resources. It checks whether variants and CpGs overlap
regulatory features such as promoters, enhancers, ATAC-seq peaks, H3K27ac,
H3K4me3, CTCF, or chromatin states.

Inputs:
    --features        feature matrix TSV with variant_id and phenotype_id
    --cpg-bed         methylation phenotype BED.gz
    --annotation-dir  directory containing .bed/.bed.gz annotation tracks
    --out             annotated feature matrix TSV

Annotation file naming:
    The script uses the filename stem as the feature name.
    Example:
        annotations/ENCODE_ATAC_blood.bed.gz      -> overlap_ENCODE_ATAC_blood
        annotations/Roadmap_H3K27ac_blood.bed.gz  -> overlap_Roadmap_H3K27ac_blood

Output:
    Adds binary columns:
        variant_overlap_<track>
        cpg_overlap_<track>

Dependencies:
    Python only. For very large annotation sets, use bedtools in production.
"""

import argparse
import gzip
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)


def read_bed(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            chrom = p[0] if p[0].startswith("chr") else "chr" + p[0]
            try:
                start, end = int(p[1]), int(p[2])
            except ValueError:
                continue
            rows.append((chrom, start, end))
    return pd.DataFrame(rows, columns=["chr", "start", "end"])


def parse_variant_position(variant_id):
    s = str(variant_id)
    m = re.match(r"^(chr)?([0-9XYM]+)[:_](\d+)", s)
    if not m:
        return pd.Series({"variant_chr": np.nan, "variant_pos": np.nan})
    return pd.Series({"variant_chr": "chr" + m.group(2), "variant_pos": int(m.group(3))})


def load_cpg(cpg_bed):
    df = pd.read_csv(cpg_bed, sep="\t", compression="infer")
    chrom_col = "#chr" if "#chr" in df.columns else "chr"
    out = df[[chrom_col, "start", "end", "phenotype_id"]].copy()
    out = out.rename(columns={chrom_col: "cpg_chr", "start": "cpg_start", "end": "cpg_end"})
    return out


def sanitize_track_name(path):
    name = Path(path).name
    for suffix in [".bed.gz", ".bed", ".narrowPeak.gz", ".narrowPeak"]:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return re.sub(r"[^A-Za-z0-9_]+", "_", name)


def point_overlaps_bed(chroms, positions, bed):
    result = np.zeros(len(chroms), dtype=np.int8)
    bed_by_chr = {c: d[["start", "end"]].values for c, d in bed.groupby("chr")}
    for i, (chrom, pos) in enumerate(zip(chroms, positions)):
        if pd.isna(chrom) or pd.isna(pos) or chrom not in bed_by_chr:
            continue
        intervals = bed_by_chr[chrom]
        p = int(pos)
        # Good enough for pilot; for production use pyranges/bedtools
        if np.any((intervals[:, 0] <= p) & (intervals[:, 1] > p)):
            result[i] = 1
    return result


def interval_overlaps_bed(chroms, starts, ends, bed):
    result = np.zeros(len(chroms), dtype=np.int8)
    bed_by_chr = {c: d[["start", "end"]].values for c, d in bed.groupby("chr")}
    for i, (chrom, start, end) in enumerate(zip(chroms, starts, ends)):
        if pd.isna(chrom) or chrom not in bed_by_chr:
            continue
        intervals = bed_by_chr[chrom]
        s, e = int(start), int(end)
        if np.any((intervals[:, 0] < e) & (intervals[:, 1] > s)):
            result[i] = 1
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--features", required=True)
    p.add_argument("--cpg-bed", required=True)
    p.add_argument("--annotation-dir", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    features = pd.read_csv(args.features, sep="\t")
    if "variant_id" not in features.columns or "phenotype_id" not in features.columns:
        raise ValueError("Feature matrix must contain variant_id and phenotype_id")

    cpg = load_cpg(args.cpg_bed)
    features = features.merge(cpg, on="phenotype_id", how="left")

    if "variant_chr" not in features.columns or "variant_pos" not in features.columns:
        parsed = features["variant_id"].apply(parse_variant_position)
        features = pd.concat([features, parsed], axis=1)

    ann_paths = sorted(list(Path(args.annotation_dir).glob("*.bed")) + list(Path(args.annotation_dir).glob("*.bed.gz")) +
                       list(Path(args.annotation_dir).glob("*.narrowPeak")) + list(Path(args.annotation_dir).glob("*.narrowPeak.gz")))

    if not ann_paths:
        raise FileNotFoundError(f"No BED/narrowPeak files found in {args.annotation_dir}")

    for path in ann_paths:
        track = sanitize_track_name(path)
        log.info(f"Annotating {track}")
        bed = read_bed(path)
        if bed.empty:
            log.warning(f"Skipping empty annotation: {path}")
            continue

        features[f"variant_overlap_{track}"] = point_overlaps_bed(
            features["variant_chr"], features["variant_pos"], bed
        )
        features[f"cpg_overlap_{track}"] = interval_overlaps_bed(
            features["cpg_chr"], features["cpg_start"], features["cpg_end"], bed
        )

    features.to_csv(args.out, sep="\t", index=False)
    log.info(f"Written: {args.out} ({len(features):,} rows)")


if __name__ == "__main__":
    main()
