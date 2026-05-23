#!/usr/bin/env python3
"""
01_aggregate_bedmethyl.py
=========================
Aggregate per-sample modkit bedMethyl files into a sites × samples matrix
suitable for QTL analysis with tensorQTL.

Inputs:
    - Directory of <sample>_combined.bedmethyl[.gz] files
    - OR a manifest TSV: sample_id, bedmethyl_path

Outputs:
    - methylation_matrix.bed.gz   — sites × samples, beta values
    - methylation_matrix_Mval.bed.gz — sites × samples, M-values (for QTL)
    - methylation_matrix_hp1.bed.gz  — HP1 haplotype (optional)
    - methylation_matrix_hp2.bed.gz  — HP2 haplotype (optional)
    - site_qc.tsv                 — per-site coverage and missingness stats

M-value conversion: M = log2(beta / (1 - beta))
  - Used for statistical testing (more Gaussian than beta)
  - Beta values retained for visualisation/reporting

bedMethyl format (modkit v0.6.1 output):
  col 1:  chrom
  col 2:  start (0-based)
  col 3:  end
  col 4:  mod_code (m = 5mC)
  col 5:  coverage (valid calls)
  col 6:  strand (. = combined)
  col 7:  start (display)
  col 8:  end (display)
  col 9:  color
  col 10: coverage (all calls, incl filtered)
  col 11: percent_modified (beta * 100)
  col 12: n_modified
  col 13: n_canonical
  col 14: n_other_modified
  col 15: n_delete
  col 16: n_fail
  col 17: n_diff
  col 18: n_no_call

Usage:
    python 01_aggregate_bedmethyl.py \\
        --input-dir /path/to/bedmethyl/ \\
        --pattern "*_combined.bedmethyl.gz" \\
        --output-dir /path/to/output/ \\
        --min-coverage 5 \\
        --min-samples 0.8 \\
        --haplotype  # also build HP1/HP2 matrices

    python 01_aggregate_bedmethyl.py \\
        --manifest manifest.tsv \\
        --output-dir /path/to/output/

Author: Kim Navarro (k1mnav)
Date: 2026-02-26
"""

import argparse
import gzip
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

BEDMETHYL_COLS = [
    "chrom", "start", "end", "mod_code", "score", "strand",
    "display_start", "display_end", "color",
    "coverage", "pct_modified",
    "n_modified", "n_canonical", "n_other_modified",
    "n_delete", "n_fail", "n_diff", "n_no_call",
]


def read_bedmethyl(path: Path, min_coverage: int = 5) -> pd.DataFrame:
    """Read a bedMethyl file, apply coverage filter, return site:beta series."""
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=BEDMETHYL_COLS,
            dtype={
                "chrom": str, "start": int, "end": int,
                "coverage": int, "pct_modified": float,
            },
        )
    except Exception as e:
        log.error(f"Failed to read {path}: {e}")
        raise

    # Coverage filter
    df = df[df["coverage"] >= min_coverage].copy()

    # Site key: chr:start (0-based, modkit outputs 0-based)
    df["site"] = df["chrom"] + ":" + df["start"].astype(str)

    # Beta value (0–1)
    df["beta"] = df["pct_modified"] / 100.0

    return df.set_index("site")[["chrom", "start", "end", "beta", "coverage"]]


def beta_to_mvalue(beta: pd.Series, offset: float = 0.001) -> pd.Series:
    """
    Convert beta values to M-values: M = log2(beta / (1 - beta))
    Clamp beta to [offset, 1-offset] to avoid inf values.
    """
    beta_clamped = beta.clip(offset, 1 - offset)
    return np.log2(beta_clamped / (1 - beta_clamped))


def build_matrix(
    sample_files: dict,  # {sample_id: Path}
    min_coverage: int = 5,
    min_sample_frac: float = 0.8,
    value: str = "beta",  # "beta" or "mvalue"
) -> pd.DataFrame:
    """
    Build a sites × samples matrix.

    Sites present in < min_sample_frac of samples are dropped.
    Remaining missing values are imputed with the site mean.
    """
    n = len(sample_files)
    log.info(f"Building matrix from {n} samples (min_cov={min_coverage}, min_samples={min_sample_frac:.0%})")

    beta_dict = {}
    for sid, path in sample_files.items():
        log.info(f"  Reading {sid} ← {path.name}")
        df = read_bedmethyl(path, min_coverage=min_coverage)
        beta_dict[sid] = df["beta"]

    # Align on common sites
    matrix = pd.DataFrame(beta_dict)  # sites × samples (beta)
    log.info(f"  Sites before filtering: {len(matrix):,} (union across samples)")

    # Drop sites missing in too many samples
    min_obs = int(np.ceil(min_sample_frac * n))
    matrix = matrix.dropna(thresh=min_obs)
    log.info(f"  Sites after min-sample filter ({min_obs}/{n}): {len(matrix):,}")

    # Impute missing with site mean (mean across samples that have coverage)
    row_means = matrix.mean(axis=1)
    matrix = matrix.T.fillna(row_means).T
    log.info(f"  Imputed missing values with site mean")

    if value == "mvalue":
        log.info("  Converting to M-values")
        matrix = matrix.apply(beta_to_mvalue)

    return matrix


def matrix_to_phenotype_bed(matrix: pd.DataFrame, output_path: Path) -> None:
    """
    Write matrix as tensorQTL-compatible phenotype BED file.

    Format: #chr  start  end  phenotype_id  sample1  sample2 ...
    - 0-based, half-open coordinates
    - Sorted by chrom, start
    - Gzipped
    """
    # Parse site coordinates from index (chr:start)
    coords = matrix.index.str.split(":", expand=True).to_frame(index=False)
    coords.columns = ["chrom", "pos"]
    coords["pos"] = coords["pos"].astype(int)
    out = pd.DataFrame({
        "#chr": coords["chrom"].values,
        "start": coords["pos"].values,
        "end": coords["pos"].values + 1,
        "phenotype_id": matrix.index,
    })
    out = pd.concat([out, matrix.reset_index(drop=True)], axis=1)

    # Sort
    chrom_order = {f"chr{i}": i for i in list(range(1, 23)) + ["X", "Y", "M"]}
    out["_sort"] = out["#chr"].map(chrom_order).fillna(99).infer_objects(copy=False)
    out = out.sort_values(["_sort", "start"]).drop(columns=["_sort"])

    out.to_csv(output_path, sep="\t", index=False, float_format="%.6f",
               compression="gzip" if str(output_path).endswith(".gz") else None)
    log.info(f"  Wrote {len(out):,} sites × {len(matrix.columns)} samples → {output_path}")


def collect_sample_files(input_dir: Path, pattern: str) -> dict:
    """Glob files matching pattern, extract sample IDs from filename."""
    files = sorted(input_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matching {pattern} in {input_dir}")
    result = {}
    for f in files:
        # Extract sample ID: everything before first underscore in stem
        stem = f.name.replace(".bedmethyl.gz", "").replace(".bedmethyl", "")
        # Remove suffix like _combined, _hp1, _hp2
        for suffix in ("_combined", "_hp1", "_hp2"):
            stem = stem.replace(suffix, "")
        result[stem] = f
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate modkit bedMethyl files into a sites × samples matrix"
    )
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--input-dir", type=Path,
                     help="Directory containing bedMethyl files")
    src.add_argument("--manifest", type=Path,
                     help="TSV with columns: sample_id, bedmethyl_path")

    parser.add_argument("--pattern", default="*_combined.bedmethyl.gz",
                        help="Glob pattern for combined bedMethyl files (default: *_combined.bedmethyl.gz)")
    parser.add_argument("--hp1-pattern", default="*_hp1.bedmethyl.gz")
    parser.add_argument("--hp2-pattern", default="*_hp2.bedmethyl.gz")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-coverage", type=int, default=5,
                        help="Minimum read coverage per site per sample (default: 5)")
    parser.add_argument("--min-samples", type=float, default=0.8,
                        help="Minimum fraction of samples with coverage for a site to be retained (default: 0.8)")
    parser.add_argument("--haplotype", action="store_true",
                        help="Also build HP1 and HP2 haplotype matrices")
    parser.add_argument("--no-mvalue", action="store_true",
                        help="Skip M-value matrix output")

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Collect sample files
    if args.input_dir:
        combined_files = collect_sample_files(args.input_dir, args.pattern)
        if args.haplotype:
            hp1_files = collect_sample_files(args.input_dir, args.hp1_pattern)
            hp2_files = collect_sample_files(args.input_dir, args.hp2_pattern)
    else:
        manifest = pd.read_csv(args.manifest, sep="\t")
        combined_files = {row.sample_id: Path(row.bedmethyl_path)
                         for _, row in manifest.iterrows()}

    log.info(f"Found {len(combined_files)} samples")

    # --- Combined beta matrix ---
    beta_matrix = build_matrix(
        combined_files,
        min_coverage=args.min_coverage,
        min_sample_frac=args.min_samples,
        value="beta",
    )
    matrix_to_phenotype_bed(
        beta_matrix,
        args.output_dir / "methylation_beta.bed.gz",
    )

    # --- M-value matrix (for QTL) ---
    if not args.no_mvalue:
        mval_matrix = beta_matrix.apply(beta_to_mvalue)
        matrix_to_phenotype_bed(
            mval_matrix,
            args.output_dir / "methylation_Mval.bed.gz",
        )

    # --- Site QC ---
    site_qc = pd.DataFrame({
        "n_samples": beta_matrix.notna().sum(axis=1),
        "mean_beta": beta_matrix.mean(axis=1),
        "sd_beta": beta_matrix.std(axis=1),
        "min_beta": beta_matrix.min(axis=1),
        "max_beta": beta_matrix.max(axis=1),
    })
    site_qc.to_csv(args.output_dir / "site_qc.tsv", sep="\t", float_format="%.4f")
    log.info(f"Site QC written → {args.output_dir}/site_qc.tsv")

    # --- Haplotype matrices ---
    if args.haplotype and args.input_dir:
        for tag, files in [("hp1", hp1_files), ("hp2", hp2_files)]:
            log.info(f"Building {tag.upper()} matrix...")
            hp_matrix = build_matrix(
                files,
                min_coverage=args.min_coverage,
                min_sample_frac=args.min_samples,
                value="mvalue",
            )
            matrix_to_phenotype_bed(
                hp_matrix,
                args.output_dir / f"methylation_Mval_{tag}.bed.gz",
            )

    log.info("Done.")


if __name__ == "__main__":
    main()
