#!/usr/bin/env python3
"""
03_aggregate_pilot.py — Aggregate per-sample bedMethyl → tensorQTL phenotype BED

Reads: per-sample *_chr22/combined.bedmethyl.gz from pilot samples
Outputs:
  - pilot.chr22.methylation.bed.gz  (sites × samples M-values)
  - pilot.chr22.meth_pca.tsv        (top N methylation PCs)
  - qc/coverage_per_sample.tsv
  - qc/site_missingness.tsv
  - qc/sample_correlation.tsv

Phenotype transformation: percent methylation → M-values (no quantile normalisation).
M = log2(β / (1 − β)), β = pct / 100, clamped to [0.001, 0.999].
"""

import argparse
import os
import sys
import gzip
import csv
import logging
import numpy as np
import pandas as pd
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(description="Aggregate bedMethyl → methylation matrix")
    p.add_argument("--samples-file", required=True)
    p.add_argument("--bedmethyl-dir", required=True)
    p.add_argument("--out-bed", required=True)
    p.add_argument("--out-pca", required=True)
    p.add_argument("--out-qc-dir", required=True)
    p.add_argument("--min-cov", type=int, default=5)
    p.add_argument("--min-sample-frac", type=float, default=0.8)
    p.add_argument("--n-pcs", type=int, default=3)
    p.add_argument("--region", default="chr22")
    p.add_argument("--threads", type=int, default=4)
    return p.parse_args()


def read_bedmethyl(path: Path, min_cov: int) -> pd.Series:
    """Read a bedMethyl.gz file → Series indexed by 'chr:pos', values = percent_modified."""
    data = {}
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 11:
                continue
            chrom, start, end = parts[0], parts[1], parts[2]
            coverage = int(parts[9])
            pct_mod = float(parts[10])
            if coverage < min_cov:
                continue
            key = f"{chrom}:{start}"
            data[key] = pct_mod
    return pd.Series(data, dtype=float)


def pct_to_mvalue(pct: np.ndarray, offset: float = 0.001) -> np.ndarray:
    """Convert percent methylation (0–100) to M-values: M = log2(beta / (1 - beta)).
    Beta is derived as pct / 100 and clamped to [offset, 1 - offset] to avoid log(0).
    offset = 0.001 matches the production script (01_aggregate_bedmethyl.py).
    """
    beta = np.clip(pct / 100.0, offset, 1.0 - offset)
    return np.log2(beta / (1.0 - beta))



def main():
    args = parse_args()

    # --- Load sample list ---
    samples_df = pd.read_csv(args.samples_file, sep='\t')
    sample_ids = samples_df['sample_id'].tolist()
    log.info(f"Samples: {sample_ids}")

    bedmethyl_dir = Path(args.bedmethyl_dir)
    qc_dir = Path(args.out_qc_dir)
    qc_dir.mkdir(parents=True, exist_ok=True)

    # --- Load per-sample bedMethyl ---
    sample_series = {}
    cov_stats = []
    for sid in sample_ids:
        bed_path = bedmethyl_dir / f"{sid}_chr22" / "combined.bedmethyl.gz"
        if not bed_path.exists():
            log.error(f"Missing bedMethyl: {bed_path}")
            sys.exit(1)
        log.info(f"Loading {sid}...")
        s = read_bedmethyl(bed_path, args.min_cov)
        sample_series[sid] = s
        cov_stats.append({
            'sample_id': sid,
            'n_cpg_covered': len(s),
            'mean_pct_meth': s.mean(),
            'median_pct_meth': s.median(),
        })
        log.info(f"  {sid}: {len(s):,} CpG sites (coverage >= {args.min_cov}x)")

    # --- Save per-sample coverage QC ---
    pd.DataFrame(cov_stats).to_csv(qc_dir / "coverage_per_sample.tsv", sep='\t', index=False)

    # --- Build sites × samples matrix ---
    log.info("Building matrix...")
    raw_df = pd.DataFrame(sample_series)   # rows=CpG sites, cols=samples
    log.info(f"Raw matrix: {raw_df.shape[0]:,} sites × {raw_df.shape[1]} samples")

    # --- Filter by completeness ---
    min_samples = int(np.ceil(args.min_sample_frac * len(sample_ids)))
    valid_count = raw_df.notna().sum(axis=1)
    mask = valid_count >= min_samples
    filtered_df = raw_df[mask].copy()
    log.info(f"After coverage filter (>={min_samples}/{len(sample_ids)} samples): {filtered_df.shape[0]:,} sites")

    # Save missingness QC
    miss_df = pd.DataFrame({
        'cpg_id': raw_df.index,
        'n_valid': valid_count.values,
        'frac_valid': (valid_count / len(sample_ids)).values,
        'passes_filter': mask.values,
    })
    miss_df.to_csv(qc_dir / "site_missingness.tsv", sep='\t', index=False)

    # --- Fill remaining NaN with sample median (impute) ---
    for col in filtered_df.columns:
        med = filtered_df[col].median()
        filtered_df[col] = filtered_df[col].fillna(med)

    # --- Convert to M-values ---
    log.info("Converting to M-values...")
    mval_df = filtered_df.apply(lambda col: pct_to_mvalue(col.values), axis=0)
    mval_df.index = filtered_df.index

    # --- Sample correlation QC ---
    corr = mval_df.corr(method='pearson')
    corr.to_csv(qc_dir / "sample_correlation.tsv", sep='\t')

    # --- Build tensorQTL BED format ---
    # Parse index "chr22:start" → chrom, start, end
    cpg_ids = mval_df.index.tolist()
    chroms, starts = zip(*[x.split(':') for x in cpg_ids])
    starts = [int(s) for s in starts]
    ends = [s + 1 for s in starts]

    bed_df = pd.DataFrame({
        '#chr': chroms,
        'start': starts,
        'end': ends,
        'phenotype_id': cpg_ids,
    })
    for sid in sample_ids:
        bed_df[sid] = mval_df[sid].values

    # Sort by position
    bed_df = bed_df.sort_values(['#chr', 'start'])
    log.info(f"Final methylation BED: {len(bed_df):,} CpG sites × {len(sample_ids)} samples")

    # --- Write BED (bgzip) ---
    out_bed = args.out_bed
    log.info(f"Writing {out_bed}...")
    with gzip.open(out_bed, 'wt') as fh:
        bed_df.to_csv(fh, sep='\t', index=False, float_format='%.6f')
    log.info("BED written (tabix indexing done by PBS job).")

    # --- PCA on methylation matrix ---
    log.info("Running methylation PCA...")
    from sklearn.decomposition import PCA
    # PCA input: samples × sites
    X = mval_df.T.values  # shape: n_samples × n_sites
    n_pcs = min(args.n_pcs, X.shape[0] - 1)
    pca = PCA(n_components=n_pcs)
    pcs = pca.fit_transform(X)  # shape: n_samples × n_pcs
    explained = pca.explained_variance_ratio_

    pca_df = pd.DataFrame(
        pcs.T,
        index=[f"MethPC{i+1}" for i in range(n_pcs)],
        columns=sample_ids,
    )
    pca_df.to_csv(args.out_pca, sep='\t')
    for i, ev in enumerate(explained):
        log.info(f"  MethPC{i+1} explained variance: {ev:.3f}")

    log.info("Done.")


if __name__ == '__main__':
    main()
