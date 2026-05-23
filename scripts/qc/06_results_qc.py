#!/usr/bin/env python3
"""
06_results_qc.py — QTL results QC: QQ plot, Manhattan, lambda GC, hit table

Reads:
  - permutation output: single .txt.gz OR a directory of chunked .txt.gz files
  - nominal parquets: directory containing chunked parquet files

Outputs (in --out-dir):
  - {label}.qqplot.pdf
  - {label}.manhattan.pdf
  - {label}.lambda_gc.txt
  - {label}.significant_mqtl.tsv
  - {label}.permutation_merged.txt.gz   (if input was chunked)

Changes from 100-sample version:
  - --permutation can now be a file OR directory of chunk*.txt.gz files (R46)
  - --nominal-dir points at chunked/ subdir; Manhattan samples rows to stay within
    memory rather than loading all 98 GB at once (R46)
  - Memory-efficient Manhattan: stream parquets, reservoir-sample 5M pairs max
"""

import argparse
import glob
import os
import sys
import logging
import numpy as np
import pandas as pd
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--permutation", required=True,
                   help="Permutation .txt.gz file OR directory of chunked .txt.gz files")
    p.add_argument("--nominal-dir", required=True, help="Dir with nominal parquet files")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--label", default="pilot.chr22")
    p.add_argument("--fdr-threshold", type=float, default=0.05)
    p.add_argument("--manhattan-max-pairs", type=int, default=5_000_000,
                   help="Max nominal pairs to plot (sampled per chunk)")
    p.add_argument("--manhattan-only", action="store_true",
                   help="Skip permutation merge/FDR — regenerate Manhattan only")
    return p.parse_args()


def compute_lambda_gc(pvals):
    """Compute genomic inflation factor λ from p-values."""
    pvals = np.array(pvals)
    pvals = pvals[~np.isnan(pvals) & (pvals > 0) & (pvals < 1)]
    if len(pvals) < 10:
        return np.nan
    from scipy.stats import chi2
    chi2_obs = chi2.ppf(1 - pvals, df=1)
    return np.median(chi2_obs) / chi2.ppf(0.5, df=1)


def storey_qvalue(pvals, pi0_lambda=0.5):
    """Simple Storey q-value via BH FDR (approximation when n is small)."""
    from statsmodels.stats.multitest import multipletests
    _, qvals, _, _ = multipletests(pvals, method='fdr_bh')
    return qvals


def load_permutation(perm_path, out_dir, label):
    """Load permutation results — single file or directory of chunks."""
    perm_path = Path(perm_path)
    if perm_path.is_dir():
        chunks = sorted(perm_path.glob("*.txt.gz"))
        if not chunks:
            raise FileNotFoundError(f"No .txt.gz files in {perm_path}")
        log.info(f"Merging {len(chunks)} permutation chunks from {perm_path}")
        dfs = []
        for c in chunks:
            df = pd.read_csv(c, sep='\t', index_col=0)
            dfs.append(df)
            log.info(f"  {c.name}: {len(df):,} rows")
        merged = pd.concat(dfs)
        # Deduplicate — chunks may overlap on phenotype boundaries
        merged = merged[~merged.index.duplicated(keep='first')]
        log.info(f"Merged permutation: {len(merged):,} CpGs (deduplicated)")
        merged_path = os.path.join(out_dir, f"{label}.permutation_merged.txt.gz")
        merged.to_csv(merged_path, sep='\t', compression='gzip')
        log.info(f"Merged permutation written: {merged_path}")
        return merged
    else:
        log.info(f"Loading permutation: {perm_path}")
        return pd.read_csv(perm_path, sep='\t', index_col=0)


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # --- Set up matplotlib once at start ---
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        _have_plt = True
    except ImportError:
        plt = None
        _have_plt = False
        log.warning("matplotlib not available — plots will be skipped")

    # --- Short-circuit for Manhattan-only re-runs ---
    if args.manhattan_only:
        log.info("--manhattan-only: skipping permutation merge and FDR")
        perm_df = None
    else:
        # --- Load permutation results (single file or chunked directory) ---
        perm_df = load_permutation(args.permutation, args.out_dir, args.label)
        log.info(f"Permutation: {len(perm_df)} CpGs")

    if perm_df is None:
        log.info("Skipping lambda GC, QQ plot and FDR (--manhattan-only)")
        pvals = None
    else:
        pvals = perm_df['pval_beta'].dropna()
        log.info(f"Valid p-values: {len(pvals)}")

    # --- Lambda GC ---
    if pvals is None:
        lam = None
        log.info("Skipping lambda GC")
    else:
        lam = compute_lambda_gc(pvals.values)
    if lam is not None:
        log.info(f"Lambda GC: {lam:.4f}")
        lambda_path = os.path.join(args.out_dir, f"{args.label}.lambda_gc.txt")
        with open(lambda_path, 'w') as f:
            f.write(f"label\tlambda_gc\tn_tests\n")
            f.write(f"{args.label}\t{lam:.4f}\t{len(pvals)}\n")
        log.info(f"Lambda GC written: {lambda_path}")

    # --- QQ plot ---
    if pvals is None:
        log.info("Skipping QQ plot (--manhattan-only)")
    elif not _have_plt:
        log.warning("Skipping QQ plot — matplotlib unavailable")
    else:
      try:
        sorted_obs = np.sort(-np.log10(pvals.values))[::-1]
        n = len(sorted_obs)
        expected = -np.log10(np.arange(1, n + 1) / (n + 1))

        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(expected, sorted_obs, s=5, alpha=0.6, color='steelblue', label='Observed')
        ax.plot([0, max(expected)], [0, max(expected)], 'r--', lw=1, label='Expected')
        ax.set_xlabel('Expected −log₁₀(p)', fontsize=12)
        ax.set_ylabel('Observed −log₁₀(p)', fontsize=12)
        ax.set_title(f'QQ plot — {args.label}\nλ = {lam:.4f} | n = {n:,}', fontsize=11)
        ax.legend(fontsize=10)
        plt.tight_layout()
        qq_path = os.path.join(args.out_dir, f"{args.label}.qqplot.pdf")
        plt.savefig(qq_path)
        plt.close()
        log.info(f"QQ plot: {qq_path}")
      except Exception as e:
        log.warning(f"QQ plot failed: {e}")

    # --- FDR + significant hits ---
    if pvals is None:
        log.info("Skipping FDR (--manhattan-only)")
    else:
      try:
        perm_df_valid = perm_df[perm_df['pval_beta'].notna()].copy()
        qvals = storey_qvalue(perm_df_valid['pval_beta'].values)
        perm_df_valid['qval'] = qvals
        sig = perm_df_valid[perm_df_valid['qval'] < args.fdr_threshold]
        log.info(f"Significant mQTLs (FDR < {args.fdr_threshold}): {len(sig)}")
        sig_path = os.path.join(args.out_dir, f"{args.label}.significant_mqtl.tsv")
        sig.to_csv(sig_path, sep='\t')
        log.info(f"Significant hits written: {sig_path}")
      except Exception as e:
        log.warning(f"FDR calculation failed: {e}")

    # --- Manhattan plot (memory-efficient: stream + sample nominal parquets) ---
    try:
        import pyarrow.parquet as pq
        parquet_files = sorted(glob.glob(os.path.join(args.nominal_dir, "*.parquet")))
        if parquet_files:
            log.info(f"Loading nominal parquets for Manhattan ({len(parquet_files)} files)")
            max_pairs = args.manhattan_max_pairs
            cols = ['variant_id', 'pval_nominal']

            # Row-group streaming: sample each parquet without loading it fully.
            # Peak memory = one row group at a time (~few hundred MB) (R48).
            per_chunk = max(1, max_pairs // len(parquet_files))
            parts = []
            total_pairs = 0
            rng = np.random.default_rng(42)

            for f in parquet_files:
                pf_meta = pq.ParquetFile(f)
                n_rows = pf_meta.metadata.num_rows
                total_pairs += n_rows

                # Pre-select which global row indices to keep (sorted for sequential scan)
                if n_rows <= per_chunk:
                    keep_idx = np.arange(n_rows)
                else:
                    keep_idx = np.sort(rng.choice(n_rows, per_chunk, replace=False))

                # Scan row groups, read only those containing selected rows
                sampled = []
                row_offset = 0
                for rg_i in range(pf_meta.metadata.num_row_groups):
                    rg_n = pf_meta.metadata.row_group(rg_i).num_rows
                    rg_end = row_offset + rg_n
                    local_idx = keep_idx[(keep_idx >= row_offset) & (keep_idx < rg_end)] - row_offset
                    if len(local_idx) > 0:
                        rg_df = pf_meta.read_row_group(rg_i, columns=cols).to_pandas()
                        sampled.append(rg_df.iloc[local_idx])
                    row_offset = rg_end

                if sampled:
                    parts.append(pd.concat(sampled, ignore_index=True))
                kept = len(parts[-1]) if parts else 0
                log.info(f"  {Path(f).name}: {n_rows:,} rows → kept {kept:,}")

            nom_df = pd.concat(parts, ignore_index=True)
            log.info(f"Total nominal pairs: {total_pairs:,} | Plotting: {len(nom_df):,}")

            def parse_pos(vid):
                parts = str(vid).split(':')
                try:
                    return int(parts[1]) if len(parts) >= 2 else np.nan
                except (ValueError, IndexError):
                    return np.nan

            nom_df['pos'] = nom_df['variant_id'].apply(parse_pos)
            nom_df = nom_df.dropna(subset=['pval_nominal', 'pos'])
            nom_df['neg_log10p'] = -np.log10(nom_df['pval_nominal'].clip(1e-300))

            fig, ax = plt.subplots(figsize=(10, 4))
            ax.scatter(nom_df['pos'] / 1e6, nom_df['neg_log10p'],
                       s=1, alpha=0.3, color='steelblue')
            ax.set_xlabel('chr22 position (Mb)', fontsize=11)
            ax.set_ylabel('−log₁₀(p nominal)', fontsize=11)
            ax.set_title(
                f'Manhattan — {args.label}\n{total_pairs:,} total pairs '
                f'(plotting {min(len(nom_df), max_pairs):,})', fontsize=11)
            plt.tight_layout()
            man_path = os.path.join(args.out_dir, f"{args.label}.manhattan.png")
            plt.savefig(man_path, dpi=150, bbox_inches='tight')
            plt.close()
            log.info(f"Manhattan plot: {man_path}")
        else:
            log.warning("No nominal parquet files found — skipping Manhattan")
    except Exception as e:
        log.warning(f"Manhattan plot failed: {e}")

    log.info("Results QC complete.")


if __name__ == '__main__':
    main()
