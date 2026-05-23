#!/usr/bin/env python3
"""
05_run_tensorqtl_pilot.py — tensorQTL cis-mQTL: nominal + permutation pass

Decision R7 conventions:
  - plink1 bed/bim/fam format (not pgen)
  - covariates_df: file is rows=covariates, cols=samples
    → transpose to rows=samples, cols=covariates before passing
  - cis.map_nominal: prefix is positional arg #5; no threshold param
  - cis.map_cis: nperm (not n_permutations)
"""

import argparse
import gzip
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
    p.add_argument("--plink-prefix", required=True, help="PLINK1 bed/bim/fam prefix")
    p.add_argument("--phenotype-bed", required=True, help="Methylation BED (bgzip, tabix)")
    p.add_argument("--covariates", required=True, help="Covariate TSV (rows=cov, cols=samples)")
    p.add_argument("--nominal-out", required=True, help="Directory for nominal output parquet")
    p.add_argument("--permutation-out", required=True, help="Output path for permutation .txt.gz")
    p.add_argument("--cis-window", type=int, default=1_000_000)
    p.add_argument("--n-permutations", type=int, default=1000)
    p.add_argument("--region", default="chr22")
    p.add_argument("--skip-nominal", action="store_true",
                   help="Skip nominal pass (too large for full-chr pilot)")
    # Chunking args for array jobs (R26 fix: were referenced but not registered)
    p.add_argument("--n-chunks", type=int, default=1,
                   help="Total number of phenotype chunks (for array jobs)")
    p.add_argument("--chunk-idx", type=int, default=0,
                   help="0-based chunk index (PBS_ARRAY_INDEX)")
    # Trans mQTL args (added for 100-sample test — chr22 intra-chromosomal trans)
    p.add_argument("--run-trans", action="store_true",
                   help="Run trans mQTL pass after cis (pval threshold sparse output)")
    p.add_argument("--trans-pval-threshold", type=float, default=1e-5,
                   help="Sparse trans output: store only pairs with p < threshold")
    p.add_argument("--trans-batch-size", type=int, default=20000,
                   help="Genotype batch size for map_trans GPU computation. "
                        "Default 20000 safe for chunked trans (58k phenotypes: 4.6 GB/batch, "
                        "fits V100 32 GB). For monolithic trans (584k phenotypes): use 8000. "
                        "Formula: safe_batch = floor(28e9 / (n_phenotypes * 4)). See R39, R42.")
    p.add_argument("--trans-output", default=None,
                   help="Output parquet path for trans results")
    # NEW: chunked trans — for 08c_tensorqtl_trans_array.pbs (R37 fix)
    # When --run-trans AND --trans-chunk-idx are both set, subsets phenotypes to
    # the given chunk before calling map_trans. Same chunk logic as cis chunking.
    # This solves the CPU-bound walltime problem: each chunk is ~58k phenotypes (~2h).
    p.add_argument("--trans-chunk-idx", type=int, default=None,
                   help="If set with --run-trans: subset phenotypes to this chunk index "
                        "before calling map_trans (same logic as --chunk-idx for cis). "
                        "Use with --n-chunks to divide the phenotype file across array jobs.")
    # NEW: nominal-mode — for 09_tensorqtl_nominal_streaming.pbs (R36 fix)
    # When set, runs map_nominal on the phenotype chunk and writes output to
    # --nominal-out (parquet). Skips permutation pass entirely. Designed to be
    # run as a 10-chunk array job so each chunk fits within the per-node RAM budget
    # (~58k CpGs × ~12k cis variants per chunk << 7B pairs for full phenotype set).
    p.add_argument("--nominal-mode", action="store_true",
                   help="Run nominal pass only (map_nominal on phenotype chunk). "
                        "Skips permutation pass. Writes parquet to --nominal-out. "
                        "Use with --chunk-idx / --n-chunks for chunked array submission. "
                        "Solves R36 OOM: splits the full 7B-pair nominal pass into "
                        "manageable per-chunk sizes.")
    return p.parse_args()


def load_phenotype_bed(path: str, args):
    """Load tensorQTL phenotype BED for v1.0.10 (R26).
    pos_df uses chr+pos columns required by get_cis_ranges.
    """
    df = pd.read_csv(path, sep="\t",
                     dtype={"#chr": str, "start": int, "end": int})
    df = df.sort_values(["#chr", "start"]).reset_index(drop=True)
    df = df.set_index(df.columns[3])  # phenotype_id as index
    phenotype_pos = df[["#chr", "end"]].copy()
    phenotype_pos.columns = ["chr", "pos"]  # pos not tss (R26)
    phenotype_df = df.drop(columns=["#chr", "start", "end"])
    phenotype_df = phenotype_df.astype(float)
    log.info(f"Phenotypes: {phenotype_df.shape[0]} CpGs x {phenotype_df.shape[1]} samples")
    # Phenotype chunking for array jobs
    if args.n_chunks > 1:
        total = phenotype_df.shape[0]
        chunk_size = (total + args.n_chunks - 1) // args.n_chunks
        start = args.chunk_idx * chunk_size
        end = min(start + chunk_size, total)
        phenotype_df = phenotype_df.iloc[start:end]
        phenotype_pos = phenotype_pos.iloc[start:end]
        log.info(f"Chunk {args.chunk_idx}/{args.n_chunks}: {len(phenotype_df)} phenotypes ({start}-{end})")
    return phenotype_df, phenotype_pos

def main():
    args = parse_args()

    import torch
    log.info(f"PyTorch: {torch.__version__} | CUDA: {torch.cuda.is_available()}")
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    log.info(f"Using device: {device}")

    import tensorqtl
    from tensorqtl import cis, genotypeio
    log.info(f"tensorQTL: {tensorqtl.__version__}")

    # --- Load genotypes (PLINK1 format, Decision R6) ---
    log.info(f"Loading genotypes: {args.plink_prefix}")
    pr = genotypeio.PlinkReader(args.plink_prefix)
    variant_df = pr.bim.copy()
    variant_df = variant_df.rename(columns={'chrom': 'chr', 'pos': 'pos'})
    log.info(f"Variants: {len(variant_df):,}")

    # --- Load phenotypes ---
    log.info(f"Loading phenotypes: {args.phenotype_bed}")
    phenotype_df, phenotype_pos = load_phenotype_bed(args.phenotype_bed, args)

    # --- Load covariates (Decision R7: file=rows×cov, cols=samples → transpose) ---
    log.info(f"Loading covariates: {args.covariates}")
    cov_raw = pd.read_csv(args.covariates, sep='\t', index_col=0)
    log.info(f"Covariates file shape: {cov_raw.shape} (rows=covariates, cols=samples)")

    # Align samples: intersect genotype, phenotype, covariate samples
    geno_samples = pr.fam['iid'].tolist()
    pheno_samples = phenotype_df.columns.tolist()
    cov_samples = cov_raw.columns.tolist()
    common_samples = sorted(set(geno_samples) & set(pheno_samples) & set(cov_samples))
    log.info(f"Common samples: {len(common_samples)} — {common_samples}")

    if len(common_samples) < 5:
        log.error(f"Too few common samples: {len(common_samples)}")
        sys.exit(1)

    # Subset + align
    phenotype_df = phenotype_df[common_samples]
    cov_df = cov_raw[common_samples].T  # Decision R7: transpose → rows=samples, cols=covariates
    log.info(f"Covariates for tensorQTL: {cov_df.shape} (rows=samples, cols=covariates)")

    # Load genotype matrix (tensorQTL 1.0.10: load_genotypes() returns single DataFrame)
    # Use bim["snp"] for variant IDs — plink2 assigns chr:pos:ref:alt via --set-all-var-ids.
    # (Previously reconstructed as chrom:pos:a0:a1 when snp column was all '.'
    #  at n=10 pilot; now that IDs are assigned, use them directly for correct REF/ALT order.)
    bim = pr.bim.copy()
    bim["vid"] = bim["snp"]
    genotype_df = pr.load_genotypes()
    genotype_df.index = bim["vid"].values
    variant_df = bim.set_index("vid")[["chrom", "pos"]]
    # Subset to common samples
    genotype_df = genotype_df[common_samples]
    log.info(f"Genotype matrix loaded: {genotype_df.shape}")

    # =========================================================================
    # NEW: nominal-mode — chunked nominal pass (solves R36 OOM)
    # When --nominal-mode is set, we run map_nominal on the current phenotype
    # chunk (already sliced by load_phenotype_bed above) and write parquet.
    # Permutation pass is skipped entirely — nominal-mode is the full output.
    # Designed to be submitted as a 10-job array (08c approach) where each
    # chunk covers ~58k CpGs, keeping per-node memory well below the node limit.
    # =========================================================================
    if args.nominal_mode:
        log.info("=== Nominal-mode: chunked nominal pass (R36 fix) ===")
        log.info(f"Chunk {args.chunk_idx}/{args.n_chunks}: {len(phenotype_df)} phenotypes")

        if args.nominal_out is None:
            log.error("--nominal-out is required with --nominal-mode")
            sys.exit(1)

        # nominal_out is treated as a file path (parquet) in nominal-mode,
        # not a directory prefix as in the legacy inline nominal pass.
        nom_out_path = args.nominal_out
        os.makedirs(os.path.dirname(os.path.abspath(nom_out_path)), exist_ok=True)

        # Use a temporary directory prefix for map_nominal (writes parquet files
        # named <prefix>.cis_qtl_pairs.<chrom>.parquet internally), then move
        # the result to the target path.
        import tempfile, glob, shutil
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = os.path.join(tmpdir, f"nominal_chunk{args.chunk_idx}")
            cis.map_nominal(
                genotype_df, variant_df, phenotype_df, phenotype_pos, prefix,
                covariates_df=cov_df, window=args.cis_window,
            )
            # Collect all parquet files written by map_nominal and concatenate
            parquet_files = sorted(glob.glob(glob.escape(prefix) + "*.parquet"))
            if not parquet_files:
                log.error(f"map_nominal produced no parquet output in {tmpdir}")
                sys.exit(1)
            log.info(f"map_nominal wrote {len(parquet_files)} parquet shard(s)")
            if len(parquet_files) == 1:
                shutil.move(parquet_files[0], nom_out_path)
            else:
                # Multiple shards (one per chromosome) — concat and write single parquet
                dfs = [pd.read_parquet(f) for f in parquet_files]
                pd.concat(dfs, ignore_index=True).to_parquet(nom_out_path, index=False)

        n_pairs = len(pd.read_parquet(nom_out_path))
        log.info(f"Nominal chunk {args.chunk_idx}: {n_pairs:,} pairs → {nom_out_path}")
        log.info("TENSORQTL_OK")
        return  # nominal-mode exits here — no permutation pass

    # --- Nominal pass (inline, legacy — skipped for production by --skip-nominal) ---
    if not args.skip_nominal:
        log.info("=== Nominal pass ===")
        os.makedirs(args.nominal_out, exist_ok=True)
        prefix = os.path.join(args.nominal_out, "pilot.chr22")
        cis.map_nominal(
            genotype_df, variant_df, phenotype_df, phenotype_pos, prefix,
            covariates_df=cov_df, window=args.cis_window,
        )
        log.info("Nominal pass complete.")
    else:
        log.info("Skipping nominal pass (--skip-nominal).")

    # --- Permutation pass ---
    log.info("=== Permutation pass ===")
    perm_df = cis.map_cis(
        genotype_df,
        variant_df,
        phenotype_df,
        phenotype_pos,
        covariates_df=cov_df,
        nperm=args.n_permutations,
        window=args.cis_window,
    )
    log.info(f"Permutation pass complete: {len(perm_df)} CpGs tested")

    # Write permutation results
    perm_out = args.permutation_out
    perm_df.to_csv(perm_out, sep='\t', index=True,
                   compression='gzip', float_format='%.8g')
    log.info(f"Permutation results: {perm_out}")

    # Summary stats
    valid_pvals = perm_df['pval_beta'].dropna()
    log.info(f"  CpGs with valid p-value: {len(valid_pvals)}")
    if len(valid_pvals) > 0:
        log.info(f"  Min pval_beta: {valid_pvals.min():.4g}")
        log.info(f"  Median pval_beta: {valid_pvals.median():.4g}")

    # =========================================================================
    # Trans pass (optional)
    # Two modes:
    #   (A) Legacy: --run-trans without --trans-chunk-idx
    #       Runs on the full (already cis-chunked) phenotype_df.
    #       Used in 08b_tensorqtl_trans.pbs.
    #   (B) NEW chunked trans: --run-trans --trans-chunk-idx INT (R37 fix)
    #       Subsets phenotype_df to a secondary chunk before calling map_trans.
    #       This allows the full 584k-phenotype trans pass to be split across
    #       10 array jobs (~58k phenotypes each), each completing in ~2h.
    #       Used in 08c_tensorqtl_trans_array.pbs.
    # =========================================================================
    if args.run_trans:
        if args.trans_output is None:
            log.error("--trans-output required with --run-trans")
            sys.exit(1)

        # NEW: chunked trans — apply phenotype chunk subsetting for 08c
        if args.trans_chunk_idx is not None:
            log.info(f"=== Trans pass — chunk {args.trans_chunk_idx}/{args.n_chunks} "
                     f"(R37 chunked array mode) ===")
            total = phenotype_df.shape[0]
            chunk_size = (total + args.n_chunks - 1) // args.n_chunks
            t_start = args.trans_chunk_idx * chunk_size
            t_end = min(t_start + chunk_size, total)
            # Re-slice phenotype_df to the trans chunk.
            # Note: if --n-chunks was already applied via load_phenotype_bed
            # (i.e. cis chunking is also active), we slice *again* here.
            # In 08c usage, cis chunking is NOT active (no --chunk-idx), so
            # phenotype_df at this point contains all phenotypes.
            phenotype_df_trans = phenotype_df.iloc[t_start:t_end]
            log.info(f"Trans chunk {args.trans_chunk_idx}: {len(phenotype_df_trans)} "
                     f"phenotypes (rows {t_start}–{t_end})")
        else:
            # Legacy mode: use full phenotype_df (may already be cis-chunked)
            phenotype_df_trans = phenotype_df
            log.info("=== Trans pass (legacy full-phenotype mode) ===")

        log.info(f"pval threshold: {args.trans_pval_threshold}")
        from tensorqtl import trans
        log.info(f"Trans: batch_size={args.trans_batch_size} | device={'cuda' if torch.cuda.is_available() else 'cpu'}")
        trans_df = trans.map_trans(
            genotype_df,
            phenotype_df_trans,
            covariates_df=cov_df,
            batch_size=args.trans_batch_size,
            return_sparse=True,
            pval_threshold=args.trans_pval_threshold,
        )
        log.info(f"Trans pairs (p < {args.trans_pval_threshold}): {len(trans_df):,}")

        # Post-hoc distance filter: keep only pairs > 1 Mb apart
        # (tensorQTL does not natively filter by distance)
        if len(trans_df) > 0 and "pos" in trans_df.columns and "phenotype_pos" in trans_df.columns:
            trans_far = trans_df[
                (trans_df["pos"] - trans_df["phenotype_pos"]).abs() > 1_000_000
            ].copy()
            log.info(f"Trans pairs >1 Mb: {len(trans_far):,}")
        else:
            trans_far = trans_df.copy()
            log.warning("Could not apply distance filter — missing pos/phenotype_pos columns")

        # Write full sparse output
        trans_df.to_parquet(args.trans_output)
        log.info(f"Trans output: {args.trans_output}")

        # Write >1Mb filtered output alongside
        far_out = args.trans_output.replace(".parquet", "_gt1mb.parquet")
        trans_far.to_parquet(far_out)
        log.info(f"Trans >1Mb output: {far_out} ({len(trans_far):,} pairs)")

    log.info("Done.")


if __name__ == '__main__':
    main()
