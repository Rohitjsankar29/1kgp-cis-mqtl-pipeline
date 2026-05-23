#!/usr/bin/env python3
"""
04_build_covariates.py — Build tensorQTL covariate matrix

Combines:
  - Genotype PCs (from PLINK2 eigenvec)
  - Methylation PCs (from 03_aggregate_pilot.py)
  - Sex (from samples_pilot.tsv)
  - Basecaller model (from samples_pilot.tsv)

Output format: rows=covariates, cols=samples (Decision R7)
tensorQTL expects: cov[samples].T  — handled in run script
"""

import argparse
import pandas as pd
import sys
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--samples-file", required=True)
    p.add_argument("--geno-pca", required=True, help="PLINK2 .eigenvec file")
    p.add_argument("--meth-pca", required=True, help="MethPC TSV (rows=PCs, cols=samples)")
    p.add_argument("--n-geno-pcs", type=int, default=3)
    p.add_argument("--n-meth-pcs", type=int, default=3)
    p.add_argument("--out", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    # --- Sample manifest ---
    samples_df = pd.read_csv(args.samples_file, sep='\t')
    sample_ids = samples_df['sample_id'].tolist()
    log.info(f"Samples: {sample_ids}")

    rows = {}

    # --- Genotype PCs (PLINK2 eigenvec: #FID IID PC1 PC2 ...) ---
    geno_pca = pd.read_csv(args.geno_pca, sep=r'\s+')
    # Column names: #FID, IID, PC1, PC2, ...
    geno_pca = geno_pca.rename(columns={'#FID': 'FID'})
    geno_pca = geno_pca.set_index('IID')

    for i in range(1, args.n_geno_pcs + 1):
        pc_col = f"PC{i}"
        if pc_col not in geno_pca.columns:
            log.warning(f"Genotype {pc_col} not found — skipping")
            continue
        row = {}
        for sid in sample_ids:
            if sid not in geno_pca.index:
                log.error(f"Sample {sid} missing from genotype PCA")
                sys.exit(1)
            row[sid] = geno_pca.loc[sid, pc_col]
        rows[f"GenoPC{i}"] = row
        log.info(f"GenoPC{i}: OK")

    # --- Methylation PCs (rows=PCs, cols=samples) ---
    meth_pca = pd.read_csv(args.meth_pca, sep='\t', index_col=0)
    for i in range(1, args.n_meth_pcs + 1):
        row_name = f"MethPC{i}"
        if row_name not in meth_pca.index:
            log.warning(f"{row_name} not found in methylation PCA — skipping")
            continue
        row = {}
        for sid in sample_ids:
            if sid not in meth_pca.columns:
                log.error(f"Sample {sid} missing from methylation PCA")
                sys.exit(1)
            row[sid] = meth_pca.loc[row_name, sid]
        rows[row_name] = row
        log.info(f"{row_name}: OK")

    # --- Sex covariate (0=XX, 1=XY) ---
    sex_map = {row['sample_id']: (1 if row['sex'] == 'XY' else 0)
               for _, row in samples_df.iterrows()}
    rows['sex'] = {sid: sex_map[sid] for sid in sample_ids}
    log.info("sex: OK")

    # --- Basecaller model (0=R9, 1=R10) ---
    bc_map = {row['sample_id']: row['basecaller_model']
              for _, row in samples_df.iterrows()}
    rows['basecaller_model'] = {sid: bc_map[sid] for sid in sample_ids}
    log.info("basecaller_model: OK")

    # --- Build DataFrame (rows=covariates, cols=samples) ---
    cov_df = pd.DataFrame(rows).T  # shape: n_covariates × n_samples
    cov_df = cov_df[sample_ids]    # ensure column order matches samples

    log.info(f"Covariate matrix: {cov_df.shape[0]} covariates × {cov_df.shape[1]} samples")
    log.info(f"Covariates: {cov_df.index.tolist()}")

    cov_df.to_csv(args.out, sep='\t', float_format='%.6f')
    log.info(f"Written: {args.out}")


if __name__ == '__main__':
    main()
