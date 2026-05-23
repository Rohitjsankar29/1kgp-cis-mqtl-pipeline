#!/usr/bin/env python3
"""
02_prepare_covariates.py
========================
Prepare covariate matrix for tensorQTL.

Covariates included:
    1. basecaller_model   — fixed effect, categorical (one-hot encoded)
                           Addresses systematic methylation calibration
                           differences across Guppy/Dorado model versions
                           (Decision R4, 2026-02-26)
    2. sex                — inferred from X/Y coverage ratio (cramino stats)
                           or from metadata TSV if available
    3. Genotype PCs       — top N PCs from plink2 PCA output
                           Addresses population stratification
    4. Methylation PCs    — top N PCs from methylation beta matrix
                           Addresses hidden technical/biological structure

Output:
    covariates.tsv  — tensorQTL-format: rows = covariates, cols = samples
                      (transpose of what you'd expect — tensorQTL convention)

Reference:
    Taylor-Weiner et al. (2019) Scaling computational genomics to millions of
    individuals with GPUs. Genome Biology. doi:10.1186/s13059-019-1836-7

Usage:
    python 02_prepare_covariates.py \\
        --samples config/samples.tsv \\
        --methylation-pcs methylation_beta.bed.gz \\
        --geno-pcs plink/cohort_pca.eigenvec \\
        --n-geno-pcs 5 \\
        --n-methyl-pcs 10 \\
        --output covariates.tsv

Author: Kim Navarro (k1mnav)
Date: 2026-02-26
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def load_samples(samples_tsv: Path) -> pd.DataFrame:
    """Load sample manifest (config/samples.tsv)."""
    df = pd.read_csv(samples_tsv, sep="\t")
    # Normalise column names
    df.columns = df.columns.str.lower().str.replace(" ", "_")
    required = {"sample_id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"samples.tsv missing columns: {missing}")
    return df.set_index("sample_id")


def encode_basecaller(samples: pd.DataFrame) -> pd.DataFrame:
    """
    One-hot encode basecaller_model column (drop first to avoid collinearity).
    Returns DataFrame with one column per model level (minus reference level).
    """
    if "basecaller_model" not in samples.columns:
        log.warning("No 'basecaller_model' column in samples.tsv — skipping")
        return pd.DataFrame(index=samples.index)

    dummies = pd.get_dummies(
        samples["basecaller_model"],
        prefix="basecaller",
        drop_first=True,  # reference level absorbed into intercept
    ).astype(int)
    log.info(f"Basecaller covariates: {list(dummies.columns)}")
    return dummies


def encode_sex(samples: pd.DataFrame) -> pd.DataFrame:
    """
    Encode sex as 0/1 (female=0, male=1).
    Looks for 'sex' or 'predicted_sex' column in samples.
    """
    for col in ("predicted_sex", "sex", "gender"):
        if col in samples.columns:
            sex = samples[col].str.lower().map({"female": 0, "male": 1, "f": 0, "m": 1})
            if sex.isna().any():
                log.warning(f"Some sex values could not be coded — will be NaN in covariates")
            log.info(f"Sex covariate from column '{col}'")
            return pd.DataFrame({"sex_male": sex}, index=samples.index)
    log.warning("No sex column found — skipping sex covariate")
    return pd.DataFrame(index=samples.index)


def compute_methylation_pcs(
    methylation_bed: Path,
    n_pcs: int = 10,
    sample_order: list = None,
) -> pd.DataFrame:
    """
    Compute top N PCs from methylation beta matrix.
    Input: tensorQTL phenotype BED (gzipped).
    Returns: DataFrame (samples × PCs).
    """
    log.info(f"Computing {n_pcs} methylation PCs from {methylation_bed.name}")
    df = pd.read_csv(methylation_bed, sep="\t", index_col=3, compression="gzip")
    # Drop coordinate columns
    df = df.drop(columns=["#chr", "start", "end"], errors="ignore")

    # Transpose: samples × sites
    X = df.T
    if sample_order:
        X = X.reindex(sample_order)

    # Drop sites with any NaN after reindex
    X = X.dropna(axis=1)

    # Scale
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # PCA
    pca = PCA(n_components=min(n_pcs, X.shape[0] - 1, X.shape[1]))
    pcs = pca.fit_transform(X_scaled)
    explained = pca.explained_variance_ratio_
    log.info(f"  Top {n_pcs} PCs explain {explained[:n_pcs].sum():.1%} of variance")
    for i, ev in enumerate(explained[:n_pcs]):
        log.info(f"    PC{i+1}: {ev:.3%}")

    cols = [f"methyl_PC{i+1}" for i in range(pcs.shape[1])]
    return pd.DataFrame(pcs, index=X.index, columns=cols)


def load_geno_pcs(eigenvec_path: Path, n_pcs: int = 5) -> pd.DataFrame:
    """
    Load plink2 PCA output (.eigenvec).
    plink2 format: #IID  PC1  PC2 ... (header row, IID = sample ID)
    """
    log.info(f"Loading {n_pcs} genotype PCs from {eigenvec_path.name}")
    df = pd.read_csv(eigenvec_path, sep="\t", index_col=1)
    df = df.drop(columns=["#FID"], errors="ignore")
    df.index.name = "sample_id"
    pc_cols = [c for c in df.columns if c.startswith("PC")][:n_pcs]
    log.info(f"  Using PCs: {pc_cols}")
    return df[pc_cols].rename(columns=lambda c: f"geno_{c}")


def main():
    parser = argparse.ArgumentParser(
        description="Prepare covariate matrix for tensorQTL"
    )
    parser.add_argument("--samples", type=Path, required=True,
                        help="Sample manifest TSV (config/samples.tsv)")
    parser.add_argument("--methylation-pcs", type=Path,
                        help="Methylation beta matrix BED.gz for PC calculation")
    parser.add_argument("--geno-pcs", type=Path,
                        help="plink2 .eigenvec file")
    parser.add_argument("--n-geno-pcs", type=int, default=5)
    parser.add_argument("--n-methyl-pcs", type=int, default=10)
    parser.add_argument("--output", type=Path, required=True,
                        help="Output covariates.tsv (tensorQTL format)")
    args = parser.parse_args()

    samples = load_samples(args.samples)
    sample_ids = list(samples.index)
    log.info(f"Samples: {len(sample_ids)}")

    cov_parts = []

    # 1. Basecaller model
    bc = encode_basecaller(samples)
    if not bc.empty:
        cov_parts.append(bc)

    # 2. Sex
    sex = encode_sex(samples)
    if not sex.empty:
        cov_parts.append(sex)

    # 3. Genotype PCs
    if args.geno_pcs and args.geno_pcs.exists():
        geno_pcs = load_geno_pcs(args.geno_pcs, n_pcs=args.n_geno_pcs)
        cov_parts.append(geno_pcs)
    else:
        log.warning("No genotype PCs provided — omitting")

    # 4. Methylation PCs
    if args.methylation_pcs and args.methylation_pcs.exists():
        methyl_pcs = compute_methylation_pcs(
            args.methylation_pcs,
            n_pcs=args.n_methyl_pcs,
            sample_order=sample_ids,
        )
        cov_parts.append(methyl_pcs)
    else:
        log.warning("No methylation matrix provided — omitting methylation PCs")

    if not cov_parts:
        log.error("No covariates could be built — check inputs")
        sys.exit(1)

    # Combine
    covariates = pd.concat(cov_parts, axis=1)
    covariates = covariates.reindex(sample_ids)
    log.info(f"Covariate matrix: {covariates.shape[0]} samples × {covariates.shape[1]} covariates")
    log.info(f"Covariates: {list(covariates.columns)}")

    # tensorQTL expects rows = covariates, cols = samples (transposed)
    out = covariates.T
    out.index.name = "ID"
    out.to_csv(args.output, sep="\t", float_format="%.6f")
    log.info(f"Written → {args.output}")


if __name__ == "__main__":
    main()
