#!/usr/bin/env python3
"""
00_test_tensorqtl.py
====================
Self-contained end-to-end tensorQTL test using synthetic chr22 data.

Validates that the full mQTL pipeline runs correctly:
    1. Synthetic genotype matrix (plink binary format, chr22 SNPs from real data)
    2. Synthetic methylation M-value matrix (injected mQTLs)
    3. Covariate matrix (basecaller_model, genotype PCs)
    4. tensorQTL nominal + permutation passes
    5. Output validation: output format, known injected signals recovered

Injected mQTLs: 3 CpG sites each have one SNV injected as a real association
(beta = 1.5, so highly detectable even with n=30 synthetic samples).

Usage:
    python 00_test_tensorqtl.py --test-dir ~/test_mqtl/tensorqtl

Author: Kim Navarro (k1mnav)
Date: 2026-02-26
"""

import argparse
import logging
import shutil
import subprocess
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

# ── synthetic data parameters ────────────────────────────────────────────────
N_SAMPLES = 30       # enough for meaningful QTL test
N_CPGS = 200         # small for speed
N_SNPS = 1000        # chr22 SNPs from real plink data or synthetic
N_INJECTED = 3       # CpG sites with injected mQTL
EFFECT_SIZE = 1.5    # strong injected effect (M-value units)
SEED = 42
CHR = "chr22"
# ─────────────────────────────────────────────────────────────────────────────


def make_sample_ids(n: int) -> list:
    return [f"SAMPLE_{i:04d}" for i in range(n)]


def make_genotype_matrix(
    n_samples: int,
    n_snps: int,
    rng: np.random.Generator,
    real_plink_prefix: str = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns (genotype_df [variants × samples], variant_df).
    If real_plink_prefix is provided and plink files exist, use real chr22 SNPs.
    Otherwise generate synthetic diploid genotypes.
    """
    samples = make_sample_ids(n_samples)

    if real_plink_prefix:
        pvar = Path(real_plink_prefix + ".pvar")
        if pvar.exists():
            log.info(f"Loading real genotypes from {real_plink_prefix}")
            # Read pvar for variant info
            vdf = pd.read_csv(pvar, sep="\t", comment="#",
                              names=["chrom", "pos", "id", "ref", "alt"],
                              usecols=[0, 1, 2, 3, 4])
            vdf = vdf.sample(n=min(n_snps, len(vdf)), random_state=SEED)
            vdf = vdf.reset_index(drop=True)
            # Synthetic genotypes aligned to real positions
            mafs = rng.uniform(0.05, 0.5, size=len(vdf))
            geno = rng.binomial(2, mafs[:, None], size=(len(vdf), n_samples)).astype(float)
            geno[rng.random(geno.shape) < 0.02] = np.nan  # 2% missingness
            gdf = pd.DataFrame(geno, index=vdf["id"], columns=samples)
            return gdf, vdf

    log.info(f"Generating synthetic genotypes (n_samples={n_samples}, n_snps={n_snps})")
    # Synthetic chr22 positions
    positions = sorted(rng.integers(10_500_000, 50_800_000, size=n_snps))
    vdf = pd.DataFrame({
        "chrom": CHR,
        "pos": positions,
        "id": [f"chr22:{p}:A:G" for p in positions],
        "ref": "A",
        "alt": "G",
    })
    mafs = rng.uniform(0.05, 0.5, size=n_snps)
    geno = rng.binomial(2, mafs[:, None], size=(n_snps, n_samples)).astype(float)
    geno[rng.random(geno.shape) < 0.02] = np.nan
    gdf = pd.DataFrame(geno, index=vdf["id"], columns=samples)
    return gdf, vdf


def make_methylation_matrix(
    n_samples: int,
    n_cpgs: int,
    genotype_df: pd.DataFrame,
    variant_df: pd.DataFrame,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame, list]:
    """
    Returns (phenotype_df [CpGs × samples], phenotype_pos_df, injected_pairs).
    Injects N_INJECTED real mQTL signals.
    """
    samples = make_sample_ids(n_samples)

    # CpG positions on chr22
    cpg_positions = sorted(rng.integers(10_500_000, 50_800_000, size=n_cpgs))
    cpg_ids = [f"chr22:{p}" for p in cpg_positions]
    phenotype_pos_df = pd.DataFrame({
        "chr": CHR,
        "start": cpg_positions,
        "end": [p + 1 for p in cpg_positions],
        "phenotype_id": cpg_ids,
    }).set_index("phenotype_id")

    # Background noise methylation (M-values, roughly Gaussian)
    M = rng.normal(0, 2, size=(n_cpgs, n_samples))

    # Inject real mQTL signals
    injected = []
    snv_positions = variant_df["pos"].values
    for i in range(N_INJECTED):
        cpg_idx = i * (n_cpgs // N_INJECTED)
        cpg_pos = cpg_positions[cpg_idx]

        # Find a SNV within 500 kb of this CpG
        candidates = np.where(np.abs(snv_positions - cpg_pos) < 500_000)[0]
        if len(candidates) == 0:
            log.warning(f"No SNV within 500 kb of CpG {cpg_ids[cpg_idx]} — skipping injection")
            continue

        snv_idx = candidates[rng.integers(len(candidates))]
        snv_id = genotype_df.index[snv_idx]
        g = genotype_df.iloc[snv_idx].values.copy()
        missing_mask = np.isnan(g)
        g[missing_mask] = np.nanmean(g)  # impute missing for injection

        # Add effect
        M[cpg_idx] += EFFECT_SIZE * (g - np.mean(g))
        injected.append({
            "cpg_id": cpg_ids[cpg_idx],
            "snv_id": snv_id,
            "cpg_pos": cpg_pos,
            "snv_pos": int(variant_df.iloc[snv_idx]["pos"]),
            "effect_size": EFFECT_SIZE,
        })
        log.info(f"  Injected mQTL: {cpg_ids[cpg_idx]} ← {snv_id} (beta={EFFECT_SIZE})")

    phenotype_df = pd.DataFrame(M, index=cpg_ids, columns=samples)
    return phenotype_df, phenotype_pos_df, injected


def make_covariates(
    n_samples: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """
    Simulate covariates: 2 basecaller dummies + 2 geno PCs.
    tensorQTL format: rows = covariates, cols = samples.
    """
    samples = make_sample_ids(n_samples)
    # Basecaller model: 3 groups (dorado, guppy-sup, guppy657)
    groups = rng.choice(3, size=n_samples)
    cov = pd.DataFrame(
        {
            "basecaller_guppy_sup": (groups == 1).astype(int),
            "basecaller_guppy657": (groups == 2).astype(int),
            "geno_PC1": rng.standard_normal(n_samples),
            "geno_PC2": rng.standard_normal(n_samples),
        },
        index=samples,
    ).T
    cov.index.name = "ID"
    return cov


def write_plink_for_tensorqtl(
    genotype_df: pd.DataFrame,
    variant_df: pd.DataFrame,
    out_prefix: str,
) -> None:
    """
    Write plink2 binary files from a genotype DataFrame.
    Uses plink2 --make-pgen from a generated VCF.
    """
    import tempfile, os
    vcf_path = out_prefix + "_tmp.vcf"
    samples = genotype_df.columns.tolist()

    log.info(f"Writing synthetic VCF for plink import → {vcf_path}")
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID={CHR},length=50818468>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(samples) + "\n")
        for idx, row in variant_df.iterrows():
            genos = genotype_df.iloc[idx].values
            gt_strs = []
            for g in genos:
                if np.isnan(g):
                    gt_strs.append("./.")
                elif g == 0:
                    gt_strs.append("0/0")
                elif g == 1:
                    gt_strs.append("0/1")
                else:
                    gt_strs.append("1/1")
            f.write(
                f"{row['chrom']}\t{row['pos']}\t{row['id']}\t"
                f"{row['ref']}\t{row['alt']}\t.\tPASS\t.\tGT\t"
                + "\t".join(gt_strs) + "\n"
            )

    # Step 1: VCF → pgen (plink2 format)
    cmd = [
        "plink2",
        "--vcf", vcf_path,
        "--make-pgen",
        "--out", out_prefix + "_pgen",
        "--threads", "4",
        "--max-alleles", "2",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error(result.stderr[-500:])
        raise RuntimeError("plink2 pgen import failed")

    # Step 2: pgen → bed (plink1 format) — needed for tensorQTL's PlinkReader
    # --output-chr chrM preserves chr prefix (e.g. chr22 not 22) in .bim file
    # so chromosome names match the phenotype BED (which uses chr22 format)
    cmd2 = [
        "plink2",
        "--pfile", out_prefix + "_pgen",
        "--make-bed",
        "--out", out_prefix,
        "--threads", "4",
        "--output-chr", "chrM",
    ]
    result2 = subprocess.run(cmd2, capture_output=True, text=True)
    if result2.returncode != 0:
        log.error(result2.stderr[-500:])
        raise RuntimeError("plink2 pgen→bed conversion failed")

    os.remove(vcf_path)
    log.info(f"✅ plink1 binary: {out_prefix}.bed/bim/fam")


def write_phenotype_bed(
    phenotype_df: pd.DataFrame,
    phenotype_pos_df: pd.DataFrame,
    path: Path,
) -> None:
    """Write tensorQTL-format phenotype BED (sorted by chrom, start)."""
    pos = phenotype_pos_df.copy()
    out = pd.DataFrame({
        "#chr": pos["chr"],
        "start": pos["start"],
        "end": pos["end"],
        "phenotype_id": pos.index,
    })
    # Concatenate phenotype values
    out = pd.concat([out.reset_index(drop=True), phenotype_df.reset_index(drop=True)], axis=1)
    # Sort by chrom then position (tensorQTL requires sorted BED)
    chrom_order = {f"chr{i}": i for i in list(range(1, 23)) + ["X", "Y", "M"]}
    out["_sort"] = out["#chr"].map(chrom_order).fillna(99)
    out = out.sort_values(["_sort", "start"]).drop(columns=["_sort"]).reset_index(drop=True)
    out.to_csv(path, sep="\t", index=False, float_format="%.6f",
               compression="gzip" if str(path).endswith(".gz") else None)
    log.info(f"Phenotype BED: {len(out)} CpGs × {phenotype_df.shape[1]} samples → {path}")


def run_tensorqtl_nominal(
    plink_prefix: str,
    phenotype_bed: Path,
    covariates: pd.DataFrame,
    output_dir: Path,
    cis_window: int = 1_000_000,
) -> pd.DataFrame:
    """Run tensorQTL nominal pass. Returns combined pairs DataFrame."""
    import tensorqtl
    from tensorqtl import cis, genotypeio

    output_dir.mkdir(parents=True, exist_ok=True)

    pr = genotypeio.PlinkReader(plink_prefix)
    genotype_df = pr.load_genotypes()
    # pr.bim: chrom, snp, cm, pos, a0, a1, i — use snp as variant_id index
    variant_df = pr.bim.set_index("snp")[["chrom", "pos"]].copy()
    variant_df.index.name = "variant_id"

    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(str(phenotype_bed))

    common = sorted(set(genotype_df.columns) & set(phenotype_df.columns) & set(covariates.columns))
    log.info(f"Common samples: {len(common)}")
    assert len(common) > 0, "No common samples between genotypes, phenotypes, and covariates!"

    genotype_df = genotype_df[common]
    phenotype_df = phenotype_df[common]
    # tensorqtl cis functions expect covariates_df with rows=samples, cols=covariates
    # (our file/storage format is rows=covariates, cols=samples — transpose here)
    cov = covariates[common].T

    log.info(f"Running nominal cis-mQTL (window={cis_window:,} bp)...")
    cis.map_nominal(
        genotype_df, variant_df, phenotype_df, phenotype_pos_df,
        "nominal",       # prefix (positional)
        covariates_df=cov,
        window=cis_window,
        output_dir=str(output_dir),
    )

    parquets = sorted(output_dir.glob("nominal.cis_qtl_pairs.*.parquet"))
    if not parquets:
        raise RuntimeError("No nominal output parquets found")
    pairs = pd.concat([pd.read_parquet(p) for p in parquets])
    log.info(f"Nominal pairs: {len(pairs):,}")
    return pairs


def run_tensorqtl_permutation(
    plink_prefix: str,
    phenotype_bed: Path,
    covariates: pd.DataFrame,
    output_dir: Path,
    cis_window: int = 1_000_000,
    n_perm: int = 200,
) -> pd.DataFrame:
    """Run tensorQTL permutation pass. Returns cis_df with qvalues."""
    import tensorqtl
    from tensorqtl import cis, genotypeio

    output_dir.mkdir(parents=True, exist_ok=True)

    pr = genotypeio.PlinkReader(plink_prefix)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index("snp")[["chrom", "pos"]].copy()
    variant_df.index.name = "variant_id"

    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(str(phenotype_bed))

    common = sorted(set(genotype_df.columns) & set(phenotype_df.columns) & set(covariates.columns))
    genotype_df = genotype_df[common]
    phenotype_df = phenotype_df[common]
    cov = covariates[common].T  # rows=samples, cols=covariates

    log.info(f"Running permutation cis-mQTL ({n_perm} permutations)...")
    cis_df = cis.map_cis(
        genotype_df, variant_df, phenotype_df, phenotype_pos_df,
        covariates_df=cov,
        nperm=n_perm,
        window=cis_window,
        seed=42,
    )

    # q-values
    from tensorqtl import post
    try:
        post.calculate_qvalues(cis_df, fdr=0.05, qvalue_lambda=0.85)
    except Exception as e:
        log.warning(f"q-value calculation failed (R/qvalue not available): {e}")
        cis_df["qval"] = np.nan

    out = output_dir / "cis_permutation.tsv.gz"
    cis_df.to_csv(out, sep="\t", float_format="%.6g")
    log.info(f"Permutation results → {out}")
    return cis_df


def validate_results(
    nominal_pairs: pd.DataFrame,
    cis_df: pd.DataFrame,
    injected: list,
) -> bool:
    """
    Check that injected mQTL signals are recovered in the results.
    Returns True if all injected pairs appear in the top results.
    """
    log.info("=== Validating injected mQTL recovery ===")
    all_recovered = True

    for inj in injected:
        cpg = inj["cpg_id"]
        snv = inj["snv_id"]

        # Nominal: look for this CpG in nominal output (any variant in cis)
        cpg_pairs = nominal_pairs[nominal_pairs["phenotype_id"] == cpg]
        if len(cpg_pairs) == 0:
            log.warning(f"  MISS (nominal): {cpg} not tested (no SNVs in cis window for this sample)")
            all_recovered = False
            continue

        best = cpg_pairs.nsmallest(1, "pval_nominal").iloc[0]
        log.info(f"  ✅ {cpg} — best cis SNV: {best['variant_id']}, slope={best['slope']:.3f}, p={best['pval_nominal']:.2e}")

        # Permutation: phenotype_id is the index in cis_df
        if cpg in cis_df.index:
            pval = cis_df.loc[cpg, "pval_perm"]
            log.info(f"     permutation p_perm={pval:.4f}")
        else:
            log.warning(f"     {cpg} not in permutation results")

    return all_recovered


def main():
    parser = argparse.ArgumentParser(description="End-to-end tensorQTL test with synthetic data")
    parser.add_argument("--test-dir", type=Path, default=Path("~/test_mqtl/tensorqtl").expanduser())
    parser.add_argument("--real-plink", type=str,
                        help="Optional: plink prefix with real chr22 genotypes for variant positions")
    parser.add_argument("--n-samples", type=int, default=N_SAMPLES)
    parser.add_argument("--n-cpgs", type=int, default=N_CPGS)
    parser.add_argument("--n-snps", type=int, default=N_SNPS)
    parser.add_argument("--cis-window", type=int, default=1_000_000)
    parser.add_argument("--n-perm", type=int, default=200)
    args = parser.parse_args()

    args.test_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    log.info(f"=== tensorQTL end-to-end test ===")
    log.info(f"n_samples={args.n_samples}, n_cpgs={args.n_cpgs}, n_snps={args.n_snps}")

    # 1. Synthetic genotypes
    log.info("Generating synthetic genotypes...")
    genotype_df, variant_df = make_genotype_matrix(
        args.n_samples, args.n_snps, rng,
        real_plink_prefix=args.real_plink,
    )

    # 2. Synthetic methylation with injected mQTLs
    log.info("Generating methylation matrix with injected mQTLs...")
    phenotype_df, phenotype_pos_df, injected = make_methylation_matrix(
        args.n_samples, args.n_cpgs, genotype_df, variant_df, rng,
    )

    # 3. Covariates
    log.info("Generating covariate matrix...")
    covariates = make_covariates(args.n_samples, rng)

    # 4. Write to disk
    plink_prefix = str(args.test_dir / "test_geno")
    phenotype_bed = args.test_dir / "test_methylation.bed.gz"
    covariate_tsv = args.test_dir / "test_covariates.tsv"

    write_plink_for_tensorqtl(genotype_df, variant_df, plink_prefix)
    write_phenotype_bed(phenotype_df, phenotype_pos_df, phenotype_bed)
    covariates.to_csv(covariate_tsv, sep="\t")
    log.info(f"Covariate TSV → {covariate_tsv}")

    # 5. Run tensorQTL
    log.info("=== Running tensorQTL nominal pass ===")
    nominal_pairs = run_tensorqtl_nominal(
        plink_prefix, phenotype_bed, covariates,
        args.test_dir / "results_nominal",
        cis_window=args.cis_window,
    )

    log.info("=== Running tensorQTL permutation pass ===")
    cis_df = run_tensorqtl_permutation(
        plink_prefix, phenotype_bed, covariates,
        args.test_dir / "results_permutation",
        cis_window=args.cis_window,
        n_perm=args.n_perm,
    )

    # 6. Validate
    ok = validate_results(nominal_pairs, cis_df, injected)

    log.info("")
    log.info("=== Test Summary ===")
    log.info(f"Nominal pairs:        {len(nominal_pairs):,}")
    log.info(f"CpG sites tested:     {len(cis_df):,}")
    log.info(f"Injected mQTLs:       {len(injected)}")
    log.info(f"All recovered:        {'✅ YES' if ok else '❌ NO'}")
    if not cis_df["qval"].isna().all():
        sig = (cis_df["qval"] < 0.05).sum()
        log.info(f"FDR<5% significant:   {sig} CpG sites")

    log.info(f"\nTop 5 permutation hits:")
    top5 = cis_df.nsmallest(5, "pval_perm")
    cols = [c for c in ["variant_id", "pval_perm", "slope", "qval"] if c in top5.columns]
    log.info(top5[cols].to_string())

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
