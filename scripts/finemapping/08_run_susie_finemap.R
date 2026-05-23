#!/usr/bin/env Rscript
# 08_run_susie_finemap.R
# Fine-map cis-mQTL loci using SuSiE.
#
# Inputs:
#   --qtl          TensorQTL nominal/permutation table with variant_id, phenotype_id, pval/beta
#   --geno         genotype dosage matrix for one chromosome/locus (variants x samples)
#   --pheno        phenotype matrix for one chromosome/locus (CpGs x samples)
#   --loci         optional loci table: phenotype_id, chr, start, end
#   --out          output credible set table
#
# Expected use:
#   Run after TensorQTL. For each CpG, take variants in the cis-window and use SuSiE
#   to estimate posterior inclusion probabilities (PIP) and credible sets.
#
# Notes:
#   - This is a practical pilot script. For full production, run chromosome/locus chunks.
#   - Genotype and phenotype sample columns must match.
#   - If qtl has effect sizes only, this script refits per-locus regression through susieR.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(susieR)
})

option_list <- list(
  make_option("--qtl", type="character", help="TensorQTL nominal results TSV/CSV/parquet-export table"),
  make_option("--geno", type="character", help="Genotype dosage matrix TSV: variant_id + sample columns"),
  make_option("--pheno", type="character", help="Methylation phenotype matrix TSV: phenotype_id + sample columns"),
  make_option("--loci", type="character", default=NULL, help="Optional loci TSV: phenotype_id chr start end"),
  make_option("--out", type="character", help="Output TSV for SuSiE results"),
  make_option("--max-variants", type="integer", default=5000, help="Max variants per CpG locus"),
  make_option("--L", type="integer", default=10, help="Max number of causal effects for SuSiE"),
  make_option("--min-pip", type="double", default=0.01, help="Only write variants with PIP >= this value")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$qtl) || is.null(opt$geno) || is.null(opt$pheno) || is.null(opt$out)) {
  stop("Required: --qtl --geno --pheno --out")
}

message("[", Sys.time(), "] Loading inputs")
qtl <- fread(opt$qtl)
geno <- fread(opt$geno)
pheno <- fread(opt$pheno)

if (!("variant_id" %in% names(geno))) stop("Genotype matrix must contain variant_id column")
if (!("phenotype_id" %in% names(pheno))) stop("Phenotype matrix must contain phenotype_id column")
if (!all(c("variant_id", "phenotype_id") %in% names(qtl))) {
  stop("QTL table must contain variant_id and phenotype_id")
}

geno_mat <- as.data.frame(geno)
rownames(geno_mat) <- geno_mat$variant_id
geno_mat$variant_id <- NULL

pheno_mat <- as.data.frame(pheno)
rownames(pheno_mat) <- pheno_mat$phenotype_id
pheno_mat$phenotype_id <- NULL

common_samples <- intersect(colnames(geno_mat), colnames(pheno_mat))
if (length(common_samples) < 10) stop("Too few matching samples between genotype and phenotype matrices")

geno_mat <- geno_mat[, common_samples, drop=FALSE]
pheno_mat <- pheno_mat[, common_samples, drop=FALSE]

message("[", Sys.time(), "] Common samples: ", length(common_samples))

phenotypes <- unique(qtl$phenotype_id)
results <- list()
k <- 1

for (pid in phenotypes) {
  if (!(pid %in% rownames(pheno_mat))) next

  vars <- unique(qtl[phenotype_id == pid]$variant_id)
  vars <- vars[vars %in% rownames(geno_mat)]
  if (length(vars) < 2) next

  if (length(vars) > opt$max_variants) {
    if ("pval_nominal" %in% names(qtl)) {
      top <- qtl[phenotype_id == pid & variant_id %in% vars][order(pval_nominal)][1:opt$max_variants]
      vars <- top$variant_id
    } else {
      vars <- vars[1:opt$max_variants]
    }
  }

  X <- t(as.matrix(geno_mat[vars, , drop=FALSE]))
  y <- as.numeric(pheno_mat[pid, common_samples])

  keep <- complete.cases(y) & apply(X, 1, function(z) all(is.finite(z)))
  X <- X[keep, , drop=FALSE]
  y <- y[keep]

  if (nrow(X) < 10 || ncol(X) < 2) next

  # Remove monomorphic variants
  sds <- apply(X, 2, sd, na.rm=TRUE)
  X <- X[, sds > 0, drop=FALSE]
  if (ncol(X) < 2) next

  # Standardise predictors and response
  X <- scale(X)
  y <- as.numeric(scale(y))

  fit <- tryCatch(
    susie(X, y, L=opt$L, standardize=FALSE, verbose=FALSE),
    error=function(e) {
      message("SuSiE failed for ", pid, ": ", e$message)
      NULL
    }
  )
  if (is.null(fit)) next

  pip <- susie_get_pip(fit)
  cs <- susie_get_cs(fit, X=X, coverage=0.95)

  out <- data.table(
    phenotype_id = pid,
    variant_id = colnames(X),
    pip = pip
  )
  out <- out[pip >= opt$min_pip]

  # Add credible-set labels where available
  out[, credible_set := NA_character_]
  if (!is.null(cs$cs)) {
    for (cs_name in names(cs$cs)) {
      idx <- cs$cs[[cs_name]]
      cs_vars <- colnames(X)[idx]
      out[variant_id %in% cs_vars, credible_set := cs_name]
    }
  }

  if (nrow(out) > 0) {
    results[[k]] <- out
    k <- k + 1
  }

  if (k %% 50 == 0) message("[", Sys.time(), "] Processed ", k, " loci")
}

if (length(results) == 0) {
  warning("No SuSiE results produced")
  fwrite(data.table(phenotype_id=character(), variant_id=character(), pip=numeric(), credible_set=character()), opt$out, sep="\t")
} else {
  res <- rbindlist(results, fill=TRUE)
  fwrite(res, opt$out, sep="\t")
  message("[", Sys.time(), "] Written: ", opt$out, " (", nrow(res), " rows)")
}
