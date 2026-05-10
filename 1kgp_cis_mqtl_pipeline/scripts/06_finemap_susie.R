</> R

#!/usr/bin/env Rscript

# Fine-map significant cis-mQTL loci using SuSiE.
# Input: TensorQTL nominal results + genotype matrix + CpG methylation matrix.
# Output: variant-CpG fine-mapping table with PIP and credible set info.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(susieR)
})

option_list <- list(
  make_option("--nominal", type="character", help="TensorQTL nominal results TSV"),
  make_option("--phenotype", type="character", help="CpG methylation matrix TSV"),
  make_option("--genotype", type="character", help="Variant genotype dosage matrix TSV"),
  make_option("--out", type="character", help="Output fine-mapping TSV"),
  make_option("--pval-threshold", type="double", default=1e-5),
  make_option("--max-loci", type="integer", default=1000),
  make_option("--L", type="integer", default=10)
)

args <- parse_args(OptionParser(option_list=option_list))

stopifnot(file.exists(args$nominal))
stopifnot(file.exists(args$phenotype))
stopifnot(file.exists(args$genotype))

message("Loading TensorQTL nominal results...")
nom <- fread(args$nominal)

required_cols <- c("phenotype_id", "variant_id", "pval_nominal")
missing_cols <- setdiff(required_cols, names(nom))
if (length(missing_cols) > 0) {
  stop("Missing columns in nominal file: ", paste(missing_cols, collapse=", "))
}

nom <- nom[pval_nominal <= args$pval_threshold]
if (nrow(nom) == 0) stop("No associations pass p-value threshold.")

message("Loading methylation phenotype matrix...")
pheno <- fread(args$phenotype)
pheno_ids <- pheno[[4]]
pheno_mat <- as.matrix(pheno[, -(1:4), with=FALSE])
rownames(pheno_mat) <- pheno_ids

message("Loading genotype dosage matrix...")
geno <- fread(args$genotype)
variant_ids <- geno[[1]]
geno_mat <- as.matrix(geno[, -1, with=FALSE])
rownames(geno_mat) <- variant_ids

common_samples <- intersect(colnames(pheno_mat), colnames(geno_mat))
if (length(common_samples) < 10) stop("Too few overlapping samples.")

pheno_mat <- pheno_mat[, common_samples, drop=FALSE]
geno_mat <- geno_mat[, common_samples, drop=FALSE]

cpgs <- unique(nom$phenotype_id)
cpgs <- cpgs[seq_len(min(length(cpgs), args$max_loci))]

results <- list()

for (cpg in cpgs) {
  message("Fine-mapping: ", cpg)

  locus <- nom[phenotype_id == cpg]
  vars <- intersect(locus$variant_id, rownames(geno_mat))

  if (length(vars) < 2) next
  if (!(cpg %in% rownames(pheno_mat))) next

  X <- t(geno_mat[vars, , drop=FALSE])
  y <- as.numeric(pheno_mat[cpg, ])

  keep <- complete.cases(y) & complete.cases(X)
  X <- X[keep, , drop=FALSE]
  y <- y[keep]

  if (nrow(X) < 10 || ncol(X) < 2) next

  X <- scale(X)
  y <- scale(y)[, 1]

  fit <- tryCatch(
    susie(X, y, L=args$L, estimate_residual_variance=TRUE),
    error=function(e) NULL
  )

  if (is.null(fit)) next

  pip <- susie_get_pip(fit)

  cs <- rep(NA, length(vars))
  if (!is.null(fit$sets$cs)) {
    for (i in seq_along(fit$sets$cs)) {
      cs[fit$sets$cs[[i]]] <- paste0("CS", i)
    }
  }

  res <- data.table(
    phenotype_id = cpg,
    variant_id = vars,
    pip = pip,
    credible_set = cs
  )

  results[[cpg]] <- res
}

if (length(results) == 0) stop("No loci were successfully fine-mapped.")

out <- rbindlist(results, fill=TRUE)
fwrite(out, args$out, sep="\t")

message("Written: ", args$out)
