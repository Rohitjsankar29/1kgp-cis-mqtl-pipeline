#!/usr/bin/env Rscript
# 04_run_coloc_susie_ancestry.R — re-run regional coloc.susie on the SMR+HEIDI
# hits using an ancestry-restricted LD matrix, so results can be compared
# locus-by-locus against the primary multi-ancestry coloc.
#
# Assumes the primary coloc workflow (scripts/coloc/run_coloc.R) has already:
#   1) built per-locus mQTL summary-stat frames and per-locus GWAS frames;
#   2) fit SuSiE credible sets per trait per locus.
# This script re-fits the SuSiE step using the ancestry-restricted LD (computed
# in-script from the per-superpop PLINK panel) and re-invokes coloc.susie() to
# produce a parallel per-locus PP4 (and PP-vector) that can be joined with the
# primary output.
#
# Cost: coloc.susie per locus ~30-60 CPU-min; 90 loci = 45-90 CPU-hours per
# ancestry. On normalq (24 CPU) that's ~2-4 hours wall-time. ~50-90 SU per
# ancestry.
#
# Usage:
#   Rscript 04_run_coloc_susie_ancestry.R --pop EUR --loci egfr_sig_loci_hg19.tsv
#   Rscript 04_run_coloc_susie_ancestry.R --pop EUR --loci cvd_smr_heidi_hits.tsv
#   Rscript 04_run_coloc_susie_ancestry.R --pop EUR --loci t2d_smr_heidi_hits.tsv

suppressPackageStartupMessages({
    library(argparse)
    library(coloc)
    library(susieR)
    library(data.table)
    library(snpStats)      # for reading PLINK BED via read.plink()
})

ap <- ArgumentParser()
ap$add_argument("--pop",  required = TRUE, help = "superpopulation code: EUR/AFR/SAS/EAS/AMR")
ap$add_argument("--loci", required = TRUE, help = "TSV of SMR+HEIDI-significant loci to re-analyse")
ap$add_argument("--gwas-ma",  default = NULL, help = "GWAS .ma file (inferred from --loci name if omitted)")
ap$add_argument("--out-dir",  default = NULL)
ap$add_argument("--window",   type = "integer", default = 1000000)
a <- ap$parse_args()

G   <- "/g/data/cy94/rs4477/downstream/genome_452"
POP <- a$pop
BFILE <- file.path(G, "ancestry/plink", POP)
if (!file.exists(paste0(BFILE, ".bed"))) stop(sprintf("MISSING: %s.bed — run step 02 first", BFILE))

# infer trait from loci filename if the GWAS wasn't passed
if (is.null(a$gwas_ma)) {
    trait <- sub(".*(t2d|cvd|egfr|cad).*", "\\1", tolower(basename(a$loci)))
    trait <- ifelse(trait == "cvd", "cad", trait)
    a$gwas_ma <- file.path(G, "gwas", sprintf("%s.ma", trait))
}
if (is.null(a$out_dir)) {
    trait <- sub(".*(t2d|cvd|egfr|cad).*", "\\1", tolower(basename(a$loci)))
    trait <- ifelse(trait == "cvd", "cad", trait)
    a$out_dir <- file.path(G, "ancestry/coloc", POP, trait)
}
dir.create(a$out_dir, recursive = TRUE, showWarnings = FALSE)

loci   <- fread(a$loci)
gwas   <- fread(a$gwas_ma)
setnames(gwas, tolower(names(gwas)))
setkey(gwas, snp)

# read the ancestry panel once; subset per locus below
message(sprintf("[%s] reading %s panel", format(Sys.time(), "%H:%M:%S"), POP))
plink_all <- read.plink(paste0(BFILE, ".bed"), paste0(BFILE, ".bim"), paste0(BFILE, ".fam"))
bim <- as.data.table(plink_all$map)
setnames(bim, c("chr", "snp", "cm", "bp", "a1", "a2"))
G_mat <- as(plink_all$genotypes, "numeric")   # sample x SNP (0/1/2)

results <- list()
for (i in seq_len(nrow(loci))) {
    row <- loci[i]
    chr <- as.integer(row$chr); pos <- as.integer(row$bp %||% row$pos %||% row$topSNP_bp)
    if (is.na(chr) || is.na(pos)) next

    keep <- bim[chr == !!chr & bp >= pos - a$window & bp <= pos + a$window, snp]
    if (length(keep) < 20) next

    g_win <- G_mat[, keep, drop = FALSE]
    g_win[is.na(g_win)] <- matrix(colMeans(g_win, na.rm = TRUE),
                                  nrow = nrow(g_win), ncol = ncol(g_win), byrow = TRUE)[is.na(g_win)]
    LD <- cor(g_win, use = "pairwise.complete.obs")

    # trait-side SuSiE with the ancestry LD; the mQTL side already has its own
    # SuSiE fit from the primary run — reload it via run_coloc.R's helper.
    gwas_win <- gwas[snp %in% keep, .(snp, beta = as.numeric(b), varbeta = as.numeric(se)^2)]
    setkey(gwas_win, snp); gwas_win <- gwas_win[keep, nomatch = 0]
    if (nrow(gwas_win) < 20) next

    s_trait <- tryCatch(
        susie_rss(z = gwas_win$beta / sqrt(gwas_win$varbeta), R = LD, n = row$gwas_n %||% 100000, L = 10),
        error = function(e) { message(sprintf("locus %d: SuSiE failed (%s)", i, conditionMessage(e))); NULL })
    if (is.null(s_trait)) next

    # primary mQTL SuSiE fit reloaded from the flagship coloc outputs
    mqtl_susie_file <- file.path(G, "coloc/susie_fits", sprintf("%s_susie.rds", row$topSNP %||% row$snp))
    if (!file.exists(mqtl_susie_file)) next
    s_mqtl <- readRDS(mqtl_susie_file)

    co <- tryCatch(coloc.susie(s_mqtl, s_trait), error = function(e) NULL)
    if (is.null(co) || is.null(co$summary)) next
    co$summary$locus <- row$topSNP %||% row$snp
    co$summary$pop   <- POP
    results[[length(results) + 1]] <- co$summary
}

out <- rbindlist(results, fill = TRUE)
outf <- file.path(a$out_dir, sprintf("coloc_susie_%s.tsv", POP))
fwrite(out, outf, sep = "\t")
message(sprintf("[%s] DONE — %d loci → %s", format(Sys.time(), "%H:%M:%S"), nrow(out), outf))
