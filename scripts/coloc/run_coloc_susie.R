suppressMessages({library(coloc); library(susieR)})
args <- commandArgs(trailingOnly=TRUE)
# args: mqtl_file gwas_file ld_file out_prefix
d <- read.table(args[1], header=TRUE, sep="\t", stringsAsFactors=FALSE)
# d columns: snp, mqtl_beta, mqtl_se, mqtl_maf, mqtl_n, gwas_beta, gwas_se, gwas_maf, gwas_n
ld <- as.matrix(read.table(args[2]))
rownames(ld) <- colnames(ld) <- d$snp
n_mqtl <- d$mqtl_n[1]; n_gwas <- max(d$gwas_n)
D1 <- list(beta=d$mqtl_beta, varbeta=d$mqtl_se^2, snp=d$snp, MAF=d$mqtl_maf, N=n_mqtl, type="quant", LD=ld, sdY=1)
D2 <- list(beta=d$gwas_beta, varbeta=d$gwas_se^2, snp=d$snp, MAF=d$gwas_maf, N=n_gwas, type="quant", LD=ld, sdY=1)
s1 <- tryCatch(runsusie(D1), error=function(e) NULL)
s2 <- tryCatch(runsusie(D2), error=function(e) NULL)
if(is.null(s1) || is.null(s2)){ cat(args[3], "susie_failed NA\n"); quit() }
res <- tryCatch(coloc.susie(s1, s2), error=function(e) NULL)
if(is.null(res) || is.null(res$summary)){ cat(args[3], "coloc_failed NA\n"); quit() }
# best PP4 across credible-set pairs
pp4 <- max(res$summary$PP.H4.abf, na.rm=TRUE)
cat(args[3], "PP4", round(pp4,3), "\n")
