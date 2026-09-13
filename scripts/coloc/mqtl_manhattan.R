suppressMessages(library(CMplot))
setwd("/scratch/cy94/rs4477/figures_thesis")
d <- read.table("/scratch/cy94/rs4477/1kgp-cis-mqtl/tensorqtl/genome_452/significant_cpgs.tsv",
                header=TRUE, sep="\t")
# parse chr + pos from phenotype_id (chrN_pos)
d$Chr <- as.integer(sub("chr","",sub("_.*","",d$phenotype_id)))
d$BP  <- as.integer(sub(".*_","",d$phenotype_id))
d <- d[!is.na(d$Chr) & !is.na(d$BP) & !is.na(d$pval_beta),]
cm <- data.frame(SNP=d$phenotype_id, Chr=d$Chr, BP=d$BP, P=d$pval_beta)
cat("significant CpGs:", nrow(cm), "\n")
# cap extreme p (for plotting) 
cm$P[cm$P < 1e-300] <- 1e-300
CMplot(cm, plot.type="m", threshold=0.05, threshold.col="red", threshold.lty=2,
       col=c("#4197d8","#f8c120"), file="jpg", dpi=200, file.output=TRUE,
       main="Genome-wide significant cis-mQTLs", verbose=FALSE,
       ylim=c(0, 50))
cat("mQTL Manhattan done\n")
