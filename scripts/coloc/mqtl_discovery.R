suppressMessages({library(CMplot); library(data.table)})
setwd("/scratch/cy94/rs4477/figures_thesis")
cat("reading 25.9M CpGs...\n")
D <- fread("/scratch/cy94/rs4477/1kgp-cis-mqtl/tensorqtl/genome_452/genome_permutation_all.tsv",
           select=c("phenotype_id","pval_beta","qval"))
D <- D[!is.na(pval_beta)]
D[, Chr := as.integer(sub("chr","",sub("_.*","",phenotype_id)))]
D[, BP := as.integer(sub(".*_","",phenotype_id))]
D <- D[!is.na(Chr) & !is.na(BP)]
cat("total CpGs tested:", nrow(D), "\n")
sig_thresh <- max(D[qval < 0.05]$pval_beta, na.rm=TRUE)
cat("FDR<0.05 corresponds to pval_beta <", sig_thresh, "\n")
sig <- D[qval < 0.05]; nonsig <- D[qval >= 0.05]
set.seed(1); nonsig_s <- nonsig[sample(.N, as.integer(.N*0.03))]
plotD <- rbind(sig, nonsig_s)
cat("significant:", nrow(sig), "| plotted:", nrow(plotD), "\n")
cm <- data.frame(SNP=plotD$phenotype_id, Chr=plotD$Chr, BP=plotD$BP, P=plotD$pval_beta)
cm$P[cm$P < 1e-50] <- 1e-50
CMplot(cm, plot.type="m", threshold=sig_thresh, threshold.col="red", threshold.lty=2,
       col=c("#4197d8","#f8c120"), file="jpg", dpi=300, file.output=TRUE,
       main="Genome-wide cis-mQTL discovery", verbose=FALSE, cex=0.15, ylim=c(0,50))
cat("discovery Manhattan done\n")
