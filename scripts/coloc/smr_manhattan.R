suppressMessages(library(CMplot))
setwd("/scratch/cy94/rs4477/figures_thesis")
dis <- list(eGFR="smr_eGFR_all.smr", CAD="smr_cvd_all.smr", T2D="smr_t2d_all.smr")
base <- NULL; hits_all <- c()
for(nm in names(dis)){
  d <- read.table(paste0("/scratch/cy94/rs4477/coloc/",dis[[nm]]), header=TRUE, sep="\t")
  d$Chr <- as.integer(d$ProbeChr); d$BP <- as.integer(d$Probe_bp)
  d <- d[!is.na(d$Chr) & !is.na(d$p_SMR),]
  tot <- nrow(d)
  h <- d$probeID[d$p_SMR < 0.05/tot & !is.na(d$p_HEIDI) & d$p_HEIDI>0.05]
  hits_all <- c(hits_all, h)
  sub <- data.frame(SNP=d$probeID, Chr=d$Chr, BP=d$BP, P=d$p_SMR)
  colnames(sub)[4] <- nm
  if(is.null(base)) base <- sub else base <- merge(base, sub, by=c("SNP","Chr","BP"), all=TRUE)
}
cat("combined probes:", nrow(base), "| HEIDI-passing hits:", length(unique(hits_all)), "\n")
CMplot(base, plot.type="m", multracks=TRUE, threshold=0.05/nrow(base),
       threshold.col="red", threshold.lty=2,
       highlight=unique(hits_all), highlight.col="red", highlight.cex=1.2,
       col=c("grey60","grey30"), file="jpg", dpi=200, file.output=TRUE, verbose=FALSE)
cat("SMR Manhattan done\n")
