suppressMessages({library(variancePartition); library(data.table)})
# top/candidate CpGs to test
top <- fread("/scratch/cy94/rs4477/coloc/smr_candidate_cpgs.txt", header=FALSE)$V1
cat("candidate CpGs:", length(top), "\n")
# by chromosome
top_by_chr <- split(top, sub("_.*","",top))
set.seed(1)
allM <- list(); samples <- NULL
for(ch_name in names(top_by_chr)){
  ch <- sub("chr","",ch_name)
  f <- sprintf("/scratch/cy94/rs4477/1kgp-cis-mqtl/matrix/genome_452_cpgfilt/chr%s.methylation_Mval.cpgfilt.bed.gz", ch)
  if(!file.exists(f)) next
  m <- fread(cmd=paste("zcat", f))
  setnames(m, 1:4, c("chr","start","end","pid"))
  if(is.null(samples)) samples <- colnames(m)[5:ncol(m)]
  want <- top_by_chr[[ch_name]]
  # match by pid (the CpG id, e.g. chr22_10566891 -> pid may be chr22_10566891)
  sel <- m[pid %in% want]
  if(nrow(sel)==0) next
  # subsample if huge (cap per chrom for speed)
  if(nrow(sel)>500){ set.seed(1); sel <- sel[sample(nrow(sel),500)] }
  sub <- as.matrix(sel[, ..samples]); rownames(sub) <- sel$pid
  allM[[ch_name]] <- sub
  cat(ch_name, "matched", nrow(sub), "\n")
}
M <- do.call(rbind, allM)
cat("total top CpGs tested:", nrow(M), "\n")
# covariates
cov <- fread("/scratch/cy94/rs4477/1kgp-cis-mqtl/covariates/genome_452/chr20.covariates.tsv")
covm <- as.data.frame(t(cov[,-1])); colnames(covm) <- cov$ID; covm$sample <- rownames(covm)
covm <- covm[match(samples, covm$sample),]
info <- data.frame(sex=factor(covm$sex),
  platform=factor(ifelse(covm$platform_R9==1,"R9","R10")), row.names=samples)
bc_cols <- grep("basecaller", cov$ID, value=TRUE)
if(length(bc_cols)>0){ bc<-rep("ref",length(samples)); for(b in bc_cols) bc[covm[[b]]==1]<-b; info$basecaller<-factor(bc) }
pc_cols <- grep("PC", cov$ID, value=TRUE)
for(p in pc_cols) info[[p]] <- as.numeric(covm[[p]])
form <- as.formula(paste("~ (1|sex) + (1|platform)",
  if("basecaller"%in%colnames(info)) "+ (1|basecaller)" else "",
  if(length(pc_cols)>0) paste("+",paste(pc_cols,collapse=" + ")) else ""))
vp <- fitExtractVarPartModel(M, form, info)
saveRDS(vp, "/scratch/cy94/rs4477/varpart_top_result.rds")
png("/scratch/cy94/rs4477/figures_thesis/fig_varpart_top.png", width=1400, height=900, res=160)
print(plotVarPart(sortCols(vp)))
dev.off()
cat("=== TOP CpGs median %% variance explained ===\n")
print(round(sort(sapply(vp, median), decreasing=TRUE)*100, 2))
