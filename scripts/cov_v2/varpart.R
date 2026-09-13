suppressMessages({library(variancePartition); library(data.table)})
# 1. methylation matrix (sample CpGs from chr22)
meth <- fread("zcat /scratch/cy94/rs4477/1kgp-cis-mqtl/matrix/genome_452_cpgfilt/chr22.methylation_Mval.cpgfilt.bed.gz")
# bed format: #chr start end phenotype_id  sample1 sample2 ...
setnames(meth, 1:4, c("chr","start","end","pid"))
samples <- colnames(meth)[5:ncol(meth)]
# sample 5000 random CpGs
set.seed(1); idx <- sample(nrow(meth), min(5000, nrow(meth)))
M <- as.matrix(meth[idx, ..samples]); rownames(M) <- meth$pid[idx]

# 2. covariates (transposed: rows=covariate, cols=sample)
cov <- fread("/scratch/cy94/rs4477/1kgp-cis-mqtl/covariates/genome_452/chr20.covariates.tsv")
covm <- as.data.frame(t(cov[,-1])); colnames(covm) <- cov$ID
covm$sample <- rownames(covm)
# align to methylation sample order
covm <- covm[match(samples, covm$sample),]
# build metadata: sex, platform, basecaller (from dummies), + PCs
info <- data.frame(
  sex = factor(covm$sex),
  platform = factor(ifelse(covm$platform_R9==1,"R9","R10")),
  row.names = samples
)
# basecaller: find basecaller dummy columns
bc_cols <- grep("basecaller", cov$ID, value=TRUE)
cat("basecaller cols:", bc_cols, "\n")
if(length(bc_cols)>0){
  bc <- rep("ref", length(samples))
  for(b in bc_cols){ bc[covm[[b]]==1] <- b }
  info$basecaller <- factor(bc)
}
# add first 5 genotype + 5 methylation PCs as continuous
pc_cols <- grep("PC|pc", cov$ID, value=TRUE)[1:min(10,length(grep("PC|pc",cov$ID,value=TRUE)))]
for(p in pc_cols) info[[p]] <- as.numeric(covm[[p]])
cat("model covariates:", colnames(info), "\n")

# 3. variancePartition model
pc_terms <- paste(pc_cols, collapse=" + ")
form <- as.formula(paste("~ (1|sex) + (1|platform)",
                         if("basecaller" %in% colnames(info)) "+ (1|basecaller)" else "",
                         if(length(pc_cols)>0) paste("+",pc_terms) else ""))
cat("formula:", deparse(form), "\n")
vp <- fitExtractVarPartModel(M, form, info)
saveRDS(vp, "/scratch/cy94/rs4477/varpart_result.rds")
# 4. plot
png("/scratch/cy94/rs4477/figures_thesis/fig_varpart.png", width=1400, height=900, res=160)
print(plotVarPart(sortCols(vp)))
dev.off()
cat("=== median variance explained ===\n")
print(round(sortCols(vp)[, sapply(sortCols(vp), median)*100],2))
cat("saved fig_varpart.png\n")
