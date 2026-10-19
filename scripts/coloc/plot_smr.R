# SMR locus plot (Zhu et al. 2016 Nat Genet style)
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
pf <- args[1]; outpng <- args[2]; genelab <- args[3]
L <- readLines(pf)
# find section boundaries
sec <- grep("^\\$", L)
secname <- sub("^\\$(\\w+).*","\\1", L[sec])
getsec <- function(nm){
  i <- which(secname==nm); s <- sec[i]+1
  e <- if(i<length(sec)) sec[i+1]-1 else length(L)
  if(s>e) return(NULL)
  fread(text=paste(L[s:e],collapse="\n"), header=FALSE)
}
probe <- getsec("probe")   # V1 id, V2 chr, V3 bp, ... V8 p_eqtl_at_probe, V9 ??
gwas  <- getsec("GWAS")    # snp associations
eqtl  <- getsec("eQTL")
genes <- getsec("Gene")
# probe section: cols = id, chr, bp, id, NA, NA, strand, p_smr?, se? -> actually per-probe: V3=bp, V8=p
probe[, bp:=as.numeric(V3)]; probe[, p:=as.numeric(V8)]
# GWAS: V1 snp, V2 chr, V3 bp, V4 p (or beta/se). Inspect: fread gives numeric
setnames(gwas, paste0("V",1:ncol(gwas)))
gwas[, bp:=as.numeric(V3)]
# find the p-value column in GWAS (usually last numeric that's a p)
gcol <- ncol(gwas); gwas[, p:=as.numeric(get(paste0("V",gcol)))]
setnames(eqtl, paste0("V",1:ncol(eqtl)))
eqtl[, bp:=as.numeric(V3)]; ecol<-ncol(eqtl); eqtl[, p:=as.numeric(get(paste0("V",ecol)))]

png(outpng, width=1400, height=1100, res=150)
layout(matrix(c(1,2,3),3,1), heights=c(3,3,1.5))
par(mar=c(0.5,5,3,2))
# top: GWAS
plot(gwas$bp/1e6, -log10(gwas$p), pch=20, col="grey50", xaxt="n",
     ylab=expression(-log[10](italic(P))~GWAS), xlab="", main=paste0(genelab," locus"), cex=0.6)
abline(h=-log10(5e-8), lty=2, col="red")
par(mar=c(0.5,5,0.5,2))
# middle: mQTL (eQTL section = the methylation QTL for this probe)
plot(eqtl$bp/1e6, -log10(eqtl$p), pch=20, col="#2c6fbb", xaxt="n",
     ylab=expression(-log[10](italic(P))~mQTL), xlab="", cex=0.8)
par(mar=c(4,5,0.5,2))
# bottom: gene models
plot(NA, xlim=range(gwas$bp/1e6), ylim=c(0,1), yaxt="n", xlab="Chr12 position (Mb)", ylab="Genes", bty="n")
if(!is.null(genes)){
  setnames(genes, paste0("V",1:ncol(genes)))
  for(i in 1:nrow(genes)){
    x1<-as.numeric(genes$V3[i])/1e6; x2<-as.numeric(genes$V4[i])/1e6
    segments(x1,0.5,x2,0.5,lwd=3,col="darkgreen")
  }
}
dev.off()
cat("saved", outpng, "\n")
