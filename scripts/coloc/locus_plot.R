suppressMessages({library(ggplot2); library(data.table); library(patchwork)})
args <- commandArgs(trailingOnly=TRUE)
pf <- args[1]; outpng <- args[2]; genelab <- args[3]; probe_focus <- args[4]
L <- readLines(pf)
sec <- grep("^\\$", L); secname <- sub("^\\$(\\w+).*","\\1", L[sec])
getsec <- function(nm){ i<-which(secname==nm); if(length(i)==0) return(NULL)
  s<-sec[i]+1; e<-if(i<length(sec)) sec[i+1]-1 else length(L)
  fread(text=paste(L[s:e],collapse="\n"), header=FALSE, fill=TRUE) }
# GWAS: snp chr bp p (parse); use $GWAS section
gwas <- getsec("GWAS"); setnames(gwas, 1:min(4,ncol(gwas)), c("snp","chr","bp","p")[1:min(4,ncol(gwas))])
gwas <- gwas[!is.na(as.numeric(bp))]; gwas[, bp:=as.numeric(bp)/1e6]; gwas[, p:=as.numeric(p)]
# eQTL/mQTL section
eqtl <- getsec("eQTL"); setnames(eqtl, 1:min(4,ncol(eqtl)), c("snp","chr","bp","p")[1:min(4,ncol(eqtl))])
eqtl <- eqtl[!is.na(as.numeric(bp))]; eqtl[, bp:=as.numeric(bp)/1e6]; eqtl[, p:=as.numeric(p)]
# SNP section has the actual mQTL p per SNP for the focal probe - use $SNP
snp <- getsec("SNP")
# Genes
genes <- getsec("Gene"); setnames(genes, 1:5, c("chr","start","end","strand","gene"))
genes[, start:=as.numeric(start)/1e6]; genes[, end:=as.numeric(end)/1e6]
xr <- range(gwas$bp, na.rm=TRUE)

# Panel 1: GWAS
p1 <- ggplot(gwas, aes(bp, -log10(p))) +
  geom_point(color="grey60", size=0.7, alpha=0.6) +
  geom_hline(yintercept=-log10(5e-8), linetype="dashed", color="red") +
  labs(y=expression(-log[10](italic(P))~"GWAS"), x=NULL, title=paste0(genelab, " locus")) +
  xlim(xr) + theme_classic(base_size=12) +
  theme(axis.text.x=element_blank(), plot.title=element_text(face="bold"))
# Panel 2: mQTL
p2 <- ggplot(eqtl, aes(bp, -log10(p))) +
  geom_point(color="#8B2252", size=0.7, alpha=0.6) +
  labs(y=expression(-log[10](italic(P))~"mQTL"), x=NULL) +
  xlim(xr) + theme_classic(base_size=12) + theme(axis.text.x=element_blank())
# Panel 3: genes
gp <- genes[end>=xr[1] & start<=xr[2]]
gp[, y:=rep(1:4, length.out=.N)]  # stagger to avoid overlap
gp[, mid:=(start+end)/2]
gp[, hilite:=grepl(genelab, gene, ignore.case=TRUE)]
p3 <- ggplot(gp) +
  geom_segment(aes(x=start, xend=end, y=y, yend=y, color=hilite), linewidth=2) +
  geom_text(aes(x=mid, y=y, label=gene, color=hilite), vjust=-0.8, size=2.5) +
  scale_color_manual(values=c("FALSE"="grey50","TRUE"="#8B2252"), guide="none") +
  xlim(xr) + ylim(0.5,4.8) +
  labs(x=paste0("Chromosome ", gwas$chr[1], " (Mb)"), y=NULL) +
  theme_classic(base_size=12) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
comb <- p1/p2/p3 + plot_layout(heights=c(3,3,2))
ggsave(outpng, comb, width=9, height=8, dpi=200)
cat("saved", outpng, "\n")
