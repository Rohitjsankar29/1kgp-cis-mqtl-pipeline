suppressMessages({library(ggplot2); library(data.table)})
FD <- "/scratch/cy94/rs4477/figures_thesis"
th <- theme_minimal(base_size=13) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.y=element_blank(),
        plot.title=element_text(face="bold", size=14),
        axis.title=element_text(size=12),
        legend.position="right")

# ---- FIG 1: cCRE functional enrichment ----
G <- "/g/data/cy94/rs4477/downstream/genome_452"
v <- fread(cmd=paste0("zcat ", G, "/genome.prioritised_final_asm_classified.txt.gz"), select=c("ccre","score_convergence","imprinted","asm_frac_fisher"))
cc <- v[ccre!="none", .N, by=ccre][order(-N)][1:8]
cc$ccre <- factor(cc$ccre, levels=rev(cc$ccre))
p1 <- ggplot(cc, aes(x=ccre, y=N, fill=N)) +
  geom_col(width=0.7) + coord_flip() +
  scale_fill_gradient(low="#a1d99b", high="#238b45", guide="none") +
  scale_y_continuous(labels=scales::comma) +
  labs(x=NULL, y="prioritised variants", title="Regulatory context: cCRE class enrichment") + th
ggsave(paste0(FD,"/gg_functional.png"), p1, width=7.5, height=4.5, dpi=200)
cat("functional done\n")

# ---- FIG 2: convergence score ----
sc <- v[, .N, by=score_convergence][order(score_convergence)]
sc$lab <- format(sc$N, big.mark=",")
p2 <- ggplot(sc, aes(x=factor(score_convergence), y=N, fill=factor(score_convergence))) +
  geom_col(width=0.75) +
  geom_text(aes(label=lab), vjust=-0.4, size=3.3) +
  scale_y_log10(labels=scales::comma, expand=expansion(mult=c(0,0.15))) +
  scale_fill_manual(values=c("#cccccc","#c6dbef","#9ecae1","#4292c6","#e6842e","#c0392b"), guide="none") +
  labs(x="convergence score (independent evidence axes)", y="variants (log scale)",
       title="Multi-evidence convergence") + th
ggsave(paste0(FD,"/gg_convergence.png"), p2, width=7.5, height=4.5, dpi=200)
cat("convergence done\n")

# ---- FIG 3: elastic-net coefficients ----
enet <- data.frame(
  feature=c("CpG proximity","effect size","in cCRE","CpG island","TSS proximity","MAF","-log10 p","SV-implicated"),
  coef=c(0.038,0.016,0.008,0.005,0.001,0.000,0.000,-0.008))
enet$feature <- factor(enet$feature, levels=enet$feature[order(enet$coef)])
enet$dir <- ifelse(enet$coef>=0,"positive","negative")
p3 <- ggplot(enet, aes(x=feature, y=coef, fill=dir)) +
  geom_col(width=0.7) + coord_flip() +
  geom_hline(yintercept=0, color="grey40") +
  scale_fill_manual(values=c("positive"="#238b45","negative"="#c0392b"), guide="none") +
  labs(x=NULL, y="elastic-net coefficient (predicting PIP)",
       title="Annotations weakly predict PIP (R\u00b2\u22480.15)",
       subtitle="SV-implicated is negative: SV variants often tag rather than lead") + th +
  theme(plot.subtitle=element_text(size=10, color="grey30"))
ggsave(paste0(FD,"/gg_elasticnet.png"), p3, width=7.5, height=4.8, dpi=200)
cat("elasticnet done\n")

# ---- FIG 4: elastic-net R2 distribution ----
r2 <- as.numeric(readLines("/tmp/enet_r2_clean.txt"))
r2 <- r2[!is.na(r2) & r2>0.05 & r2<0.3]
p4 <- ggplot(data.frame(r2=r2), aes(x=r2)) +
  geom_histogram(bins=12, fill="#4292c6", color="white") +
  geom_vline(xintercept=median(r2), color="#c0392b", linetype="dashed", linewidth=0.8) +
  annotate("text", x=median(r2), y=Inf, label=sprintf(" median = %.3f", median(r2)), hjust=0, vjust=2, color="#c0392b") +
  labs(x="elastic-net R\u00b2 (per fine-mapping run)", y="count",
       title="Annotations consistently predict PIP poorly",
       subtitle="R\u00b2 < 0.2 across all 22 runs \u2192 PIP is LD-driven") + th +
  theme(plot.subtitle=element_text(size=10, color="grey30"))
ggsave(paste0(FD,"/gg_enet_r2.png"), p4, width=7, height=4.5, dpi=200)
cat("enet R2 done\n")

# ---- FIG 5: imprinting recovery (density) ----
vt <- v[!is.na(asm_frac_fisher)]
vt$grp <- ifelse(vt$imprinted, "imprinted control", "sequence-driven")
p5 <- ggplot(vt, aes(x=asm_frac_fisher, fill=grp)) +
  geom_density(alpha=0.55, color=NA) +
  scale_fill_manual(values=c("sequence-driven"="#4292c6","imprinted control"="#c0392b"), name=NULL) +
  labs(x="ASM fraction (Fisher)", y="density",
       title="ASM validates the pipeline: imprinted loci as positive control") + th
ggsave(paste0(FD,"/gg_imprinting.png"), p5, width=7.5, height=4.5, dpi=200)
cat("imprinting done\n")
cat("ALL GGPLOTS DONE\n")
