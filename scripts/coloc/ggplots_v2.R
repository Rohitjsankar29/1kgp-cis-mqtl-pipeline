suppressMessages({library(ggplot2); library(data.table)})
FD <- "/scratch/cy94/rs4477/figures_thesis"
G <- "/g/data/cy94/rs4477/downstream/genome_452"

# robust base theme: generous margins, no minor grid, clean
th <- theme_bw(base_size=13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color="grey92", linewidth=0.3),
    panel.border = element_rect(color="grey70", linewidth=0.4),
    plot.title = element_text(face="bold", size=13, margin=margin(b=4)),
    plot.subtitle = element_text(size=9.5, color="grey35", margin=margin(b=8)),
    plot.margin = margin(t=14, r=18, b=10, l=12),
    axis.title = element_text(size=11),
    axis.text = element_text(size=10, color="grey20"),
    legend.position = "top", legend.title = element_blank()
  )

v <- fread(cmd=paste0("zcat ", G, "/genome.prioritised_final_asm_classified.txt.gz"),
           select=c("ccre","score_convergence","imprinted","asm_frac_fisher"))

# ==== 1. FUNCTIONAL ENRICHMENT ====
cc <- v[ccre!="none", .N, by=ccre][order(-N)][1:8]
cc$ccre <- factor(cc$ccre, levels=rev(cc$ccre))
p1 <- ggplot(cc, aes(x=ccre, y=N)) +
  geom_col(fill="#2a9d54", width=0.68) + coord_flip() +
  scale_y_continuous(labels=scales::comma, expand=expansion(mult=c(0,0.08))) +
  labs(x=NULL, y="prioritised variants", title="Regulatory context: cCRE class enrichment") + th
ggsave(paste0(FD,"/gg_functional.png"), p1, width=7, height=4, dpi=200)

# ==== 2. CONVERGENCE (log y, labels above with headroom) ====
sc <- v[, .N, by=score_convergence][order(score_convergence)]
sc$lab <- format(sc$N, big.mark=",")
p2 <- ggplot(sc, aes(x=factor(score_convergence), y=N)) +
  geom_col(aes(fill=factor(score_convergence)), width=0.72) +
  geom_text(aes(label=lab), vjust=-0.5, size=3.1, color="grey25") +
  scale_y_log10(labels=scales::comma, expand=expansion(mult=c(0,0.25))) +
  scale_fill_manual(values=c("#bdbdbd","#c6dbef","#9ecae1","#4292c6","#e6842e","#c0392b")) +
  guides(fill="none") +
  labs(x="convergence score (independent evidence axes)", y="variants (log scale)",
       title="Multi-evidence convergence") + th
ggsave(paste0(FD,"/gg_convergence.png"), p2, width=7, height=4.2, dpi=200)

# ==== 3. ELASTIC-NET COEFFICIENTS (diverging, symmetric x) ====
enet <- data.frame(
  feature=c("CpG proximity","effect size","in cCRE","CpG island","TSS proximity","MAF","-log10 p","SV-implicated"),
  coef=c(0.038,0.016,0.008,0.005,0.001,0.000,0.000,-0.008))
enet$feature <- factor(enet$feature, levels=enet$feature[order(enet$coef)])
enet$dir <- ifelse(enet$coef>=0,"positive","negative")
xlim <- max(abs(enet$coef))*1.15
p3 <- ggplot(enet, aes(x=feature, y=coef, fill=dir)) +
  geom_col(width=0.68) + coord_flip() +
  geom_hline(yintercept=0, color="grey40", linewidth=0.4) +
  scale_y_continuous(limits=c(-xlim, xlim)) +
  scale_fill_manual(values=c("positive"="#2a9d54","negative"="#c0392b")) + guides(fill="none") +
  labs(x=NULL, y="elastic-net coefficient (predicting PIP)",
       title="Annotations weakly predict fine-mapping PIP",
       subtitle="SV-implicated is negative: SV variants often tag rather than lead (R\u00b2\u22480.15)") + th
ggsave(paste0(FD,"/gg_elasticnet.png"), p3, width=7, height=4.2, dpi=200)

# ==== 4. R2 (strip plot, n=22) ====
r2 <- as.numeric(readLines("/tmp/enet_r2_clean.txt")); r2 <- r2[!is.na(r2) & r2>0.05 & r2<0.3]
med <- median(r2)
p4 <- ggplot(data.frame(r2=r2), aes(x=r2, y=0)) +
  geom_vline(xintercept=med, color="#c0392b", linetype="dashed", linewidth=0.6) +
  geom_jitter(height=0.08, width=0, size=2.8, alpha=0.55, color="#2c6fbb") +
  annotate("text", x=med, y=0.16, label=sprintf("median = %.3f", med), color="#c0392b", size=3.6) +
  scale_x_continuous(limits=c(0.10,0.20), breaks=seq(0.10,0.20,0.02), expand=expansion(mult=0.02)) +
  scale_y_continuous(limits=c(-0.14,0.22), breaks=NULL) +
  labs(x=expression("elastic-net "*R^2*" (per fine-mapping run, n=22)"), y=NULL,
       title="Annotations consistently predict PIP poorly",
       subtitle="R-squared < 0.2 across all runs: PIP is LD-driven, not annotation-driven") + th +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank())
ggsave(paste0(FD,"/gg_enet_r2.png"), p4, width=7, height=3.2, dpi=200)

# ==== 5. IMPRINTING (density) ====
vt <- v[!is.na(asm_frac_fisher)]
vt$grp <- ifelse(vt$imprinted, "imprinted control", "sequence-driven")
p5 <- ggplot(vt, aes(x=asm_frac_fisher, fill=grp)) +
  geom_density(alpha=0.5, color=NA) +
  scale_fill_manual(values=c("sequence-driven"="#4292c6","imprinted control"="#c0392b")) +
  scale_x_continuous(expand=expansion(mult=0.01)) + scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
  labs(x="ASM fraction (Fisher)", y="density",
       title="ASM validates the pipeline: imprinted loci as positive control") + th
ggsave(paste0(FD,"/gg_imprinting.png"), p5, width=7, height=4, dpi=200)

cat("ALL 5 REBUILT\n")
