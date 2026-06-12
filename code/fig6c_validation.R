#!/usr/bin/env Rscript
# Fig6 validation: CEBPA de-repression vs STAT3 activation of SPP1/FN1
library(Seurat); library(ggplot2); library(dplyr); library(patchwork)

seu <- readRDS("results/Hofbauer_Atlas_Final.rds")
expr <- GetAssayData(seu, assay="RNA", layer="data")

# Get expression for key genes
genes <- c("CEBPA","STAT3","NFKB1","SPP1","FN1","PTPRM","CD44","CD47")
df <- data.frame(
  CEBPA = expr["CEBPA",],
  STAT3 = expr["STAT3",],
  NFKB1 = expr["NFKB1",],
  SPP1  = expr["SPP1",],
  FN1   = expr["FN1",],
  PTPRM = expr["PTPRM",],
  CD44  = expr["CD44",],
  CD47  = expr["CD47",]
)

# Binary split: top 30% vs bottom 30% — with gene names in labels
bin_tf <- function(x, name) {
  pct <- ntile(x, 100)
  grp <- rep(NA, length(x))
  grp[pct <= 30] <- paste0(name, "-low")
  grp[pct > 70] <- paste0(name, "-high")
  factor(grp, levels=c(paste0(name, "-low"), paste0(name, "-high")))
}
df$CEBPA_group <- bin_tf(df$CEBPA, "CEBPA")
df$STAT3_group <- bin_tf(df$STAT3, "STAT3")
df$NFKB1_group <- bin_tf(df$NFKB1, "NFKB1")

# Plot: SPP1/FN1 by CEBPA group vs STAT3 group
make_box <- function(data, gene, group_var, group_label, color) {
  df_plot <- data.frame(
    expr = data[[gene]],
    group = data[[group_var]]
  ) %>% filter(!is.na(group)) %>%
 mutate(group=factor(group))  # drop unused levels

  # Wilcoxon p-value
 lo <- df_plot$expr[df_plot$group == levels(df_plot$group)[1]]
 hi <- df_plot$expr[df_plot$group == levels(df_plot$group)[2]]
 pval <- wilcox.test(lo, hi, exact=FALSE)$p.value
 fc <- mean(lo)/mean(hi)
 if (fc < 1) {
   dir <- "↑"
   fc_show <- 1/fc
 } else {
   dir <- "↓"
   fc_show <- fc
 }
 p_str <- ifelse(pval<0.0001, sprintf("%.1e", pval), sprintf("%.4f", pval))
 p_label <- sprintf("%s%.2f  p=%s", dir, fc_show, p_str)
 ymax <- max(df_plot$expr, na.rm=TRUE) * 1.25
  
 p <- ggplot(df_plot, aes(x=group, y=expr, fill=group)) +
   geom_violin(scale="width", trim=FALSE, linewidth=0.3, alpha=0.8) +
   geom_boxplot(width=0.05, outlier.size=0.03, alpha=0.5, fill="white", linewidth=0.3) +
   scale_fill_manual(values=setNames(c("#4575B4","#D73027"), levels(df_plot$group)), na.translate=FALSE) +
   annotate("segment", x=1, xend=2, y=ymax*0.95, yend=ymax*0.95, color="black", linewidth=0.4) +
   annotate("text", x=1.5, y=ymax*0.98, label=p_label, size=3.5, fontface="plain") +
 labs(title=paste(gene, "by", group_label), y="Expression", x="") +
 theme_classic(base_size=10) +
 theme(legend.position="none", plot.title=element_text(face="bold",size=12,hjust=0.5),
       axis.text.x=element_text(size=9,color="black",face="bold"))
 p
}

p1 <- make_box(df, "SPP1", "CEBPA_group", "CEBPA")
p2 <- make_box(df, "SPP1", "STAT3_group", "STAT3")
p3 <- make_box(df, "SPP1", "NFKB1_group", "NFKB1")
p4 <- make_box(df, "FN1", "CEBPA_group", "CEBPA")
p5 <- make_box(df, "FN1", "STAT3_group", "STAT3")
p6 <- make_box(df, "FN1", "NFKB1_group", "NFKB1")
p7 <- make_box(df, "PTPRM", "NFKB1_group", "NFKB1")
p8 <- make_box(df, "CD44", "NFKB1_group", "NFKB1")
p9 <- make_box(df, "CD47", "NFKB1_group", "NFKB1")

combined <- (p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | p9) +
  plot_annotation(title="TF regulation mode: de-repression vs activation",
    theme=theme(plot.title=element_text(face="bold",size=14,hjust=0.5)))

ggsave("figures/Fig6/Fig6c_validation.png", combined, w=11, h=10, dpi=300, bg="white")

# Stats — binary comparison
cebpa_lo <- df$SPP1[df$CEBPA_group == "CEBPA-low" & !is.na(df$CEBPA_group)]
cebpa_hi <- df$SPP1[df$CEBPA_group == "CEBPA-high" & !is.na(df$CEBPA_group)]
cat(sprintf("CEBPA-low mean: %.3f  CEBPA-high mean: %.3f  fold: %.2f  p=%.2e\n",
  mean(cebpa_lo), mean(cebpa_hi), mean(cebpa_lo)/mean(cebpa_hi),
  wilcox.test(cebpa_lo, cebpa_hi, exact=FALSE)$p.value))

stat3_lo <- df$SPP1[df$STAT3_group == "STAT3-low" & !is.na(df$STAT3_group)]
stat3_hi <- df$SPP1[df$STAT3_group == "STAT3-high" & !is.na(df$STAT3_group)]
cat(sprintf("STAT3-low mean: %.3f  STAT3-high mean: %.3f  fold: %.2f  p=%.2e\n",
  mean(stat3_lo), mean(stat3_hi), mean(stat3_lo)/mean(stat3_hi),
  wilcox.test(stat3_lo, stat3_hi, exact=FALSE)$p.value))

cat("\n=== FN1: CEBPA-low vs CEBPA-high ===\n")
fn1_ceb_lo <- df$FN1[df$CEBPA_group == "CEBPA-low" & !is.na(df$CEBPA_group)]
fn1_ceb_hi <- df$FN1[df$CEBPA_group == "CEBPA-high" & !is.na(df$CEBPA_group)]
cat(sprintf("CEBPA-low mean: %.3f  CEBPA-high mean: %.3f  fold: %.2f\n",
  mean(fn1_ceb_lo), mean(fn1_ceb_hi), mean(fn1_ceb_lo)/mean(fn1_ceb_hi)))

cat("
=== PTPRM: NFKB1-low vs NFKB1-high ===
")
nfkb_lo <- df$PTPRM[df$NFKB1_group == "NFKB1-low" & !is.na(df$NFKB1_group)]
nfkb_hi <- df$PTPRM[df$NFKB1_group == "NFKB1-high" & !is.na(df$NFKB1_group)]
cat(sprintf("NFKB1-low mean: %.3f  NFKB1-high mean: %.3f  fold: %.2f  p=%.2e
",
  mean(nfkb_lo), mean(nfkb_hi), mean(nfkb_hi)/mean(nfkb_lo),
  wilcox.test(nfkb_lo, nfkb_hi, exact=FALSE)$p.value))

cat("
=== CD44: NFKB1-low vs NFKB1-high ===
")
cd44_lo <- df$CD44[df$NFKB1_group == "NFKB1-low" & !is.na(df$NFKB1_group)]
cd44_hi <- df$CD44[df$NFKB1_group == "NFKB1-high" & !is.na(df$NFKB1_group)]
cat(sprintf("NFKB1-low mean: %.3f  NFKB1-high mean: %.3f  fold: %.2f  p=%.2e
",
  mean(cd44_lo), mean(cd44_hi), mean(cd44_hi)/mean(cd44_lo),
  wilcox.test(cd44_lo, cd44_hi, exact=FALSE)$p.value))

cat("
=== CD47: NFKB1-low vs NFKB1-high ===
")
cd47_lo <- df$CD47[df$NFKB1_group == "NFKB1-low" & !is.na(df$NFKB1_group)]
cd47_hi <- df$CD47[df$NFKB1_group == "NFKB1-high" & !is.na(df$NFKB1_group)]
cat(sprintf("NFKB1-low mean: %.3f  NFKB1-high mean: %.3f  fold: %.2f  p=%.2e
",
  mean(cd47_lo), mean(cd47_hi), mean(cd47_hi)/mean(cd47_lo),
  wilcox.test(cd47_lo, cd47_hi, exact=FALSE)$p.value))

cat("done\n")
