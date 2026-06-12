#!/usr/bin/env Rscript
# Fig6 validation: CEBPA de-repression vs STAT3 activation of SPP1/FN1
library(Seurat); library(ggplot2); library(dplyr); library(patchwork)

seu <- readRDS("results/Hofbauer_Atlas_Final.rds")
expr <- GetAssayData(seu, assay="RNA", layer="data")

# Get expression for key genes
genes <- c("CEBPA","STAT3","NFKB1","SPP1","FN1","PTPRM","CD44")
df <- data.frame(
  CEBPA = expr["CEBPA",],
  STAT3 = expr["STAT3",],
  NFKB1 = expr["NFKB1",],
  SPP1  = expr["SPP1",],
  FN1   = expr["FN1",],
  PTPRM = expr["PTPRM",],
  CD44  = expr["CD44",]
)

df$CEBPA_group <- factor(ntile(df$CEBPA, 3), levels=1:3, labels=c("low","mid","high"))
df$STAT3_group <- factor(ntile(df$STAT3, 3), levels=1:3, labels=c("low","mid","high"))
df$NFKB1_group <- factor(ntile(df$NFKB1, 3), levels=1:3, labels=c("low","mid","high"))

# Plot: SPP1/FN1 by CEBPA group vs STAT3 group
make_box <- function(data, gene, group_var, group_label, color) {
  df_plot <- data.frame(
    expr = data[[gene]],
    group = data[[group_var]]
  ) %>% filter(!is.na(group))
  
  ggplot(df_plot, aes(x=group, y=expr, fill=group)) +
    geom_violin(scale="width", trim=TRUE, linewidth=0.3, alpha=0.8) +
    geom_boxplot(width=0.15, outlier.size=0.1, alpha=0.5, fill="white", linewidth=0.3) +
    scale_fill_manual(values=c("low"="#4575B4","mid"="grey80","high"="#D73027")) +
    labs(title=paste(gene, "by", group_label), y="Expression", x="") +
    theme_classic(base_size=10) +
    theme(legend.position="none", plot.title=element_text(face="bold",size=12,hjust=0.5))
}

p1 <- make_box(df, "SPP1", "CEBPA_group", "CEBPA")
p2 <- make_box(df, "SPP1", "STAT3_group", "STAT3")
p3 <- make_box(df, "SPP1", "NFKB1_group", "NFKB1")
p4 <- make_box(df, "FN1", "CEBPA_group", "CEBPA")
p5 <- make_box(df, "FN1", "STAT3_group", "STAT3")
p6 <- make_box(df, "FN1", "NFKB1_group", "NFKB1")
p7 <- make_box(df, "PTPRM", "NFKB1_group", "NFKB1")
p8 <- make_box(df, "CD44", "NFKB1_group", "NFKB1")

combined <- (p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | plot_spacer()) +
  plot_annotation(title="TF regulation mode: de-repression vs activation",
    theme=theme(plot.title=element_text(face="bold",size=14,hjust=0.5)))

ggsave("figures/Fig6/Fig6_de-repression_validation.png", combined, w=14, h=10, dpi=300, bg="white")

# Stats
cat("\n=== SPP1: CEBPA-low vs CEBPA-high ===\n")
cebpa_lo <- df$SPP1[df$CEBPA_group == "low"]
cebpa_hi <- df$SPP1[df$CEBPA_group == "high"]
cat(sprintf("CEBPA-low mean: %.3f  CEBPA-high mean: %.3f  fold: %.2f  p=%.2e\n",
  mean(cebpa_lo), mean(cebpa_hi), mean(cebpa_lo)/mean(cebpa_hi),
  wilcox.test(cebpa_lo, cebpa_hi, exact=FALSE)$p.value))

stat3_lo <- df$SPP1[df$STAT3_group == "low"]
stat3_hi <- df$SPP1[df$STAT3_group == "high"]
cat(sprintf("STAT3-low mean: %.3f  STAT3-high mean: %.3f  fold: %.2f  p=%.2e\n",
  mean(stat3_lo), mean(stat3_hi), mean(stat3_lo)/mean(stat3_hi),
  wilcox.test(stat3_lo, stat3_hi, exact=FALSE)$p.value))

cat("\n=== FN1: CEBPA-low vs CEBPA-high ===\n")
fn1_ceb_lo <- df$FN1[df$CEBPA_group == "low"]
fn1_ceb_hi <- df$FN1[df$CEBPA_group == "high"]
cat(sprintf("CEBPA-low mean: %.3f  CEBPA-high mean: %.3f  fold: %.2f\n",
  mean(fn1_ceb_lo), mean(fn1_ceb_hi), mean(fn1_ceb_lo)/mean(fn1_ceb_hi)))

cat("\n=== PTPRM: NFKB1-low vs NFKB1-high ===\n")
nfkb_lo <- df$PTPRM[df$NFKB1_group == "low"]
nfkb_hi <- df$PTPRM[df$NFKB1_group == "high"]
cat(sprintf("NFKB1-low mean: %.3f  NFKB1-high mean: %.3f  fold: %.2f  p=%.2e\n",
  mean(nfkb_lo), mean(nfkb_hi), mean(nfkb_hi)/mean(nfkb_lo),
  wilcox.test(nfkb_lo, nfkb_hi, exact=FALSE)$p.value))

cat("done\n")
