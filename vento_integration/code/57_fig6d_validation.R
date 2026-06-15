#!/usr/bin/env Rscript
# Fig6d: TF→target gene violin validation — updated with KLF4, IRF1
library(Seurat); library(ggplot2); library(dplyr); library(patchwork)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
expr <- GetAssayData(seu, assay="RNA", layer="data")

genes <- c("CEBPA","STAT3","STAT1","NFKB1","SPP1","FN1","PTPRM","CD44","CD47")
df <- data.frame(
  CEBPA = expr["CEBPA",], STAT3 = expr["STAT3",],
  STAT1 = expr["STAT1",], NFKB1 = expr["NFKB1",],
  SPP1  = expr["SPP1",],  FN1   = expr["FN1",],
  PTPRM = expr["PTPRM",], CD44  = expr["CD44",], CD47  = expr["CD47",],
  check.names=FALSE)

bin_tf <- function(x, name) {
  pct <- ntile(x, 100)
  grp <- rep(NA, length(x))
  grp[pct <= 30] <- paste0(name, "-low")
  grp[pct > 70] <- paste0(name, "-high")
  factor(grp, levels=c(paste0(name, "-low"), paste0(name, "-high")))
}
df$CEBPA_group <- bin_tf(df$CEBPA, "CEBPA")
df$STAT3_group <- bin_tf(df$STAT3, "STAT3")
df$STAT1_group <- bin_tf(df$STAT1, "STAT1")
df$NFKB1_group <- bin_tf(df$NFKB1, "NFKB1")

make_box <- function(data, gene, group_var, group_label, color_lo="#4575B4", color_hi="#D73027") {
  dfp <- data.frame(expr=data[[gene]], group=data[[group_var]]) %>%
    filter(!is.na(group)) %>% mutate(group=factor(group))
  lo <- dfp$expr[dfp$group==levels(dfp$group)[1]]
  hi <- dfp$expr[dfp$group==levels(dfp$group)[2]]
  pv <- wilcox.test(lo, hi, exact=FALSE)$p.value
  fc <- mean(hi)/mean(lo)
  p_str <- ifelse(pv<0.0001, sprintf("%.1e", pv), sprintf("%.4f", pv))
  p_label <- sprintf("FC=%.2f  p=%s", fc, p_str)
  ymax <- max(dfp$expr, na.rm=TRUE) * 1.3
  
  ggplot(dfp, aes(x=group, y=expr, fill=group)) +
    geom_violin(scale="width", trim=FALSE, linewidth=0.3, alpha=0.8) +
    geom_boxplot(width=0.05, outlier.size=0.03, alpha=0.5, fill="white", linewidth=0.3) +
    scale_fill_manual(values=setNames(c(color_lo, color_hi), levels(dfp$group))) +
    annotate("segment", x=1, xend=2, y=ymax*0.92, yend=ymax*0.92, color="black", linewidth=0.4) +
    annotate("text", x=1.5, y=ymax*0.96, label=p_label, size=3.2, fontface="plain") +
    labs(title=paste(gene, "by", group_label), y="Expression", x="") +
    theme_classic(10) +
    theme(legend.position="none", plot.title=element_text(face="bold",size=12,hjust=0.5),
      axis.text.x=element_text(size=9,color="black",face="bold"))
}

# Row1: CEBPA → SPP1/FN1, STAT3 → SPP1/FN1
p1 <- make_box(df, "SPP1", "CEBPA_group", "CEBPA")
p2 <- make_box(df, "FN1",  "CEBPA_group", "CEBPA")
p3 <- make_box(df, "SPP1", "STAT3_group", "STAT3")
p4 <- make_box(df, "FN1",  "STAT3_group", "STAT3")

# Row2-3: NFKB1 → SPP1/FN1/PTPRM/CD44/CD47
p5 <- make_box(df, "SPP1",  "NFKB1_group", "NFKB1")
p6 <- make_box(df, "FN1",   "NFKB1_group", "NFKB1")
p7 <- make_box(df, "PTPRM", "NFKB1_group", "NFKB1")
p8 <- make_box(df, "CD44",  "NFKB1_group", "NFKB1")
p9 <- make_box(df, "CD47",  "NFKB1_group", "NFKB1")

combined <- (p1 | p2 | p3 | p4) / (p5 | p6 | p7 | p8 | p9) +
  plot_annotation(title="TF → Target Gene Validation",
    theme=theme(plot.title=element_text(face="bold",size=14,hjust=0.5)))

ggsave(file.path(FIGDIR,"Fig6c_validation_updated.png"), combined, w=16, h=8, dpi=300, bg="white")
cat("Fig6d restored to developmental TFs\n")
