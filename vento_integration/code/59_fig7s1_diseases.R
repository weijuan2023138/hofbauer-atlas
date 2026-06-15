#!/usr/bin/env Rscript
library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
set.seed(42)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
detail <- labels$disease_detail[1:ncol(seu)]

dc <- rep(NA_character_, ncol(seu))
dc[seu$dataset %in% c("E-MTAB-12421","E-MTAB-6701") | (seu$dataset=="E-MTAB-12795" & detail=="normal")] <- "Normal_Early"
dc[seu$dataset=="GSE214607"] <- "Miscarriage"
dc[detail %in% c("toxoplasmosis","listeriosis","Plasmodium malariae malaria")] <- "Infection"
dc[(seu$dataset=="GSE290578" & detail=="Normal") | (seu$dataset=="GSE298602" & detail=="Control")] <- "Normal_Late"
dc[detail %in% c("PTL","PTNL")] <- "Preterm"

disease_colors <- c("Normal_Early"="#4575B4","Miscarriage"="#D73027","Infection"="#FDB462",
  "Normal_Late"="#66C2A5","Preterm"="#E41A1C")

# Build individual panels for each disease vs its matched normal
ecm_genes <- c("FN1","SPP1","COL1A1","MMP9")
imm_genes <- c("IL1B","TNF","CXCL8","CD44","CD47")

plot_module <- function(seu, dc, group_labels, group_colors, mod_genes, mod_name) {
  keep <- dc %in% group_labels
  sub <- subset(seu, cells=colnames(seu)[keep])
  mod_score <- rowMeans(FetchData(sub, vars=intersect(mod_genes, rownames(sub))))
  df <- data.frame(score=mod_score, group=factor(dc[keep], levels=group_labels))
  pval_df <- data.frame()
  base <- df$score[df$group==group_labels[1]]
  for(d in group_labels[-1]) {
    dis <- df$score[df$group==d]
    pv <- wilcox.test(base, dis, exact=FALSE)$p.value
    pval_df <- rbind(pval_df, data.frame(disease=d, pval=pv))
  }
  pval_df$label <- ifelse(pval_df$pval<0.001,"***",ifelse(pval_df$pval<0.01,"**",ifelse(pval_df$pval<0.05,"*","ns")))
  pval_df$ypos <- max(df$score) * (1.05 + 0.08*(1:nrow(pval_df)))
  
  p <- ggplot(df, aes(x=group, y=score, fill=group)) +
    geom_violin(scale="width", linewidth=0.3, alpha=0.7) +
    geom_boxplot(width=0.04, outlier.size=0.02, alpha=0.5, fill="white", linewidth=0.2) +
    scale_fill_manual(values=group_colors[group_labels]) +
    labs(title=sprintf("%s module", mod_name), y=sprintf("%s score", mod_name), x="") +
    theme_classic(11) + theme(plot.title=element_text(face="bold",size=13,hjust=0.5), legend.position="none",
      axis.text.x=element_text(angle=0,hjust=0.5,size=10,face="bold",color="black"))
  if(nrow(pval_df)>0) {
    p <- p + geom_segment(data=pval_df, aes(x=1, xend=1+1:nrow(pval_df), y=ypos, yend=ypos), inherit.aes=FALSE, linewidth=0.5, color="black") +
      geom_text(data=pval_df, aes(x=1+0.5, y=ypos*1.01, label=label), inherit.aes=FALSE, size=3, fontface="bold")
  }
  p
}

# Early gestation: Normal_Early vs Miscarriage, Normal_Early vs Infection
p1 <- plot_module(seu, dc, c("Normal_Early","Miscarriage"), disease_colors, ecm_genes, "ECM — Miscarriage")
p2 <- plot_module(seu, dc, c("Normal_Early","Miscarriage"), disease_colors, imm_genes, "Immune — Miscarriage")
p3 <- plot_module(seu, dc, c("Normal_Early","Infection"), disease_colors, ecm_genes, "ECM — Infection")
p4 <- plot_module(seu, dc, c("Normal_Early","Infection"), disease_colors, imm_genes, "Immune — Infection")

# Late gestation: Normal_Late vs Preterm
p5 <- plot_module(seu, dc, c("Normal_Late","Preterm"), disease_colors, ecm_genes, "ECM — Preterm")
p6 <- plot_module(seu, dc, c("Normal_Late","Preterm"), disease_colors, imm_genes, "Immune — Preterm")

combined <- (p1 | p2 | p3 | p4) / (p5 | p6 | plot_spacer() | plot_spacer()) +
  plot_annotation(title="Fig7S1: ECM & Immune modules — other diseases",
    theme=theme(plot.title=element_text(face="bold",size=15,hjust=0.5)))
ggsave(file.path(FIGDIR,"Fig7S1_other_diseases.png"), combined, w=14, h=8, dpi=300, bg="white")
cat("Fig7S1 done\n")
