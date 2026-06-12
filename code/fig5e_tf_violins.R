#!/usr/bin/env Rscript
# Fig5e: TF violin plots — 12 TFs across 6 disease groups (synced with Fig7)
setwd("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration")
library(Seurat); library(ggplot2); library(dplyr)

seu <- readRDS("results/Hofbauer_Atlas_Final.rds")
mask_gse <- seu$dataset == "GSE290578"; bcs <- colnames(seu)
is_norm_gse <- grepl("_Norm_", bcs) & mask_gse; is_pt_gse <- grepl("_Pt_", bcs) & mask_gse
seu$disease_clean <- NA
seu$disease_clean[seu$disease == "Normal_1st"] <- "Normal_Early"
seu$disease_clean[seu$disease == "Normal"] <- "Normal_Late"
seu$disease_clean[seu$disease == "RM/NC"] <- "Miscarriage"
seu$disease_clean[seu$disease == "Normal/Listeria/Toxoplasma/Malaria"] <- "Infection"
seu$disease_clean[is_norm_gse] <- "Normal_Late"; seu$disease_clean[is_pt_gse] <- "PE"
seu$disease_clean[seu$disease %in% c("PTL","PTNL")] <- "Preterm"
seu$disease_clean[seu$disease == "TL"] <- "Normal_Late"
keep <- seu$disease_clean %in% c("Normal_Early","Normal_Late","Miscarriage","PE","Infection","Preterm")
seu_sub <- subset(seu, cells=colnames(seu)[keep])

groups <- c("Normal_Early","Miscarriage","Infection","Normal_Late","PE","Preterm")
fill_cols <- c("#4575B4","#D73027","#FDB462","#66C2A5","#FC8D59","#E41A1C")

all_tfs <- c("STAT1","STAT3","NFKB1","RELB","CEBPA","MAFB","ID2","KLF4","JUN","FOS","IRF1","IRF8")

mis <- read.csv("results/deg_Miscarriage_vs_NormalEarly.csv", row.names=1)
inf <- read.csv("results/deg_Infection_vs_NormalEarly.csv", row.names=1)
pe  <- read.csv("results/deg_PE_vs_NormalLate.csv", row.names=1)
pt  <- read.csv("results/deg_Preterm_vs_NormalLate.csv", row.names=1)
get_padj <- function(deg, gene) {
  if(gene %in% rownames(deg) && !is.na(deg[gene,"p_val_adj"])) deg[gene,"p_val_adj"] else 1
}

expr <- GetAssayData(seu_sub, assay="RNA", layer="data"); meta <- seu_sub@meta.data

for(g in all_tfs) {
  df <- data.frame(expr=expr[g,], Group=factor(meta$disease_clean, levels=groups))
  max_y <- max(df$expr) * 1.15
  ys <- max_y * c(1.06, 1.12)  # staggered for early and late groups
  
  p <- ggplot(df, aes(x=Group, y=expr, fill=Group)) +
    geom_violin(scale="width", trim=TRUE, linewidth=0.3, alpha=0.8) +
    geom_boxplot(width=0.04, outlier.size=0.02, alpha=0.5, fill="white", linewidth=0.2, fatten=1.5) +
    scale_fill_manual(values=fill_cols)
  
  # Significance lines: horizontal only, no vertical bars
  comps <- list(
    list(1,2,1.5,ys[1],get_padj(mis,g)),  
    list(1,3,2.5,ys[2],get_padj(inf,g)),  
    list(4,5,4.5,ys[1],get_padj(pe,g)),   
    list(4,6,5,ys[2],get_padj(pt,g)))    
  for(cp in comps) {
    pv <- cp[[5]]
    if(pv < 0.05) {
      lbl <- ifelse(pv<0.001,"***",ifelse(pv<0.01,"**","*"))
    } else {
      lbl <- "ns"
    }
    p <- p + 
      annotate("segment", x=cp[[1]], xend=cp[[2]], y=cp[[4]], yend=cp[[4]],
               color="black", linewidth=0.5) +
      annotate("text", x=cp[[3]], y=cp[[4]]*1.005, label=lbl, 
               size=2.8, fontface="bold", vjust=0)
  }
  
  p <- p + labs(title=g, x="", y="Expression") + theme_classic(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=10, color="black", face="bold"),
          axis.text.y=element_text(size=8, color="black"),
          axis.line=element_line(color="black", linewidth=0.4),
          plot.title=element_text(face="bold", size=13, hjust=0.5, color="black"),
          legend.position="none")
  
  if(g %in% c("STAT1","NFKB1","CEBPA","JUN")) {
    ggsave(sprintf("figures/Fig5/Fig5e_TF_%s.png", g), p, w=5, h=4, dpi=300, bg="white")
  } else {
    ggsave(sprintf("figures/Fig5/补充图Fig5D_TF_%s.png", g), p, w=5, h=4, dpi=300, bg="white")
  }
}
cat("Fig5e: 12 TFs synced\n")
