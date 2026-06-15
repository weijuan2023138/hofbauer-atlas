#!/usr/bin/env Rscript
library(Seurat); library(ggplot2)

seu <- readRDS("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds")
labels <- read.csv("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/per_cell_disease_labels.csv")
detail <- labels$disease_detail[1:ncol(seu)]

dc <- rep("Other", ncol(seu))
dc[seu$dataset=="GSE290578" & detail=="PE"] <- "Early PE"
dc[detail %in% c("PreE_SF","gHTN","GSE173193","GSE298119")] <- "Late PE"
dc[(seu$dataset=="GSE290578" & detail=="Normal") | (seu$dataset=="GSE298602" & detail=="Control")] <- "Normal Late"

keep <- dc %in% c("Normal Late","Early PE","Late PE")
seu <- subset(seu, cells=colnames(seu)[keep])
seu$group <- factor(dc[keep], levels=c("Normal Late","Early PE","Late PE"))
seu <- JoinLayers(seu)

deg_ep <- read.csv("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/deg_Early_PE_vs_NormalLate.csv", row.names=1)
deg_lp <- read.csv("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/deg_Late_PE_vs_NormalLate.csv", row.names=1)
get_padj <- function(deg,g){if(g%in%rownames(deg)&&!is.na(deg[g,"p_val_adj"]))deg[g,"p_val_adj"]else 1}

tfs <- c("CEBPA","STAT1","NFKB1","JUN","IRF1","IRF8","FOS","STAT3")
fill_cols <- c("Normal Late"="#66C2A5","Early PE"="#FC8D59","Late PE"="#8DA0CB")
expr_mat <- GetAssayData(seu, assay="RNA", layer="data")
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

for(g in tfs) {
  if(!g %in% rownames(expr_mat)) next
  df <- data.frame(expr=expr_mat[g,], Group=seu$group)
  max_y <- max(df$expr) * 1.2
  ys <- max_y * c(1.07, 1.13)
  
  p <- ggplot(df, aes(x=Group, y=expr, fill=Group)) +
    geom_violin(scale="width", trim=TRUE, linewidth=0.3, alpha=0.8) +
    geom_boxplot(width=0.04, outlier.size=0.02, alpha=0.5, fill="white", linewidth=0.2) +
    scale_fill_manual(values=fill_cols)
  
  comps <- list(c(1,2,ys[1],get_padj(deg_ep,g)), c(1,3,ys[2],get_padj(deg_lp,g)))
  for(cp in comps) {
    pv <- cp[4]; lbl <- ifelse(pv<0.001,"***",ifelse(pv<0.01,"**",ifelse(pv<0.05,"*","ns")))
    p <- p + annotate("segment",x=cp[1],xend=cp[2],y=cp[3],yend=cp[3],color="black",linewidth=0.5) +
      annotate("text",x=mean(as.numeric(cp[1:2])),y=cp[3]*1.01,label=lbl,size=3.2,fontface="bold",vjust=0)
  }
  
  p <- p + labs(title=g, x="", y="Expression") + theme_classic(10) +
    theme(axis.text.x=element_text(angle=0,hjust=0.5,size=11,color="black",face="bold"),
      axis.text.y=element_text(size=9,color="black"),
      axis.line=element_line(color="black",linewidth=0.4),
      plot.title=element_text(face="bold",size=14,hjust=0.5),legend.position="none")
  
  prefix <- if(g %in% c("CEBPA","STAT1","NFKB1","JUN")) "Fig5e_TF_PE_" else "补充Fig5D_TF_PE_"
  ggsave(file.path(FIGDIR,paste0(prefix,g,"_3grp.png")), p, w=4, h=4, dpi=300, bg="white")
}
cat("Done\n")
