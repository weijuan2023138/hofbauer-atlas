#!/usr/bin/env Rscript
# Fig5e TF violins — corrected groups, PE split into Early/Late
library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
set.seed(42)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
detail <- labels$disease_detail[1:ncol(seu)]
tri_map <- c("E-MTAB-12421"="Early","E-MTAB-6701"="Early","E-MTAB-12795"="Early","GSE214607"="Early","UCSF Li 2026"="Mid","GSE290578"="Late","GSE333257"="Late","GSE298602"="Late","GSE298119"="Late","GSE173193"="Late")
tri <- tri_map[as.character(seu$dataset)]

# Group 1: Normal_Early
# Group 2: Miscarriage = GSE214607
# Group 3: Infection = Hoo infected only (toxoplasmosis, listeriosis, malaria)  
# Group 4: Normal_Late = GSE290578 Norm + GSE298602 Control + UCSF Late (exclude TL)
# Group 5: Early_PE = GSE290578 PE
# Group 6: Late_PE = PreE_SF + gHTN + GSE173193 + GSE298119

dc <- rep("Other", ncol(seu))
dc[seu$dataset %in% c("E-MTAB-12421","E-MTAB-6701") | (seu$dataset=="E-MTAB-12795" & detail=="normal")] <- "Normal_Early"
dc[seu$dataset=="GSE214607"] <- "Miscarriage"
dc[detail %in% c("toxoplasmosis","listeriosis","Plasmodium malariae malaria")] <- "Infection"
dc[(seu$dataset=="GSE290578" & detail=="Normal") | (seu$dataset=="GSE298602" & detail=="Control")] <- "Normal_Late"
dc[seu$dataset=="GSE290578" & detail=="PE"] <- "Early_PE"
dc[detail %in% c("PreE_SF","gHTN","GSE173193","GSE298119")] <- "Late_PE"
dc[detail=="TL"] <- "Excluded"

target <- c("Normal_Early","Miscarriage","Infection","Normal_Late","Early_PE","Late_PE")
cells <- which(dc %in% target)
seu <- subset(seu, cells=colnames(seu)[cells])
dc <- factor(dc[cells], levels=target)

for(g in target) cat(sprintf("%s: %d  ", g, sum(dc==g))); cat("\n")

# DEG for significance: each disease vs its matched normal
seu$grp <- dc
names(dc) <- colnames(seu)

run_deg <- function(disease_grp, control_grp) {
  sub <- subset(seu, grp %in% c(disease_grp, control_grp))
  sub$comparison <- ifelse(sub$grp==disease_grp, "Disease", "Control")
  Idents(sub) <- "comparison"; sub <- JoinLayers(sub)
  FindMarkers(sub, ident.1="Disease", ident.2="Control", logfc.threshold=0, min.pct=0.1, test.use="wilcox")
}

deg_mis <- run_deg("Miscarriage","Normal_Early")
deg_inf <- run_deg("Infection","Normal_Early")
deg_ep  <- run_deg("Early_PE","Normal_Late")
deg_lp  <- run_deg("Late_PE","Normal_Late")

get_padj <- function(deg, gene) {
  if(gene %in% rownames(deg) && !is.na(deg[gene,"p_val_adj"])) deg[gene,"p_val_adj"] else 1
}

# TF violins
all_tfs <- c("STAT1","STAT3","NFKB1","RELB","CEBPA","MAFB","ID2","KLF4","JUN","FOS","IRF1","IRF8")
fill_cols <- c("Normal_Early"="#4575B4","Miscarriage"="#D73027","Infection"="#FDB462",
               "Normal_Late"="#66C2A5","Early_PE"="#FC8D59","Late_PE"="#E41A1C")

seu <- JoinLayers(seu)
expr_mat <- GetAssayData(seu, assay="RNA", layer="data")
meta <- data.frame(group=dc, row.names=colnames(seu))

for(g in all_tfs) {
  if(!g %in% rownames(expr_mat)) next
  df <- data.frame(expr=expr_mat[g,], Group=factor(dc, levels=target))
  max_y <- max(df$expr) * 1.2; ys <- max_y * c(1.08, 1.14)
  
  p <- ggplot(df, aes(x=Group, y=expr, fill=Group)) +
    geom_violin(scale="width", trim=TRUE, linewidth=0.3, alpha=0.8) +
    geom_boxplot(width=0.04, outlier.size=0.02, alpha=0.5, fill="white", linewidth=0.2, fatten=1.5) +
    scale_fill_manual(values=fill_cols)
  
  # Significance: Early PE vs NL, Late PE vs NL, Miscarriage vs NE, Infection vs NE
  comps <- list(
    list(1,2,1.5,ys[1],get_padj(deg_mis,g)),
    list(1,3,2.5,ys[2],get_padj(deg_inf,g)),
    list(4,5,4.5,ys[1],get_padj(deg_ep,g)),
    list(4,6,5.5,ys[2],get_padj(deg_lp,g)))
  comps <- list(list(1,2,ys[1],get_padj(deg_mis,g)), list(1,3,ys[2],get_padj(deg_inf,g)),
                list(4,5,ys[1],get_padj(deg_pe,g)),  list(4,6,ys[2],get_padj(deg_pt,g)))
  for(cp in comps) {
    pv <- cp[[4]]; lbl <- ifelse(pv<0.001,"***",ifelse(pv<0.01,"**",ifelse(pv<0.05,"*","ns")))
    mid_x <- mean(c(cp[[1]], cp[[2]]))
    p <- p + annotate("segment",x=cp[[1]],xend=cp[[2]],y=cp[[3]],yend=cp[[3]],color="black",linewidth=0.5) +
      annotate("text",x=mid_x,y=cp[[3]]*1.005,label=lbl,size=2.5,fontface="bold",vjust=0)
  }
  
  p <- p + labs(title=g, x="", y="Expression") + theme_classic(10) +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=10,color="black",face="bold"),
      axis.text.y=element_text(size=8,color="black"),
      axis.line=element_line(color="black",linewidth=0.4),
      plot.title=element_text(face="bold",size=13,hjust=0.5,color="black"),
      legend.position="none")
  
  if(g %in% c("STAT1","NFKB1","CEBPA","JUN")) {
    ggsave(file.path(FIGDIR,sprintf("Fig5e_TF_%s_corrected.png",g)), p, w=5, h=4, dpi=300, bg="white")
  } else {
    ggsave(file.path(FIGDIR,sprintf("补充Fig5D_TF_%s_corrected.png",g)), p, w=5, h=4, dpi=300, bg="white")
  }
}
cat("All 12 TFs done\n")
