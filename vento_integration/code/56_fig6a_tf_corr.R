#!/usr/bin/env Rscript
# Fig6a: TF-communication gene correlation — updated with KLF4, IRF1
library(Seurat); library(ComplexHeatmap); library(circlize)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
expr <- GetAssayData(seu, assay="RNA", layer="data")

# TFs matching Fig6c motif enrichment (developmental)
tfs <- c("CEBPA","STAT1","STAT3","NFKB1")
comm_genes <- c("SPP1","FN1","COL1A2","COL1A1","TGFB1","IGF1","PTPRM",
  "CD44","ITGAV","ITGB1","ITGA1","ITGA2","ITGA5","ITGA9","ITGB5",
  "SDC1","SDC2","SDC4","CD47","THBS1","LAMC1","COL4A1","COL4A2",
  "VEGFA","VEGFB","BMP2","BMP4","GDF5","MIF","MDK")

tfs <- intersect(tfs, rownames(expr))
comm_genes <- intersect(comm_genes, rownames(expr))

cor_mat <- matrix(NA, length(tfs), length(comm_genes), dimnames=list(tfs, comm_genes))
p_mat <- matrix(NA, length(tfs), length(comm_genes), dimnames=list(tfs, comm_genes))
for(i in seq_along(tfs)) for(j in seq_along(comm_genes)) {
  ct <- cor.test(expr[tfs[i],], expr[comm_genes[j],], method="spearman")
  cor_mat[i,j] <- ct$estimate; p_mat[i,j] <- ct$p.value
}

p_adj <- matrix(p.adjust(p_mat, method="fdr"), nrow=nrow(p_mat), dimnames=dimnames(p_mat))

max_abs <- max(abs(cor_mat))
ht <- Heatmap(cor_mat, name="Spearman r",
  col=colorRamp2(c(-max_abs, 0, max_abs), c("#4575B4","white","#D73027")),
  cluster_rows=FALSE, cluster_columns=FALSE,
  row_names_gp=gpar(fontsize=11, fontface="bold"),
  column_names_gp=gpar(fontsize=9, fontface="bold.italic"),
  column_names_rot=45,
  heatmap_legend_param=list(title_gp=gpar(fontsize=10,fontface="bold"), labels_gp=gpar(fontsize=9)),
  width=unit(length(comm_genes)*0.35, "inch"), height=unit(length(tfs)*0.35, "inch"))

png(file.path(FIGDIR,"Fig6a_TF_comm_correlation.png"), w=14, h=4, units="in", res=300, bg="white")
draw(ht); dev.off()
cat(sprintf("Fig6a done: %d TFs x %d genes\n", length(tfs), length(comm_genes)))
