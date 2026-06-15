#!/usr/bin/env Rscript
# 10-dataset Hofbauer Atlas: cluster markers for manual annotation (optimized)
library(Seurat); library(dplyr)

INPUT  <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_harmony_10datasets.rds"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"

seu <- readRDS(INPUT)
cat(sprintf("Loaded: %d cells, %d clusters (res=0.15)\n", ncol(seu), length(unique(seu$RNA_snn_res.0.15))))

# Only use variable features for speed
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

# Find markers with subsampling (max 500 cells per cluster for speed)
cat("\nFinding cluster markers (wilcox, max 500 cells/cluster)...\n")
markers <- FindAllMarkers(seu, only.pos=TRUE, min.pct=0.3, logfc.threshold=0.5,
                          max.cells.per.ident=500, features=VariableFeatures(seu))
markers <- markers[!grepl("^(RPL|RPS|MT-|MALAT1|LINC|AC0|AL|AP0|RP11|RP1-|CTD-)", markers$gene), ]

write.csv(markers, file.path(OUTDIR, "cluster_markers_10datasets.csv"), row.names=FALSE)

# Per-cluster summary
cat("\n============================================================\n")
cat("CLUSTER SUMMARY — Top 5 markers + QC\n")
cat("============================================================\n\n")

mat_check <- c("HLA-DRA", "HLA-DRB1", "CD74", "FCGR3A")
hb_check  <- c("FOLR2", "CD163", "DAB2", "TREM2")

for(cl_char in sort(unique(Idents(seu)))) {
  cl <- as.numeric(as.character(cl_char))
  cells <- WhichCells(seu, idents=cl_char)
  n <- length(cells)
  
  mat_expr <- mean(colMeans(FetchData(seu, vars=intersect(mat_check, rownames(seu)), cells=cells)))
  hb_expr  <- mean(colMeans(FetchData(seu, vars=intersect(hb_check, rownames(seu)), cells=cells)))
  
  top5 <- head(markers[markers$cluster==cl_char, ], 5)
  top_genes <- if(nrow(top5)>0) paste(top5$gene[1:min(5,nrow(top5))], collapse=", ") else "none"
  
  cat(sprintf("C%-2s  %5d cells  HLA-DRA=%.2f  FOLR2=%.2f  %s\n",
      cl_char, n, mat_expr, hb_expr, top_genes))
  
  ds_tab <- sort(table(seu$dataset[cells]), decreasing=TRUE)
  cat(sprintf("      datasets: %s\n", paste(names(ds_tab)[1:4], ds_tab[1:4], sep="=", collapse=", ")))
  
  dis_tab <- sort(table(seu$disease_group[cells]), decreasing=TRUE)
  cat(sprintf("      diseases: %s\n\n", paste(names(dis_tab)[1:4], dis_tab[1:4], sep="=", collapse=", ")))
}

cat(sprintf("Done! %d markers saved to cluster_markers_10datasets.csv\n", nrow(markers)))
