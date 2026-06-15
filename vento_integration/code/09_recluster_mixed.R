#!/usr/bin/env Rscript
# Re-cluster at higher resolution for uniform dataset mixing (matching old pipeline)
library(Seurat); library(dplyr)

INPUT  <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_harmony_10datasets.rds"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"

seu <- readRDS(INPUT)

# Remove 7 contaminant clusters
hb_clusters <- c("0","1","2","4","5","6","7","8","9","16")
seu <- subset(seu, RNA_snn_res.0.15 %in% hb_clusters)
cat(sprintf("Hofbauer: %d cells\n", ncol(seu)))

# Re-cluster at higher resolutions (matching old pipeline: 0.5, 0.8, 1.0)
cat("\nRe-clustering at higher resolutions...\n")
for (res in c(0.5, 0.8, 1.0)) {
  seu <- FindClusters(seu, resolution=res, verbose=FALSE)
  nc <- length(unique(seu[[paste0("RNA_snn_res.", res)]]))
  cat(sprintf("  res=%.1f: %d clusters\n", res, nc))
}

# Check dataset mixing at each resolution
for (res in c(0.5, 0.8, 1.0)) {
  col_name <- paste0("RNA_snn_res.", res)
  tab <- table(seu[[col_name]], seu$dataset)
  tab_pct <- round(prop.table(tab, 1)*100, 1)
  
  cat(sprintf("\n=== res=%.1f (%d clusters) ===\n", res, nrow(tab)))
  
  # Show max single-dataset per cluster
  dominated <- 0
  for(cl in rownames(tab_pct)) {
    max_pct <- max(tab_pct[cl,])
    max_ds <- colnames(tab_pct)[which.max(tab_pct[cl,])]
    n_cells <- tab[cl, max_ds]
    cat(sprintf("  C%s: %5d cells | %.0f%% %s\n", cl, sum(tab[cl,]), max_pct, max_ds))
    if(max_pct > 70) dominated <- dominated + 1
  }
  cat(sprintf("  Clusters >70%% single-dataset: %d\n", dominated))
}

# Save with best resolution's cluster IDs
best_res <- "0.5"  # will check and adjust
seu$cluster_new <- seu[[paste0("RNA_snn_res.", best_res)]]
Idents(seu) <- "cluster_new"

# Find markers
seu <- JoinLayers(seu)
cat("\nFinding markers at res=0.5...\n")
markers <- FindAllMarkers(seu, only.pos=TRUE, min.pct=0.3, logfc.threshold=0.5,
                          max.cells.per.ident=500, features=VariableFeatures(seu))
markers <- markers[!grepl("^(RPL|RPS|MT-|MALAT1|LINC|AC0|AL|AP0|RP11|RP1-|CTD-)", markers$gene), ]
write.csv(markers, file.path(OUTDIR, "cluster_markers_res0.5.csv"), row.names=FALSE)

# Per-cluster top markers
cat("\n=== Top 5 markers per cluster (res=0.5) ===\n")
for(cl in sort(unique(seu$cluster_new))) {
  top <- head(markers[markers$cluster==cl, "gene"], 5)
  n <- sum(seu$cluster_new==cl)
  cat(sprintf("  C%s  %4d cells: %s\n", cl, n, paste(top, collapse=", ")))
}

# Save
saveRDS(seu, file.path(OUTDIR, "seurat_hofbauer_reclustered.rds"))
cat(sprintf("\nSaved: seurat_hofbauer_reclustered.rds\n"))
