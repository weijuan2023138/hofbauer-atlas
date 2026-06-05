#!/usr/bin/env Rscript
# Final: keep only confirmed Hofbauer clusters (C0, C1, C2, C3, C4, C6)
# Re-number to 0-5, export final marker tables

library(Seurat)
library(dplyr)

INPUT <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/hofbauer_corrected_clustered.rds'
OUTDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results'

cat("============================================================\n")
cat("Final Hofbauer Atlas: Keep C0,C1,C2,C3,C4,C6 only\n")
cat("============================================================\n\n")

seu <- readRDS(INPUT)
cat("Original:", ncol(seu), "cells,", length(unique(seu$seurat_clusters)), "clusters\n")

# Keep only confirmed Hofbauer
keep <- c(0, 1, 2, 3, 4, 6)
seu <- subset(seu, seurat_clusters %in% keep)
cat("Remaining:", ncol(seu), "cells\n\n")

# Renumber to 0-5
old_labels <- sort(unique(seu$seurat_clusters))
cat("Renumbering:\n")
for(i in seq_along(old_labels)) {
  n <- sum(seu$seurat_clusters == old_labels[i])
  cat(sprintf("  C%s -> C%d (%d cells)\n", old_labels[i], i-1, n))
}

# Create new cluster factor directly
new_clusters <- rep(NA, ncol(seu))
for(i in seq_along(old_labels)) {
  new_clusters[seu$seurat_clusters == old_labels[i]] <- i - 1
}
new_clusters <- as.factor(new_clusters)
names(new_clusters) <- colnames(seu)
seu$cluster_final <- new_clusters
Idents(seu) <- "cluster_final"

# Export DEG
cat("\nRunning FindAllMarkers...\n")
seu <- JoinLayers(seu)
am <- FindAllMarkers(seu, only.pos=TRUE, min.pct=0.3, logfc.threshold=0.5, test.use='t')
am <- am[!grepl("^(RPL|RPS|MT-|MALAT1|LINC|AC0|AL|AP0|RP11|RP1-|CTD-|XXbac|Z9|THUMPD)", am$gene), ]
write.csv(am, file.path(OUTDIR, 'Hofbauer_Atlas_Final_all_markers.csv'), row.names=FALSE)
cat("Saved all_markers (", nrow(am), "rows)\n")

# Top 10
cluster_info <- list(
  "0"=list(name="Pro-inflammatory", desc="IL1B, NR4A3, EGR3"),
  "1"=list(name="MHCII+", desc="HLA-DQA1, HLA-DQB1, HLA-DRB5"),
  "2"=list(name="Homeostatic", desc="TREM2, CD68, PAGE4"),
  "3"=list(name="PRKN+ Autophagy", desc="PRKN, C9, SOX5"),
  "4"=list(name="Vascular remodeling", desc="PAPPA, NOTUM, HLA-G"),
  "5"=list(name="MKI67+ Proliferating", desc="MKI67, BUB1B, HMMR")
)

sink(file.path(OUTDIR, 'Hofbauer_Atlas_Final_top10_markers.txt'))
cat("Hofbauer Atlas Final - Top 10 Markers per Cluster\n\n")
for(cl in 0:5) {
  info <- cluster_info[[as.character(cl)]]
  n <- sum(seu$cluster_final == cl)
  cat(sprintf("Cluster %d: %s (%s)\n", cl, info$name, info$desc))
  cat(sprintf("Cells: %d (%.1f%%)\n\n", n, n/ncol(seu)*100))
  
  dd <- am[am$cluster==cl,]; dd <- dd[order(-dd$avg_log2FC),]; dd <- head(dd, 10)
  for(i in 1:nrow(dd)) {
    cat(sprintf("  %-15s  log2FC=%-5.2f  pct=%.0f%%\n", dd$gene[i], dd$avg_log2FC[i], dd$pct.1[i]*100))
  }
  cat("\n")
}
sink()
cat("Saved top10_markers\n")

# Save
saveRDS(seu, file.path(OUTDIR, 'Hofbauer_Atlas_Final.rds'))
cat("Saved Hofbauer_Atlas_Final.rds (", round(file.size(file.path(OUTDIR, 'Hofbauer_Atlas_Final.rds'))/1e6,1), "MB)\n\n")

# Summary
cat("============================================================\n")
cat("FINAL SUMMARY\n")
cat("============================================================\n")
cat("Total Hofbauer cells:", ncol(seu), "\n")
cat("Clusters: 6\n\n")
for(cl in 0:5) {
  info <- cluster_info[[as.character(cl)]]
  n <- sum(seu$cluster_final==cl)
  cat(sprintf("  C%d %-25s: %5d cells (%.1f%%)\n", cl, info$name, n, n/ncol(seu)*100))
}
cat("\nDataset distribution:\n")
print(table(seu$dataset))
cat("\nDisease group distribution:\n")
print(table(seu$disease_group))
cat("\nDone!\n")
