#!/usr/bin/env Rscript
# Remove C5, C8, C12 contaminants and export final markers
library(Seurat)
library(dplyr)

INPUT <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/hofbauer_corrected_clustered.rds'
OUTPUT_DIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results'

cat("============================================================\n")
cat("Remove C5, C8, C12 + Export Final Markers\n")
cat("============================================================\n\n")

seu <- readRDS(INPUT)
cat("Original cells:", ncol(seu), "\n")
cat("Original clusters:", length(unique(seu$seurat_clusters)), "\n\n")

# Remove contaminant clusters
contaminants <- c(5, 8, 12)
cat("Removing clusters:", paste(contaminants, collapse=", "), "\n")
seu_clean <- subset(seu, !(seurat_clusters %in% contaminants))
cat("Remaining cells:", ncol(seu_clean), "\n")
cat("Removed:", ncol(seu) - ncol(seu_clean), "cells\n\n")

# Re-run DEG
cat("Running FindAllMarkers...\n")
seu_clean <- JoinLayers(seu_clean)
am <- FindAllMarkers(seu_clean, only.pos=TRUE, min.pct=0.3, logfc.threshold=0.5, test.use='t')
am <- am[!grepl("^(RPL|RPS|MT-|MALAT1|LINC|AC0|AL|AP0|RP11|RP1-|CTD-|XXbac|Z9|THUMPD)", am$gene), ]

write.csv(am, file.path(OUTPUT_DIR, 'all_markers_final.csv'), row.names=FALSE)
cat("Saved all_markers_final.csv (", nrow(am), "rows)\n")

# Top 10 per cluster
sink(file.path(OUTPUT_DIR, 'top10_markers_final.txt'))
for (cl in sort(unique(am$cluster))) {
  dd <- am[am$cluster==cl,]; dd <- dd[order(-dd$avg_log2FC),]; dd <- head(dd, 10)
  cat(sprintf("\nC%s (n=%d cells):\n", cl, sum(seu_clean$seurat_clusters==cl)))
  cat(rep("-", 50), "\n", sep="")
  for (i in 1:nrow(dd)) {
    cat(sprintf("  %-15s  log2FC=%-5.2f  pct=%.0f%%\n", dd$gene[i], dd$avg_log2FC[i], dd$pct.1[i]*100))
  }
}
sink()
cat("Saved top10_markers_final.txt\n")

# Save clean Seurat
saveRDS(seu_clean, file.path(OUTPUT_DIR, 'hofbauer_final_clean.rds'))
cat("Saved hofbauer_final_clean.rds (", round(file.size(file.path(OUTPUT_DIR, 'hofbauer_final_clean.rds'))/1e6, 1), "MB)\n\n")

# Summary
cat("============================================================\n")
cat("Final Summary\n")
cat("============================================================\n")
cat("Total Hofbauer cells:", ncol(seu_clean), "\n")
cat("Clusters:", length(unique(seu_clean$seurat_clusters)), "\n\n")
cat("Dataset distribution:\n")
print(table(seu_clean$dataset))
cat("\nDisease group distribution:\n")
print(table(seu_clean$disease_group))
cat("\nCluster distribution:\n")
print(table(seu_clean$seurat_clusters))
cat("\n============================================================\n")
cat("Done!\n")
