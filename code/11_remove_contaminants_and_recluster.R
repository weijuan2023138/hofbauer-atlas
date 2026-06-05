#!/usr/bin/env Rscript
# Remove contaminants and re-cluster with lower resolution (6-8 clusters)
# Input: all_hofbauer_clustered_final.rds
# Output: Clean Hofbauer cells with 6-8 clusters + DEG tables

library(Seurat)
library(Matrix)
library(harmony)
library(dplyr)

# Paths
INPUT <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/all_hofbauer_clustered_final.rds'
OUTPUT_DIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results'

cat("=" , rep("=", 59), "\n", sep="")
cat("Remove Contaminants + Re-cluster (6-8 clusters)\n")
cat("=" , rep("=", 59), "\n", sep="")

# 1. Load clustered data
cat("\n[1/8] Loading clustered data...\n")
seu <- readRDS(INPUT)
cat("  Total cells:", ncol(seu), "\n")
cat("  Original clusters:", length(unique(seu$seurat_clusters)), "\n")

# 2. Define contaminant clusters
cat("\n[2/8] Removing contaminants...\n")

# Confirmed contaminants
confirmed_contaminants <- c(15, 18, 19, 20)

# Highly suspected contaminants
suspected_contaminants <- c(4, 6, 11, 14, 16, 21, 22, 23, 25, 26, 27)

# All contaminants to remove
all_contaminants <- c(confirmed_contaminants, suspected_contaminants)

# Keep true Hofbauer clusters
true_hofbauer_clusters <- c(0, 1, 2, 3, 5, 7, 8, 9, 10, 12, 13, 17, 24)

cat("  Confirmed contaminants:", paste(confirmed_contaminants, collapse=", "), "\n")
cat("  Suspected contaminants:", paste(suspected_contaminants, collapse=", "), "\n")
cat("  True Hofbauer clusters:", paste(true_hofbauer_clusters, collapse=", "), "\n")

# 3. Subset to true Hofbauer cells
cat("\n[3/8] Subsetting to true Hofbauer cells...\n")
seu_clean <- subset(seu, seurat_clusters %in% true_hofbauer_clusters)
cat("  Clean Hofbauer cells:", ncol(seu_clean), "\n")
cat("  Removed:", ncol(seu) - ncol(seu_clean), "contaminant cells\n")

# 4. Re-run Harmony on clean cells
cat("\n[4/8] Re-running Harmony on clean cells...\n")
seu_clean <- NormalizeData(seu_clean)
seu_clean <- FindVariableFeatures(seu_clean, nfeatures = 2000)
seu_clean <- ScaleData(seu_clean)
seu_clean <- RunPCA(seu_clean, npcs = 50, verbose = FALSE)

seu_clean <- RunHarmony(
  seu_clean,
  group.by.vars = 'dataset',
  reduction.use = 'pca',
  dims.use = 1:30,
  reduction.save = 'harmony',
  project.dim = TRUE,
  verbose = FALSE
)

# 5. UMAP + Clustering with lower resolution
cat("\n[5/8] UMAP + Clustering with lower resolution...\n")
seu_clean <- RunUMAP(seu_clean, reduction = 'harmony', dims = 1:30, verbose = FALSE)
seu_clean <- FindNeighbors(seu_clean, reduction = 'harmony', dims = 1:30, verbose = FALSE)

# Try multiple low resolutions
for (res in c(0.2, 0.3, 0.4, 0.5)) {
  seu_clean <- FindClusters(seu_clean, resolution = res, verbose = FALSE)
  cat(sprintf("  Resolution %.1f: %d clusters\n", res, length(unique(seu_clean$seurat_clusters))))
}

# Use resolution 0.3 as default (should give 6-8 clusters)
seu_clean <- FindClusters(seu_clean, resolution = 0.3, verbose = FALSE)

# 6. Save clean clustered data
cat("\n[6/8] Saving clean clustered data...\n")
output_file <- file.path(OUTPUT_DIR, 'hofbauer_clean_clustered.rds')
saveRDS(seu_clean, output_file)
cat("  Saved:", output_file, "\n")
cat("  File size:", round(file.size(output_file) / 1e6, 1), "MB\n")

# 7. Export DEG tables
cat("\n[7/8] Exporting DEG tables...\n")
cat("  This may take 5-10 minutes...\n")

# Join layers first
seu_clean <- JoinLayers(seu_clean)

# Find markers for each cluster
am <- FindAllMarkers(seu_clean, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.5, test.use = 't')

# Filter non-coding/ribosomal/mitochondrial
am <- am[!grepl("^(RPL|RPS|MT-|MALAT1|LINC|AC0|AL|AP0|RP11|RP1-|CTD-|XXbac|Z9|THUMPD)", am$gene), ]

# Save full DEG table
write.csv(am, file.path(OUTPUT_DIR, 'all_markers_clean.csv'), row.names = FALSE)
cat("  Saved all_markers_clean.csv\n")

# Output top 15 markers per cluster
cat("\nTop 15 markers per cluster:\n")
cat(rep("=", 80), "\n", sep="")
sink(file.path(OUTPUT_DIR, 'top_markers_per_cluster_clean.txt'))
for (cl in sort(unique(am$cluster))) {
  dd <- am[am$cluster == cl, ]
  dd <- dd[order(-dd$avg_log2FC), ]
  dd <- head(dd, 15)
  
  cat(sprintf("\nC%s (n=%d cells):\n", cl, sum(seu_clean$seurat_clusters == cl)))
  cat(rep("-", 60), "\n", sep="")
  
  # Dataset distribution
  ds_dist <- table(seu_clean$dataset[seu_clean$seurat_clusters == cl])
  cat("  Dataset distribution:\n")
  for (ds in names(ds_dist)) {
    cat(sprintf("    %s: %d\n", ds, ds_dist[ds]))
  }
  
  # Disease group distribution
  dg_dist <- table(seu_clean$disease_group[seu_clean$seurat_clusters == cl])
  cat("\n  Disease group distribution:\n")
  for (dg in names(dg_dist)) {
    cat(sprintf("    %s: %d\n", dg, dg_dist[dg]))
  }
  
  cat("\n  Top markers:\n")
  for (i in 1:nrow(dd)) {
    cat(sprintf("  %-15s  log2FC=%-6.2f  pct=%.0f%%  padj=%.1e\n",
        dd$gene[i], dd$avg_log2FC[i], dd$pct.1[i] * 100, dd$p_val_adj[i]))
  }
}
sink()
cat("  Saved top_markers_per_cluster_clean.txt\n")

# Generate key gene expression matrix
cat("\n  Generating key gene expression matrix...\n")
key_genes <- c('FOLR2', 'CD163', 'MRC1', 'TNF', 'IL1B', 'CXCL8', 
               'CCL13', 'AIF1', 'CLIC1',
               'C1QA', 'C1QB', 'FCGR3A', 'HLA-DRA',
               'PRKN', 'C9', 'SOX5',
               'MKI67', 'BUB1B', 'KIF4A',
               'SPP1', 'DAB2', 'S100A4',
               'NFKBIZ', 'PTGS2', 'DUSP2',
               'CD36', 'LYVE1', 'TREM2', 'MAF')

key_genes_present <- key_genes[key_genes %in% rownames(seu_clean)]
cat("  Key genes present:", length(key_genes_present), "/", length(key_genes), "\n")

avg_expr <- AverageExpression(seu_clean, features = key_genes_present, assay = 'RNA')
avg_df <- as.data.frame(avg_expr$RNA)
avg_df$Gene <- rownames(avg_df)

write.csv(avg_df, file.path(OUTPUT_DIR, 'key_gene_expression_clean.csv'), row.names = FALSE)
cat("  Saved key_gene_expression_clean.csv\n")

# Print key gene expression
cat("\nKey gene expression by cluster:\n")
cat(rep("=", 80), "\n", sep="")
print(round(avg_expr$RNA, 2))

# 8. Summary
cat("\n[8/8] Summary\n")
cat(rep("=", 60), "\n", sep="")
cat("  Total clean Hofbauer cells:", ncol(seu_clean), "\n")
cat("  Clusters:", length(unique(seu_clean$seurat_clusters)), "\n")
cat("  Genes:", nrow(seu_clean), "\n")

cat("\nDataset distribution:\n")
print(table(seu_clean$dataset))

cat("\nDisease group distribution:\n")
print(table(seu_clean$disease_group))

cat("\nCluster distribution:\n")
print(table(seu_clean$seurat_clusters))

cat("\nOutput files:\n")
cat("  1. hofbauer_clean_clustered.rds - Clean Hofbauer Seurat object\n")
cat("  2. top_markers_per_cluster_clean.txt - DEG tables for annotation\n")
cat("  3. key_gene_expression_clean.csv - Key gene expression matrix\n")
cat("  4. all_markers_clean.csv - Full DEG table\n")

cat("\n", rep("=", 60), "\n", sep="")
cat("Done! Please review the DEG tables and annotate clusters.\n")
cat("Expected: 6-8 clusters matching existing 7 subtypes.\n")
cat(rep("=", 60), "\n", sep="")
