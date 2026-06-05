#!/usr/bin/env Rscript
# Cluster all 25,961 Hofbauer cells (fixed metadata) and export DEG tables
# Input: all_hofbauer_final_fixed.h5ad
# Output: Clustered RDS + DEG tables for manual annotation

library(Seurat)
library(Matrix)
library(harmony)
library(anndata)
library(dplyr)

# Paths
INPUT <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/classification/all_hofbauer_final_fixed.h5ad'
OUTPUT_DIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results'

cat("=" , rep("=", 59), "\n", sep="")
cat("Cluster All Hofbauer Cells (Fixed Metadata) + Export DEG Tables\n")
cat("=" , rep("=", 59), "\n", sep="")

# 1. Load data
cat("\n[1/7] Loading all Hofbauer cells...\n")
adata <- read_h5ad(INPUT)
cat("  Total cells:", nrow(adata$obs), "\n")
cat("  Genes:", ncol(adata$X), "\n")

# 2. Convert to Seurat
cat("\n[2/7] Converting to Seurat object...\n")
mat <- t(as.matrix(adata$X))
rownames(mat) <- adata$var_names
colnames(mat) <- adata$obs_names

seu <- CreateSeuratObject(counts = mat, project = 'Hofbauer_Atlas_Final')

# Add metadata
seu$dataset <- as.character(adata$obs$dataset)
seu$disease <- as.character(adata$obs$disease)
seu$disease_group <- as.character(adata$obs$disease_group)

cat("  Seurat object created:", ncol(seu), "cells\n")
cat("  Dataset distribution:\n")
print(table(seu$dataset))

cat("\n  Disease group distribution:\n")
print(table(seu$disease_group))

# 3. Normalize and scale
cat("\n[3/7] Normalizing and scaling...\n")
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 2000)
seu <- ScaleData(seu)

# 4. PCA + Harmony
cat("\n[4/7] Running PCA + Harmony...\n")
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

cat("  Running Harmony batch correction...\n")
seu <- RunHarmony(
  seu,
  group.by.vars = 'dataset',
  reduction.use = 'pca',
  dims.use = 1:30,
  reduction.save = 'harmony',
  project.dim = TRUE,
  verbose = FALSE
)

# 5. UMAP + Clustering
cat("\n[5/7] UMAP + Clustering...\n")
seu <- RunUMAP(seu, reduction = 'harmony', dims = 1:30, verbose = FALSE)
seu <- FindNeighbors(seu, reduction = 'harmony', dims = 1:30, verbose = FALSE)

# Try multiple resolutions
for (res in c(0.5, 0.8, 1.0, 1.2)) {
  seu <- FindClusters(seu, resolution = res, verbose = FALSE)
  cat(sprintf("  Resolution %.1f: %d clusters\n", res, length(unique(seu$seurat_clusters))))
}

# Use resolution 0.8 as default
seu <- FindClusters(seu, resolution = 0.8, verbose = FALSE)

# 6. Save clustered data
cat("\n[6/7] Saving clustered data...\n")
output_file <- file.path(OUTPUT_DIR, 'all_hofbauer_clustered_final.rds')
saveRDS(seu, output_file)
cat("  Saved:", output_file, "\n")
cat("  File size:", round(file.size(output_file) / 1e6, 1), "MB\n")

# 7. Export DEG tables
cat("\n[7/7] Exporting DEG tables...\n")
cat("  This may take 10-15 minutes...\n")

# Join layers first (Seurat v5 requirement)
seu <- JoinLayers(seu)

# Find markers for each cluster
am <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.5, test.use = 't')

# Filter non-coding/ribosomal/mitochondrial
am <- am[!grepl("^(RPL|RPS|MT-|MALAT1|LINC|AC0|AL|AP0|RP11|RP1-|CTD-|XXbac|Z9|THUMPD)", am$gene), ]

# Save full DEG table
write.csv(am, file.path(OUTPUT_DIR, 'all_markers_final_fixed.csv'), row.names = FALSE)
cat("  Saved all_markers_final_fixed.csv\n")

# Output top 15 markers per cluster
cat("\nTop 15 markers per cluster:\n")
cat(rep("=", 80), "\n", sep="")
sink(file.path(OUTPUT_DIR, 'top_markers_per_cluster_final_fixed.txt'))
for (cl in sort(unique(am$cluster))) {
  dd <- am[am$cluster == cl, ]
  dd <- dd[order(-dd$avg_log2FC), ]
  dd <- head(dd, 15)
  
  cat(sprintf("\nC%s (n=%d cells):\n", cl, sum(seu$seurat_clusters == cl)))
  cat(rep("-", 60), "\n", sep="")
  
  # Dataset distribution for this cluster
  ds_dist <- table(seu$dataset[seu$seurat_clusters == cl])
  cat("  Dataset distribution:\n")
  for (ds in names(ds_dist)) {
    cat(sprintf("    %s: %d\n", ds, ds_dist[ds]))
  }
  
  # Disease group distribution
  dg_dist <- table(seu$disease_group[seu$seurat_clusters == cl])
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
cat("  Saved top_markers_per_cluster_final_fixed.txt\n")

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

# Filter to genes present in the data
key_genes_present <- key_genes[key_genes %in% rownames(seu)]
cat("  Key genes present:", length(key_genes_present), "/", length(key_genes), "\n")

# Calculate average expression per cluster
avg_expr <- AverageExpression(seu, features = key_genes_present, assay = 'RNA')
avg_df <- as.data.frame(avg_expr$RNA)
avg_df$Gene <- rownames(avg_df)

# Save
write.csv(avg_df, file.path(OUTPUT_DIR, 'key_gene_expression_final_fixed.csv'), row.names = FALSE)
cat("  Saved key_gene_expression_final_fixed.csv\n")

# Print key gene expression
cat("\nKey gene expression by cluster:\n")
cat(rep("=", 80), "\n", sep="")
print(round(avg_expr$RNA, 2))

# Summary
cat("\n" , rep("=", 60), "\n", sep="")
cat("Summary:\n")
cat("  Total cells:", ncol(seu), "\n")
cat("  Clusters:", length(unique(seu$seurat_clusters)), "\n")
cat("  Genes:", nrow(seu), "\n")

cat("\nDataset distribution:\n")
print(table(seu$dataset))

cat("\nDisease group distribution:\n")
print(table(seu$disease_group))

cat("\nCluster distribution:\n")
print(table(seu$seurat_clusters))

cat("\nOutput files:\n")
cat("  1. all_hofbauer_clustered_final.rds - Clustered Seurat object\n")
cat("  2. top_markers_per_cluster_final_fixed.txt - DEG tables for annotation\n")
cat("  3. key_gene_expression_final_fixed.csv - Key gene expression matrix\n")
cat("  4. all_markers_final_fixed.csv - Full DEG table\n")

cat("\n", rep("=", 60), "\n", sep="")
cat("Done! Please review the DEG tables and annotate clusters.\n")
cat("Next step: Create cluster_assignments_final.csv and apply annotation.\n")
cat(rep("=", 60), "\n", sep="")
