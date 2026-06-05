#!/usr/bin/env Rscript
# Harmony batch correction for all Hofbauer cells
# Input: all_hofbauer_combined.h5ad
# Output: Harmony-corrected UMAP + clustering

library(Seurat)
library(Matrix)
library(harmony)
library(anndata)
library(dplyr)

# Paths
INPUT <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/classification/all_hofbauer_combined.h5ad'
OUTPUT_DIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results'
FIGURES_DIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'

cat("=" , rep("=", 59), "\n", sep="")
cat("Harmony Batch Correction for All Hofbauer Cells\n")
cat("=" , rep("=", 59), "\n", sep="")

# 1. Load data
cat("\n[1/6] Loading combined Hofbauer cells...\n")
adata <- read_h5ad(INPUT)
cat("  Total cells:", nrow(adata$obs), "\n")
cat("  Genes:", ncol(adata$X), "\n")

# 2. Convert to Seurat
cat("\n[2/6] Converting to Seurat object...\n")
mat <- t(as.matrix(adata$X))
rownames(mat) <- adata$var_names
colnames(mat) <- adata$obs_names

seu <- CreateSeuratObject(counts = mat, project = 'Hofbauer_Atlas')

# Add metadata
seu$dataset <- as.character(adata$obs$dataset)
seu$disease <- as.character(adata$obs$disease)
seu$disease_group <- as.character(adata$obs$disease_group)

cat("  Seurat object created:", ncol(seu), "cells\n")
cat("  Dataset distribution:\n")
print(table(seu$dataset))

# 3. Normalize and scale
cat("\n[3/6] Normalizing and scaling...\n")
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 2000)
seu <- ScaleData(seu)

# 4. PCA
cat("\n[4/6] Running PCA...\n")
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

# 5. Harmony batch correction
cat("\n[5/6] Running Harmony batch correction...\n")
cat("  group.by.vars = 'dataset'\n")
cat("  This may take 5-10 minutes...\n")

seu <- RunHarmony(
  seu,
  group.by.vars = 'dataset',
  reduction.use = 'pca',
  dims.use = 1:30,
  reduction.save = 'harmony',
  project.dim = TRUE,
  verbose = FALSE
)

# 6. UMAP + Clustering
cat("\n[6/6] UMAP + Clustering...\n")
seu <- RunUMAP(seu, reduction = 'harmony', dims = 1:30, verbose = FALSE)
seu <- FindNeighbors(seu, reduction = 'harmony', dims = 1:30, verbose = FALSE)

# Try multiple resolutions
for (res in c(0.5, 0.8, 1.0)) {
  seu <- FindClusters(seu, resolution = res, verbose = FALSE)
  cat(sprintf("  Resolution %.1f: %d clusters\n", res, length(unique(seu$seurat_clusters))))
}

# Use resolution 0.8 as default
seu <- FindClusters(seu, resolution = 0.8, verbose = FALSE)

# 7. Save
cat("\n[7/7] Saving...\n")
output_file <- file.path(OUTPUT_DIR, 'all_hofbauer_harmony.rds')
saveRDS(seu, output_file)
cat("  Saved:", output_file, "\n")
cat("  File size:", round(file.size(output_file) / 1e6, 1), "MB\n")

# 8. Summary
cat("\n" , rep("=", 60), "\n", sep="")
cat("Summary:\n")
cat("  Total cells:", ncol(seu), "\n")
cat("  Clusters:", length(unique(seu$seurat_clusters)), "\n")
cat("  Genes:", nrow(seu), "\n")

cat("\nDataset distribution:\n")
print(table(seu$dataset))

cat("\nCluster distribution:\n")
print(table(seu$seurat_clusters))

cat("\n", rep("=", 60), "\n", sep="")
cat("Done! Next step: Assess Harmony-corrected batch effects\n")
cat(rep("=", 60), "\n", sep="")
