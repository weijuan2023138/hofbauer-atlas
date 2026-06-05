#!/usr/bin/env Rscript
# Step 2: Merge existing Atlas with UCSF HB cells + Harmony integration
# Input: existing Atlas RDS + ucsf_hb.h5ad
# Output: unified_with_ucsf.rds (Harmony-corrected, re-clustered)

library(Seurat)
library(Matrix)
library(harmony)
library(anndata)
library(dplyr)

# Paths
EXISTING_ATLAS <- '/home/weijuan/文档/胎盘单细胞数据/shiny_data/meta.rds'
UCSF_HB <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/ucsf_hb.h5ad'
OUTPUT_DIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results'

cat("=" , rep("=", 59), "\n", sep="")
cat("Step 2: Merge + Harmony Integration\n")
cat("=" , rep("=", 59), "\n", sep="")

# 1. Load existing Atlas
cat("\n[1/6] Loading existing Atlas...\n")
existing_meta <- readRDS(EXISTING_ATLAS)
cat("  Existing Atlas cells:", nrow(existing_meta), "\n")
cat("  Columns:", paste(colnames(existing_meta), collapse=", "), "\n")

# Load the full Seurat object (need expression data)
existing_rds <- '/home/weijuan/文档/胎盘单细胞数据/results/unified/pure_macrophages_annotated.rds'
if (file.exists(existing_rds)) {
  existing <- readRDS(existing_rds)
  cat("  Loaded full Seurat object:", ncol(existing), "cells\n")
} else {
  stop("Existing Atlas RDS not found at: ", existing_rds)
}

# 2. Load UCSF HB cells
cat("\n[2/6] Loading UCSF HB cells...\n")
ucsf <- read_h5ad(UCSF_HB)
cat("  UCSF HB cells:", nrow(ucsf$obs), "\n")
cat("  Genes:", ncol(ucsf$X), "\n")

# Convert to Seurat
# ucsf$X is log-normalized, need to convert to counts for Seurat
cat("  Converting to Seurat object...\n")
ucsf_mat <- t(as.matrix(ucsf$X))
rownames(ucsf_mat) <- ucsf$var_names
colnames(ucsf_mat) <- ucsf$obs_names

# Create Seurat object
ucsf_seu <- CreateSeuratObject(counts = ucsf_mat, project = 'UCSF_Li_2026')

# Add metadata
ucsf_seu$dataset <- 'UCSF_Li_2026'
ucsf_seu$disease_group <- 'Normal 1st trimester'
ucsf_seu$gestational_week <- as.character(ucsf$obs$gestational_week)
ucsf_seu$gestational_age_group <- as.character(ucsf$obs$gestational_age_group)
ucsf_seu$sample_id <- as.character(ucsf$obs$sample_id)
ucsf_seu$origin <- as.character(ucsf$obs$origin)
ucsf_seu$fetal_sex <- as.character(ucsf$obs$fetal_sex)
ucsf_seu$major_class <- 'HB'
ucsf_seu$celltype_fullname <- as.character(ucsf$obs$celltype_fullname)

cat("  UCSF Seurat object created:", ncol(ucsf_seu), "cells\n")

# 3. Merge datasets
cat("\n[3/6] Merging datasets...\n")
# Ensure both have the same assay
DefaultAssay(existing) <- 'RNA'
DefaultAssay(ucsf_seu) <- 'RNA'

# Find common genes
common_genes <- intersect(rownames(existing), rownames(ucsf_seu))
cat("  Common genes:", length(common_genes), "\n")

# Subset to common genes
existing_sub <- subset(existing, features = common_genes)
ucsf_sub <- subset(ucsf_seu, features = common_genes)

# Merge
combined <- merge(existing_sub, ucsf_sub)
cat("  Combined cells:", ncol(combined), "\n")
cat("  Dataset distribution:\n")
print(table(combined$dataset))

# 4. Normalize + Scale
cat("\n[4/6] Normalizing and scaling...\n")
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, nfeatures = 2000)
combined <- ScaleData(combined)

# 5. PCA + Harmony
cat("\n[5/6] Running PCA + Harmony...\n")
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

cat("  Running Harmony batch correction (group.by.vars='dataset')...\n")
cat("  This may take 5-10 minutes...\n")

# Fix for Harmony v2 API - use explicit named parameters
combined <- RunHarmony(
  combined,
  group.by.vars = 'dataset',
  reduction.use = 'pca',
  dims.use = 1:30,
  reduction.save = 'harmony',
  project.dim = TRUE,
  verbose = FALSE
)

# 6. UMAP + Clustering
cat("\n[6/6] UMAP + Clustering...\n")
combined <- RunUMAP(combined, reduction = 'harmony', dims = 1:30, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = 'harmony', dims = 1:30, verbose = FALSE)

# Try multiple resolutions
for (res in c(0.5, 0.8, 1.0)) {
  combined <- FindClusters(combined, resolution = res, verbose = FALSE)
  cat(sprintf("  Resolution %.1f: %d clusters\n", res, length(unique(combined$seurat_clusters))))
}

# Use resolution 0.8 as default
combined <- FindClusters(combined, resolution = 0.8, verbose = FALSE)

# 7. Save
cat("\n[7/7] Saving...\n")
output_file <- file.path(OUTPUT_DIR, 'unified_with_ucsf.rds')
saveRDS(combined, output_file)
cat("  Saved:", output_file, "\n")
cat("  File size:", round(file.size(output_file) / 1e6, 1), "MB\n")

# 8. Summary
cat("\n" , rep("=", 60), "\n", sep="")
cat("Summary:\n")
cat("  Total cells:", ncol(combined), "\n")
cat("  Existing Atlas:", sum(combined$dataset != 'UCSF_Li_2026'), "\n")
cat("  UCSF HB:", sum(combined$dataset == 'UCSF_Li_2026'), "\n")
cat("  Clusters:", length(unique(combined$seurat_clusters)), "\n")
cat("  Genes:", nrow(combined), "\n")

cat("\nDataset distribution:\n")
print(table(combined$dataset))

cat("\nCluster distribution:\n")
print(table(combined$seurat_clusters))

cat("\n", rep("=", 60), "\n", sep="")
cat("Done! Next step: 03_validate_mapping.R\n")
cat(rep("=", 60), "\n", sep="")
