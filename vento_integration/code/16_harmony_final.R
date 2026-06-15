#!/usr/bin/env Rscript
# Harmony integration on 10-dataset merged h5ad
library(Seurat); library(harmony); library(anndata); library(Matrix)

OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"
INPUT <- file.path(OUTDIR, "final_10datasets_merged.h5ad")

cat("Loading merged h5ad...\n")
ad <- read_h5ad(INPUT)
cat(sprintf("Loaded: %d cells, %d genes\n", nrow(ad$obs), ncol(ad$X)))

# Convert to Seurat
mat <- t(as.matrix(ad$X))
rownames(mat) <- ad$var_names
colnames(mat) <- make.unique(ad$obs_names, sep="_")
seu <- CreateSeuratObject(counts=mat, project="Hofbauer_Atlas_10datasets",
                          meta.data=ad$obs)
# Ensure dataset column exists
seu$dataset <- ad$obs$dataset
cat(sprintf("Datasets: %s\n", paste(names(table(seu$dataset)), collapse=", ")))

# Normalize + Harmony
seu <- NormalizeData(seu, normalization.method="LogNormalize", scale.factor=10000)
seu <- FindVariableFeatures(seu, nfeatures=2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs=30)

cat("Running Harmony (theta=2, npcs=30)...\n")
seu <- RunHarmony(seu, group.by.vars="dataset", dims.use=1:30, theta=2)
seu <- RunUMAP(seu, reduction="harmony", dims=1:30)
seu <- FindNeighbors(seu, reduction="harmony", dims=1:30)

# Clear any old clustering columns
for(col in grep("snn_res|seurat_clusters", colnames(seu@meta.data), value=T)) {
  seu[[col]] <- NULL
}

for(res in c(0.5, 0.8, 1.0)) {
  seu <- FindClusters(seu, resolution=res, verbose=FALSE)
  n <- length(unique(seu[[paste0("RNA_snn_res.",res)]]))
  cat(sprintf("  res=%.1f: %d clusters\n", res, n))
}

# Dataset mixing (use seurat_clusters from last resolution)
cat("\n=== Dataset mixing ===\n")
tab <- table(seu$seurat_clusters, seu$dataset)
tab_pct <- round(prop.table(tab, 1)*100, 1)
for(cl in rownames(tab_pct)) {
  max_pct <- max(tab_pct[cl,])
  max_ds <- colnames(tab_pct)[which.max(tab_pct[cl,])]
  cat(sprintf("  C%s: %5d cells, max %.0f%% %s\n", cl, sum(tab[cl,]), max_pct, max_ds))
}

# Save
seu$umap_1 <- Embeddings(seu, "umap")[,1]
seu$umap_2 <- Embeddings(seu, "umap")[,2]
saveRDS(seu, file.path(OUTDIR, "seurat_final_10datasets.rds"))
cat(sprintf("\nSaved: %d cells, %d clusters\n", ncol(seu), length(unique(seu$RNA_snn_res.0.5))))
