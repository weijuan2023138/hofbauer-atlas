#!/usr/bin/env Rscript
# Harmony with npcs=20 (fewer PCs = fewer clusters)
library(Seurat); library(harmony); library(anndata); library(Matrix)

OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"
INPUT <- file.path(OUTDIR, "final_10datasets_merged.h5ad")

cat("Loading...\n")
ad <- read_h5ad(INPUT)
mat <- t(as.matrix(ad$X))
rownames(mat) <- ad$var_names
colnames(mat) <- make.unique(ad$obs_names, sep="_")
seu <- CreateSeuratObject(counts=mat, meta.data=ad$obs)
seu$dataset <- ad$obs$dataset

seu <- NormalizeData(seu, normalization.method="LogNormalize", scale.factor=10000)
seu <- FindVariableFeatures(seu, nfeatures=2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs=20)

cat("Running Harmony (theta=2, npcs=20)...\n")
seu <- RunHarmony(seu, group.by.vars="dataset", dims.use=1:20, theta=2)
seu <- RunUMAP(seu, reduction="harmony", dims=1:20)
seu <- FindNeighbors(seu, reduction="harmony", dims=1:20)

for(col in grep("snn_res|seurat_clusters", colnames(seu@meta.data), value=T)) seu[[col]] <- NULL

for(res in c(0.5, 0.8, 1.0)) {
  seu <- FindClusters(seu, resolution=res, verbose=FALSE)
  n <- length(unique(seu[[paste0("RNA_snn_res.",res)]]))
  cat(sprintf("  res=%.1f: %d clusters\n", res, n))
}

cat("\n=== Dataset mixing ===\n")
tab <- table(seu$seurat_clusters, seu$dataset)
tab_pct <- round(prop.table(tab, 1)*100, 1)
for(cl in rownames(tab_pct)) {
  max_pct <- max(tab_pct[cl,])
  max_ds <- colnames(tab_pct)[which.max(tab_pct[cl,])]
  cat(sprintf("  C%s: %5d, max %.0f%% %s\n", cl, sum(tab[cl,]), max_pct, max_ds))
}

seu$umap_1 <- Embeddings(seu, "umap")[,1]
seu$umap_2 <- Embeddings(seu, "umap")[,2]
saveRDS(seu, file.path(OUTDIR, "seurat_npcs20.rds"))
cat(sprintf("\nSaved: seurat_npcs20.rds (%d cells, %d clusters)\n",
    ncol(seu), length(unique(seu$seurat_clusters))))
