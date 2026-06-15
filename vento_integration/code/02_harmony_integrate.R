#!/usr/bin/env Rscript
# Harmony integration for 10 datasets (-GSE183338, +GSE329173, +GSE298602)
library(Seurat); library(harmony); library(anndata); library(Matrix)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/all_hofbauer_10datasets.h5ad"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

message("Loading...")
adata <- read_h5ad(INPUT)
message(sprintf("Loaded: %d cells, %d genes", nrow(adata$obs), ncol(adata$X)))

message("Converting to Seurat...")
mat <- t(as.matrix(adata$X))
rownames(mat) <- adata$var_names; colnames(mat) <- adata$obs_names
seu <- CreateSeuratObject(counts=mat, project='Hofbauer_Atlas_10datasets')
for (col in colnames(adata$obs)) seu[[col]] <- adata$obs[[col]]

seu <- NormalizeData(seu, normalization.method="LogNormalize", scale.factor=10000)
seu <- FindVariableFeatures(seu, nfeatures=2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs=50)

message("Running Harmony...")
seu <- RunHarmony(seu, group.by.vars="dataset", dims.use=1:30, theta=2)
seu <- RunUMAP(seu, reduction="harmony", dims=1:30)
seu <- FindNeighbors(seu, reduction="harmony", dims=1:30)
for (res in c(0.1, 0.15, 0.2, 0.3, 0.5)) seu <- FindClusters(seu, resolution=res)

saveRDS(seu, file.path(OUTDIR, "seurat_harmony_10datasets.rds"))
message("Saved: seurat_harmony_10datasets.rds")
message(sprintf("Clusters (res=0.15): %d", length(unique(seu$RNA_snn_res.0.15))))
print(table(seu$dataset))
message("Done!")
