#!/usr/bin/env Rscript
# Final Atlas: extract matrices → merge → single normalization → Harmony
library(Seurat); library(harmony); library(anndata); library(Matrix); library(dplyr)

OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

# ── 1. Old Atlas ──
cat("Loading old Atlas...\n")
old <- readRDS("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/Hofbauer_Atlas_Final.rds")
old <- subset(old, dataset != "gse183338")
old$dataset <- recode(old$dataset,
  "Arutyunyan"="E-MTAB-12421", "GSE290578"="GSE290578",
  "gse214607"="GSE214607", "hoo_2024"="E-MTAB-12795",
  "gse173193"="GSE173193", "gse298119"="GSE298119",
  "my_preterm_cohort"="GSE333257", "UCSF_Li_2026"="UCSF Li 2026"
)
# Get raw counts
old_counts <- GetAssayData(old, layer="counts")
old_meta <- old@meta.data[, c("dataset"), drop=FALSE]
cat(sprintf("Old: %d cells\n", ncol(old)))

# ── 2. New datasets (use expm1 to get approx raw counts) ──
new_paths <- c(
  "E-MTAB-6701"="/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/classification/vento_tormo_reclassified_hofbauer.h5ad",
  "GSE298602"="/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/classification/gse298602_all_hofbauer.h5ad"
)

new_mats <- list()
new_metas <- list()

for(ds in names(new_paths)) {
  cat(sprintf("Loading %s...\n", ds))
  ad <- read_h5ad(new_paths[ds])
  mat <- expm1(as.matrix(ad$X))
  colnames(mat) <- paste0(ds, "_", 1:ncol(mat))
  # Use original gene names, make unique if needed
  gn <- ad$var_names
  gn <- make.unique(gn, sep=".")
  rownames(mat) <- gn[1:nrow(mat)]
  
  new_mats[[ds]] <- mat
  meta <- data.frame(dataset=rep(ds, ncol(mat)), row.names=NULL)
  rownames(meta) <- colnames(mat)
  new_metas[[ds]] <- meta
  cat(sprintf("  %d cells\n", ncol(mat)))
}

# ── 3. Merge on common genes ──
cat("\nFinding common genes...\n")
all_genes <- Reduce(intersect, c(list(rownames(old_counts)), lapply(new_mats, rownames)))
cat(sprintf("Common genes: %d\n", length(all_genes)))

old_sub <- old_counts[all_genes, ]
new_subs <- lapply(new_mats, function(m) m[all_genes, , drop=FALSE])

combined_mat <- cbind(old_sub, do.call(cbind, new_subs))
combined_meta <- rbind(old_meta, do.call(rbind, new_metas))
cat(sprintf("Combined: %d cells, %d genes\n", ncol(combined_mat), nrow(combined_mat)))

# ── 4. Create Seurat + normalize once ──
cat("\nCreating Seurat object...\n")
seu <- CreateSeuratObject(counts=combined_mat, project="Hofbauer_Atlas_10datasets",
                          meta.data=combined_meta)
seu <- NormalizeData(seu, normalization.method="LogNormalize", scale.factor=10000)
seu <- FindVariableFeatures(seu, nfeatures=2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs=30)

# ── 5. Harmony ──
cat("Running Harmony...\n")
seu <- RunHarmony(seu, group.by.vars="dataset", dims.use=1:30, theta=2)
seu <- RunUMAP(seu, reduction="harmony", dims=1:30)
seu <- FindNeighbors(seu, reduction="harmony", dims=1:30)

for(res in c(0.5, 0.8, 1.0)) {
  seu <- FindClusters(seu, resolution=res, verbose=FALSE)
  n <- length(unique(seu[[paste0("RNA_snn_res.",res)]]))
  cat(sprintf("  res=%.1f: %d clusters\n", res, n))
}

# ── 6. Mixing check ──
cat("\n=== Dataset mixing (res=0.5) ===\n")
tab <- table(seu$RNA_snn_res.0.5, seu$dataset)
tab_pct <- round(prop.table(tab, 1)*100, 1)
for(cl in rownames(tab_pct)) {
  max_pct <- max(tab_pct[cl,])
  max_ds <- colnames(tab_pct)[which.max(tab_pct[cl,])]
  cat(sprintf("  C%s: %5d cells, max %.0f%% %s\n", cl, sum(tab[cl,]), max_pct, max_ds))
}

# ── 7. Total per dataset ──
cat("\n=== Per-dataset count ===\n")
print(sort(table(seu$dataset), decreasing=TRUE))

# ── 8. Save ──
seu$umap_1 <- Embeddings(seu, "umap")[,1]
seu$umap_2 <- Embeddings(seu, "umap")[,2]
saveRDS(seu, file.path(OUTDIR, "seurat_final_10datasets.rds"))
cat(sprintf("\nSaved: %d cells\n", ncol(seu)))
