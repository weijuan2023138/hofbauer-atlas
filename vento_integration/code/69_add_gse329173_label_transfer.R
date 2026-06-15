#!/usr/bin/env Rscript
# Add GSE329173 (337 cells, Severe PE) to labeled 10-dataset Atlas via label transfer
# Produces 11-dataset complete object with subtype, disease, trimester
library(Seurat); library(anndata); library(Matrix); library(dplyr)

OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"

# ── 1. Load reference (labeled 10-dataset) ──
cat("Loading reference...\n")
ref <- readRDS(file.path(OUTDIR, "seurat_labeled_10datasets.rds"))
cat(sprintf("Reference: %d cells, %d genes, %d subtypes\n", 
    ncol(ref), nrow(ref), length(unique(ref$subtype_pred))))

# ── 2. Load GSE329173 from h5ad ──
cat("Loading GSE329173 from h5ad...\n")
all_ad <- read_h5ad(file.path(OUTDIR, "all_hofbauer_10datasets.h5ad"))
gse_mask <- all_ad$obs$dataset == "gse329173"
query_ad <- all_ad[gse_mask, ]
cat(sprintf("GSE329173: %d cells\n", nrow(query_ad$obs)))

# Convert to Seurat
qmat <- t(as.matrix(query_ad$X))
rownames(qmat) <- query_ad$var_names
colnames(qmat) <- make.unique(paste0("GSE329173_", query_ad$obs_names), sep="_")
query <- CreateSeuratObject(counts=qmat, project="GSE329173")
query$dataset <- "GSE329173"
cat(sprintf("Query Seurat: %d cells, %d genes\n", ncol(query), nrow(query)))

# ── 3. Find common genes ──
common_genes <- intersect(rownames(ref), rownames(query))
cat(sprintf("Common genes: %d\n", length(common_genes)))
ref_sub <- ref[common_genes, ]
query_sub <- query[common_genes, ]

# ── 4. Normalize both ──
ref_sub <- NormalizeData(ref_sub, normalization.method="LogNormalize", scale.factor=10000)
ref_sub <- FindVariableFeatures(ref_sub, nfeatures=2000)
query_sub <- NormalizeData(query_sub, normalization.method="LogNormalize", scale.factor=10000)
query_sub <- FindVariableFeatures(query_sub, nfeatures=2000)

# ── 5. Label transfer ──
n_dims <- min(20, ncol(Embeddings(ref_sub, "pca")))
cat(sprintf("Using %d PCA dims for transfer...\n", n_dims))
cat("Finding transfer anchors...\n")
anchors <- FindTransferAnchors(reference=ref_sub, query=query_sub, 
    dims=1:n_dims, reference.reduction="pca", verbose=FALSE)
cat(sprintf("Anchors: %d\n", nrow(anchors@anchors)))

cat("Transferring subtype labels...\n")
predictions <- TransferData(anchorset=anchors, refdata=ref$subtype_pred, dims=1:n_dims, verbose=FALSE)
query$subtype_pred <- predictions$predicted.id
query$subtype_score <- predictions$prediction.score.max

cat("GSE329173 subtype distribution:\n")
print(table(query$subtype_pred))

# ── 6. Project onto reference UMAP ──
cat("Projecting onto reference UMAP...\n")
# Increase future globals limit for MapQuery
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB
ref_sub <- RunUMAP(ref_sub, reduction="harmony", dims=1:n_dims, return.model=TRUE, verbose=FALSE)
query_sub <- MapQuery(anchorset=anchors, reference=ref_sub, query=query_sub,
    refdata=ref$subtype_pred, reference.reduction="pca", 
    reduction.model="umap", verbose=FALSE)

query$umap_1 <- Embeddings(query_sub, "ref.umap")[,1]
query$umap_2 <- Embeddings(query_sub, "ref.umap")[,2]
cat(sprintf("UMAP projected: range x=[%.2f, %.2f], y=[%.2f, %.2f]\n",
    min(query$umap_1), max(query$umap_1), min(query$umap_2), max(query$umap_2)))

# ── 7. Add disease/disease_group metadata ──
query$disease <- "Severe_PE"
query$disease_group <- "Severe Preeclampsia"
query$trimester <- "Late"

# Add same to reference (was missing)
ref$disease <- recode(ref$dataset,
    "E-MTAB-6701"="Normal_1st", "E-MTAB-12421"="Normal_1st",
    "E-MTAB-12795"="Normal/Listeria/Toxoplasma/Malaria",
    "GSE214607"="RM/NC", "UCSF Li 2026"="Normal",
    "GSE290578"="Normal/PE", "GSE298602"="PE/Control",
    "GSE333257"="PTL/TL", "GSE298119"="PE", "GSE173193"="PE"
)
ref$disease_group <- recode(ref$dataset,
    "E-MTAB-6701"="Normal 1st trimester", "E-MTAB-12421"="Normal 1st trimester",
    "E-MTAB-12795"="Infection", "GSE214607"="Miscarriage / Normal",
    "UCSF Li 2026"="Normal 1st/2nd/Term",
    "GSE290578"="Normal 3rd trimester / Preeclampsia",
    "GSE298602"="Preeclampsia / Control",
    "GSE333257"="Preterm Labor / Term Labor",
    "GSE298119"="Preeclampsia", "GSE173193"="Preeclampsia"
)
ref$trimester <- recode(ref$dataset,
    "E-MTAB-6701"="Early", "E-MTAB-12421"="Early",
    "E-MTAB-12795"="Early", "GSE214607"="Early",
    "UCSF Li 2026"="Mid",
    "GSE290578"="Late", "GSE298602"="Late", "GSE333257"="Late",
    "GSE298119"="Late", "GSE173193"="Late"
)

# ── 8. Ensure consistent metadata columns ──
ref$subtype <- ref$subtype_pred
query$subtype <- query$subtype_pred

keep_cols <- c("orig.ident","nCount_RNA","nFeature_RNA","dataset",
    "disease","disease_group","trimester","subtype","subtype_score",
    "umap_1","umap_2")
ref@meta.data <- ref@meta.data[, keep_cols]
query@meta.data <- query@meta.data[, keep_cols]

# ── 9. Merge ──
cat("\nMerging reference + GSE329173...\n")
combined <- merge(ref, y=query, add.cell.ids=c("ref","gse329173"))
cat(sprintf("Combined: %d cells, %d datasets\n", ncol(combined), length(unique(combined$dataset))))
cat("Per dataset:\n")
print(sort(table(combined$dataset)))

# ── 10. Save ──
outpath <- file.path(OUTDIR, "seurat_labeled_11datasets.rds")
saveRDS(combined, outpath)
cat(sprintf("\nSaved: %s (%d cells, %d datasets)\n", outpath, ncol(combined), length(unique(combined$dataset))))

# ── 11. Quick UMAP sanity check ──
library(ggplot2)
set.seed(42); idx <- sample(ncol(combined), min(8000, ncol(combined)))
p <- ggplot(combined@meta.data[idx,], aes(x=umap_1, y=umap_2, color=subtype)) +
    geom_point(size=0.3, alpha=0.7) +
    scale_color_manual(values=c(
        "Pro-inflammatory"="#C62828","MHCII+ Antigen-presenting"="#E65100",
        "Homeostatic"="#1565C0","PRKN+ Autophagy"="#6A1B9A",
        "Vascular remodeling"="#2E7D32","MKI67+ Proliferating"="#455A64")) +
    theme_minimal() + labs(title="11-dataset Atlas — Subtype UMAP")
ggsave(file.path(OUTDIR, "figures/check_11dataset_subtype_UMAP.png"), p, w=9, h=7, dpi=150, bg="white")
cat("Sanity check UMAP saved\n")

cat("\nDone!\n")
