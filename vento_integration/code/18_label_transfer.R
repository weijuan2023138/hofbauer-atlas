#!/usr/bin/env Rscript
# Label transfer: old 6 subtypes → new 10-dataset Atlas
library(Seurat); library(dplyr)

OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"

# Load
cat("Loading reference (old Atlas)...\n")
ref <- readRDS("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/Hofbauer_Atlas_Final.rds")
ref <- RenameCells(ref, add.cell.id="ref")
cat(sprintf("Reference: %d cells, %d subtypes\n", ncol(ref), length(unique(ref$subtype))))

cat("Loading query (new Atlas)...\n")
query <- readRDS(file.path(OUTDIR, "seurat_npcs20.rds"))
cat(sprintf("Query: %d cells\n", ncol(query)))

# Find common genes
common <- intersect(rownames(ref), rownames(query))
cat(sprintf("Common genes: %d\n", length(common)))

ref <- ref[common, ]
query <- query[common, ]

# Normalize both the same way
ref <- NormalizeData(ref, normalization.method="LogNormalize", scale.factor=10000)
ref <- FindVariableFeatures(ref, nfeatures=2000)
query <- NormalizeData(query, normalization.method="LogNormalize", scale.factor=10000)
query <- FindVariableFeatures(query, nfeatures=2000)

# Transfer anchors
cat("Finding transfer anchors...\n")
anchors <- FindTransferAnchors(reference=ref, query=query, dims=1:30,
                                reference.reduction="pca", verbose=FALSE)
cat(sprintf("Anchors found: %d\n", nrow(anchors@anchors)))

# Predict subtypes
cat("Predicting subtypes...\n")
predictions <- TransferData(anchorset=anchors, refdata=ref$subtype, dims=1:30, verbose=FALSE)
query$subtype_pred <- predictions$predicted.id
query$subtype_score <- predictions$prediction.score.max

# Summary
cat("\n=== Predicted subtypes ===\n")
for(s in sort(unique(query$subtype_pred))) {
  n <- sum(query$subtype_pred==s)
  cat(sprintf("  %-30s %5d (%.0f%%)\n", s, n, n/ncol(query)*100))
}

# Score distribution
cat(sprintf("\nPrediction score: mean=%.3f, median=%.3f\n",
    mean(query$subtype_score), median(query$subtype_score)))

# Compare with old
cat("\n=== OLD subtypes ===\n")
for(s in sort(unique(ref$subtype))) {
  n <- sum(ref$subtype==s)
  cat(sprintf("  %-30s %5d (%.0f%%)\n", s, n, n/ncol(ref)*100))
}

# Save
saveRDS(query, file.path(OUTDIR, "seurat_labeled_10datasets.rds"))
cat(sprintf("\nSaved: seurat_labeled_10datasets.rds (%d cells)\n", ncol(query)))
