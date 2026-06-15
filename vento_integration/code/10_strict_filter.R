#!/usr/bin/env Rscript
# Re-compute MAT_score and DIFF, apply strict Hofbauer filter
library(Seurat); library(dplyr)

INPUT  <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_harmony_10datasets.rds"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"

seu <- readRDS(INPUT)
seu <- JoinLayers(seu)
cat(sprintf("Raw: %d cells\n", ncol(seu)))

# Load classifier gene lists
gene_df <- read.csv("/home/weijuan/文档/胎盘单细胞数据/results/phase1_classifier/classifier_genes.csv")
hb_genes <- gene_df$gene[gene_df$direction=="HB_up"]
mat_genes <- gene_df$gene[gene_df$direction=="MAT_up"]

# Re-compute MAT_score
hb_found <- hb_genes[hb_genes %in% rownames(seu)]
mat_found <- mat_genes[mat_genes %in% rownames(seu)]
cat(sprintf("Classifier genes: HB=%d/%d, MAT=%d/%d\n", 
    length(hb_found), length(hb_genes), length(mat_found), length(mat_genes)))

seu <- AddModuleScore(seu, features=list(hb_found), name="HBC_recalc")
seu <- AddModuleScore(seu, features=list(mat_found), name="MAT_recalc")

# DIFF = HBC - MAT
seu$DIFF_recalc <- seu$HBC_recalc1 - seu$MAT_recalc1

# Also add contaminant scores
contam <- c("CGA","CSH1","KRT8","KRT19","PECAM1","VIM","IGKC","JCHAIN","GNLY","LYZ","S100A8")
contam_found <- contam[contam %in% rownames(seu)]
if(length(contam_found) > 0) {
  seu <- AddModuleScore(seu, features=list(contam_found), name="Contam")
}

# Apply filters: strict Hofbauer
# DIFF > 0.32 AND HBC > MAT AND Mac_score > 0
strict_mask <- seu$DIFF_recalc > 0.32 & seu$Mac_score > 0
cat(sprintf("\nDIFF_recalc > 0.32: %d (%.1f%%)\n", sum(strict_mask), sum(strict_mask)/ncol(seu)*100))

# Also filter by contaminant score if available
if("Contam1" %in% colnames(seu@meta.data)) {
  contam_high <- seu$Contam1 > 0.5
  cat(sprintf("Contam_score > 0.5: %d\n", sum(contam_high)))
  strict_mask <- strict_mask & (!contam_high)
}

# Apply
seu_clean <- subset(seu, cells=colnames(seu)[strict_mask])
cat(sprintf("\nAfter strict filtering: %d cells\n", ncol(seu_clean)))

# Check contaminant markers in clean set
cat("\nContaminant markers in clean Hofbauer:\n")
for(g in c("FOLR2","CD163","CGA","KRT8","PECAM1","IGKC","HLA-DRA")) {
  if(g %in% rownames(seu_clean)) {
    expr <- FetchData(seu_clean, vars=g)[,1]
    cat(sprintf("  %s: mean=%.3f, max=%.3f\n", g, mean(expr), max(expr)))
  }
}

# Save
saveRDS(seu_clean, file.path(OUTDIR, "seurat_hofbauer_strict.rds"))
cat(sprintf("\nSaved: seurat_hofbauer_strict.rds (%d cells)\n", ncol(seu_clean)))

# Per-dataset retention
cat("\nPer-dataset retention:\n")
for(ds in sort(unique(seu$dataset))) {
  total <- sum(seu$dataset == ds)
  kept <- sum(seu_clean$dataset == ds)
  cat(sprintf("  %-25s %5d -> %5d (%.0f%% kept)\n", ds, total, kept, kept/total*100))
}
