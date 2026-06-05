#!/usr/bin/env Rscript
# Step 3: Validate mapping between new clusters and existing 7 subtypes
# Input: unified_with_ucsf.rds
# Output: mapping_validation.csv + DEG tables for manual annotation

library(Seurat)
library(dplyr)

# Paths
INPUT <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/unified_with_ucsf.rds'
OUTPUT_DIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results'

cat("=" , rep("=", 59), "\n", sep="")
cat("Step 3: Validate Cluster-Subtype Mapping\n")
cat("=" , rep("=", 59), "\n", sep="")

# 1. Load integrated data
cat("\n[1/5] Loading integrated data...\n")
combined <- readRDS(INPUT)
cat("  Total cells:", ncol(combined), "\n")
cat("  Clusters:", length(unique(combined$seurat_clusters)), "\n")

# 2. Check mapping consistency (only for existing Atlas cells)
cat("\n[2/5] Checking cluster-subtype mapping...\n")

# Separate existing and UCSF cells
existing_mask <- combined$dataset != 'UCSF_Li_2026'
ucsf_mask <- combined$dataset == 'UCSF_Li_2026'

# Create mapping table for existing Atlas cells only
existing_clusters <- combined$seurat_clusters[existing_mask]
existing_subtypes <- combined$subtype_final[existing_mask]

# Remove NA subtypes
valid_mask <- !is.na(existing_subtypes)
existing_clusters <- existing_clusters[valid_mask]
existing_subtypes <- existing_subtypes[valid_mask]

cat("  Existing Atlas cells with valid subtypes:", sum(valid_mask), "\n")

# Create mapping table
mapping <- table(Cluster = existing_clusters, Subtype = existing_subtypes)
mapping_pct <- prop.table(mapping, margin = 1)

# Build results for ALL clusters
all_clusters <- sort(unique(as.character(combined$seurat_clusters)))
results_list <- list()

for (cl in all_clusters) {
  total <- sum(combined$seurat_clusters == cl)
  ucsf <- sum(combined$seurat_clusters == cl & ucsf_mask)
  existing <- sum(combined$seurat_clusters == cl & existing_mask)
  
  if (cl %in% rownames(mapping_pct)) {
    row_vals <- mapping_pct[cl, ]
    # Check if all values are NaN (0 existing cells with valid subtypes)
    if (all(is.nan(row_vals))) {
      best_match <- "UCSF-only"
      match_pct <- NA
    } else {
      best_idx <- which.max(row_vals)
      if (length(best_idx) == 0) {
        best_match <- "UCSF-only"
        match_pct <- NA
      } else {
        best_match <- names(row_vals)[best_idx]
        match_pct <- row_vals[best_idx]
      }
    }
  } else {
    best_match <- "UCSF-only"
    match_pct <- NA
  }
  
  results_list[[cl]] <- data.frame(
    Cluster = cl,
    Best_Match = best_match,
    Match_Pct = match_pct,
    Total_Cells = total,
    UCSF_Cells = ucsf,
    Existing_Cells = existing,
    stringsAsFactors = FALSE
  )
}

mapping_results <- do.call(rbind, results_list)
rownames(mapping_results) <- NULL

# Print results
cat("\nCluster-Subtype Mapping:\n")
cat(rep("-", 85), "\n", sep="")
cat(sprintf("%-8s %-25s %-10s %-8s %-8s %-8s\n", 
    "Cluster", "Best Match", "Match%", "Total", "UCSF", "Existing"))
cat(rep("-", 85), "\n", sep="")
for (i in 1:nrow(mapping_results)) {
  match_str <- ifelse(is.na(mapping_results$Match_Pct[i]), "N/A", 
                      sprintf("%.1f%%", mapping_results$Match_Pct[i] * 100))
  cat(sprintf("%-8s %-25s %-10s %-8d %-8d %-8d\n",
      mapping_results$Cluster[i],
      mapping_results$Best_Match[i],
      match_str,
      mapping_results$Total_Cells[i],
      mapping_results$UCSF_Cells[i],
      mapping_results$Existing_Cells[i]))
}

# 3. Identify clusters needing attention
cat("\n[3/5] Identifying clusters needing attention...\n")
unclear <- mapping_results$Cluster[!is.na(mapping_results$Match_Pct) & mapping_results$Match_Pct < 0.8]
ucsf_only <- mapping_results$Cluster[mapping_results$Best_Match == "UCSF-only"]

if (length(unclear) > 0) {
  cat("  Clusters with <80% mapping:", paste(unclear, collapse=", "), "\n")
}
if (length(ucsf_only) > 0) {
  cat("  UCSF-only clusters:", paste(ucsf_only, collapse=", "), "\n")
}
if (length(unclear) == 0 && length(ucsf_only) == 0) {
  cat("  All clusters map clearly (>=80%).\n")
}

# 4. Generate DEG tables for manual annotation
cat("\n[4/5] Generating DEG tables for manual annotation...\n")
cat("  This may take 5-10 minutes...\n")

# Join layers first (Seurat v5 requirement)
combined <- JoinLayers(combined)

# Find markers for each cluster
am <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.5, test.use = 't')

# Filter non-coding/ribosomal/mitochondrial
am <- am[!grepl("^(RPL|RPS|MT-|MALAT1|LINC|AC0|AL|AP0|RP11|RP1-|CTD-|XXbac|Z9|THUMPD)", am$gene), ]

# Save full DEG table
write.csv(am, file.path(OUTPUT_DIR, 'all_markers.csv'), row.names = FALSE)
cat("  Saved all_markers.csv\n")

# Output top 15 markers per cluster
cat("\nTop 15 markers per cluster:\n")
cat(rep("=", 80), "\n", sep="")
sink(file.path(OUTPUT_DIR, 'top_markers_per_cluster.txt'))
for (cl in sort(unique(am$cluster))) {
  dd <- am[am$cluster == cl, ]
  dd <- dd[order(-dd$avg_log2FC), ]
  dd <- head(dd, 15)
  
  cat(sprintf("\nC%s (n=%d cells, %d UCSF):\n", 
      cl, 
      sum(combined$seurat_clusters == cl),
      sum(combined$seurat_clusters == cl & ucsf_mask)))
  cat(rep("-", 60), "\n", sep="")
  for (i in 1:nrow(dd)) {
    cat(sprintf("  %-15s  log2FC=%-6.2f  pct=%.0f%%  padj=%.1e\n",
        dd$gene[i], dd$avg_log2FC[i], dd$pct.1[i] * 100, dd$p_val_adj[i]))
  }
}
sink()
cat("  Saved top_markers_per_cluster.txt\n")

# 5. Generate key gene expression matrix
cat("\n[5/5] Generating key gene expression matrix...\n")
key_genes <- c('FOLR2', 'CD163', 'MRC1', 'TNF', 'IL1B', 'CXCL8', 
               'CCL13', 'AIF1', 'CLIC1',
               'C1QA', 'C1QB', 'FCGR3A', 'HLA-DRA',
               'PRKN', 'C9', 'SOX5',
               'MKI67', 'BUB1B', 'KIF4A',
               'SPP1', 'DAB2', 'S100A4',
               'NFKBIZ', 'PTGS2', 'DUSP2')

# Filter to genes present in the data
key_genes_present <- key_genes[key_genes %in% rownames(combined)]
cat("  Key genes present:", length(key_genes_present), "/", length(key_genes), "\n")

# Calculate average expression per cluster
avg_expr <- AverageExpression(combined, features = key_genes_present, assay = 'RNA')
avg_df <- as.data.frame(avg_expr$RNA)
avg_df$Gene <- rownames(avg_df)

# Save
write.csv(avg_df, file.path(OUTPUT_DIR, 'key_gene_expression.csv'), row.names = FALSE)
cat("  Saved key_gene_expression.csv\n")

# Print key gene expression
cat("\nKey gene expression by cluster:\n")
cat(rep("=", 80), "\n", sep="")
print(round(avg_expr$RNA, 2))

# 6. Save mapping results
write.csv(mapping_results, file.path(OUTPUT_DIR, 'mapping_validation.csv'), row.names = FALSE)
cat("\n  Saved mapping_validation.csv\n")

# 7. Summary
cat("\n" , rep("=", 60), "\n", sep="")
cat("Summary:\n")
cat("  Total cells:", ncol(combined), "\n")
cat("  UCSF cells:", sum(ucsf_mask), "\n")
cat("  Clusters:", length(unique(combined$seurat_clusters)), "\n")
cat("  Clear mapping (>=80%):", sum(!is.na(mapping_results$Match_Pct) & mapping_results$Match_Pct >= 0.8), "\n")
cat("  Unclear mapping (<80%):", sum(!is.na(mapping_results$Match_Pct) & mapping_results$Match_Pct < 0.8), "\n")
cat("  UCSF-only clusters:", sum(mapping_results$Best_Match == "UCSF-only"), "\n")

cat("\nOutput files:\n")
cat("  1. mapping_validation.csv - Cluster-subtype mapping\n")
cat("  2. top_markers_per_cluster.txt - DEG tables for annotation\n")
cat("  3. key_gene_expression.csv - Key gene expression matrix\n")
cat("  4. all_markers.csv - Full DEG table\n")

cat("\n", rep("=", 60), "\n", sep="")
cat("Done! Please review the DEG tables and annotate clusters.\n")
cat("Next step: 04_apply_annotation.R (after manual annotation)\n")
cat(rep("=", 60), "\n", sep="")
