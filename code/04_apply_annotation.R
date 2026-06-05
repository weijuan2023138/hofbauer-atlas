#!/usr/bin/env Rscript
# Step 4: Apply manual annotation and update Atlas
# Input: unified_with_ucsf.rds + user-annotated cluster_assignments.csv
# Output: unified_with_ucsf_annotated.rds (final Atlas)

library(Seurat)
library(dplyr)

# Paths
INPUT <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/unified_with_ucsf.rds'
ANNOTATION_FILE <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/cluster_assignments.csv'
OUTPUT_DIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results'

cat("=" , rep("=", 59), "\n", sep="")
cat("Step 4: Apply Manual Annotation\n")
cat("=" , rep("=", 59), "\n", sep="")

# 1. Load integrated data
cat("\n[1/4] Loading integrated data...\n")
combined <- readRDS(INPUT)
cat("  Total cells:", ncol(combined), "\n")
cat("  Clusters:", length(unique(combined$seurat_clusters)), "\n")

# 2. Load annotation
cat("\n[2/4] Loading annotation...\n")
if (!file.exists(ANNOTATION_FILE)) {
  cat("  ERROR: Annotation file not found!\n")
  cat("  Please create:", ANNOTATION_FILE, "\n")
  cat("  Format: CSV with columns 'cluster' and 'subtype'\n")
  cat("  Example:\n")
  cat("    cluster,subtype\n")
  cat("    0,Homeostatic\n")
  cat("    1,CCL13+ Regulatory\n")
  cat("    2,C1Q+ Complement\n")
  cat("    ...\n")
  stop("Annotation file not found")
}

annotation <- read.csv(ANNOTATION_FILE, stringsAsFactors = FALSE)
cat("  Annotation entries:", nrow(annotation), "\n")
cat("  Subtypes:", paste(unique(annotation$subtype), collapse=", "), "\n")

# 3. Apply annotation
cat("\n[3/4] Applying annotation...\n")

# Create named vector
annot_vec <- setNames(annotation$subtype, annotation$cluster)

# Map clusters to subtypes
subtype_vec <- annot_vec[as.character(combined$seurat_clusters)]

# Check for unmapped clusters
unmapped <- is.na(subtype_vec)
if (any(unmapped)) {
  cat("  WARNING:", sum(unmapped), "cells could not be mapped!\n")
  cat("  Unmapped clusters:", paste(unique(combined$seurat_clusters[unmapped]), collapse=", "), "\n")
  cat("  Setting unmapped cells to 'Unassigned'\n")
  subtype_vec[unmapped] <- 'Unassigned'
}

# Assign to Seurat object
names(subtype_vec) <- colnames(combined)
combined$subtype_final <- subtype_vec

# Verify
cat("\n  Subtype distribution:\n")
print(table(combined$subtype_final))

# 4. Save
cat("\n[4/4] Saving...\n")
output_file <- file.path(OUTPUT_DIR, 'unified_with_ucsf_annotated.rds')
saveRDS(combined, output_file)
cat("  Saved:", output_file, "\n")
cat("  File size:", round(file.size(output_file) / 1e6, 1), "MB\n")

# 5. Summary
cat("\n" , rep("=", 60), "\n", sep="")
cat("Summary:\n")
cat("  Total cells:", ncol(combined), "\n")
cat("  UCSF cells:", sum(combined$dataset == 'UCSF_Li_2026'), "\n")
cat("  Existing Atlas cells:", sum(combined$dataset != 'UCSF_Li_2026'), "\n")
cat("  Subtypes:", length(unique(combined$subtype_final)), "\n")

cat("\nSubtype distribution:\n")
subtype_counts <- table(combined$subtype_final)
for (st in names(subtype_counts)) {
  cat(sprintf("  %-25s: %d (%.1f%%)\n", st, subtype_counts[st], 
      subtype_counts[st] / sum(subtype_counts) * 100))
}

cat("\nUCSF subtype distribution:\n")
ucsf_mask <- combined$dataset == 'UCSF_Li_2026'
ucsf_counts <- table(combined$subtype_final[ucsf_mask])
for (st in names(ucsf_counts)) {
  cat(sprintf("  %-25s: %d (%.1f%%)\n", st, ucsf_counts[st], 
      ucsf_counts[st] / sum(ucsf_counts) * 100))
}

cat("\n", rep("=", 60), "\n", sep="")
cat("Done! Atlas updated with UCSF data.\n")
cat("Next steps:\n")
cat("  1. Update Shiny app with new Atlas\n")
cat("  2. Re-generate figures\n")
cat("  3. Write Methods section\n")
cat(rep("=", 60), "\n", sep="")
