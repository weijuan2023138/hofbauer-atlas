#!/usr/bin/env Rscript
# Generate shiny_data/ files from 11-dataset labeled object
library(Seurat); library(Matrix); library(dplyr)

INPUT  <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_11datasets.rds"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/shiny/shiny_data"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

cat("Loading 11-dataset object...\n")
seu <- readRDS(INPUT)
cat(sprintf("Cells: %d, Datasets: %d\n", ncol(seu), length(unique(seu$dataset))))

# ── 1. umap_meta.csv ──
cat("\n=== 1. Generating umap_meta.csv ===\n")
meta <- seu@meta.data
# Keep columns the shiny app expects
out_cols <- c("dataset","disease","disease_group","subtype","umap_1","umap_2",
              "trimester","nCount_RNA","nFeature_RNA","subtype_score")
out_cols <- intersect(out_cols, colnames(meta))
meta_out <- meta[, out_cols, drop=FALSE]
# Rename UMAP columns to match old format
colnames(meta_out)[colnames(meta_out)=="umap_1"] <- "UMAP_1"
colnames(meta_out)[colnames(meta_out)=="umap_2"] <- "UMAP_2"
write.csv(meta_out, file.path(OUTDIR, "umap_meta.csv"), row.names=TRUE)
cat(sprintf("  umap_meta.csv: %d rows, cols=%s\n", nrow(meta_out), 
    paste(colnames(meta_out), collapse=", ")))

# ── 2. expr_full.rds ──
cat("\n=== 2. Generating expr_full.rds ===\n")
DefaultAssay(seu) <- "RNA"
# Get log-normalized data
expr <- GetAssayData(seu, layer="data")
cat(sprintf("  Expression matrix: %d genes x %d cells\n", nrow(expr), ncol(expr)))
saveRDS(expr, file.path(OUTDIR, "expr_full.rds"))
cat("  expr_full.rds saved\n")

# ── 3. DEG files (Wilcoxon, same contrasts as old) ──
cat("\n=== 3. Generating DEG files ===\n")

# Define disease groups and their comparisons (matching old shiny logic)
# Old comparisons:
# - Preterm_vs_NormalLate: Preterm Labor / Term Labor vs Normal 3rd trimester / Preeclampsia
# - Preterm_vs_NormalLate_noTL: same but excluding Term Labor cells
# - PE_vs_NormalLate_noTL: Preeclampsia vs Normal 3rd trimester / Preeclampsia
# - Miscarriage_vs_NormalEarly: Miscarriage / Normal vs Normal 1st trimester
# - Infection_vs_NormalEarly: Infection vs Normal 1st trimester

# Need to set up identity classes for comparisons
# We'll use disease_group as identity

Idents(seu) <- seu$disease_group

# Helper function for DEG
run_DEG <- function(seu, ident.1, ident.2, outfile, label.1="Disease", label.2="Control") {
  if(!ident.1 %in% levels(Idents(seu)) || !ident.2 %in% levels(Idents(seu))) {
    cat(sprintf("  SKIP %s: missing group (%s or %s)\n", outfile, ident.1, ident.2))
    return(NULL)
  }
  n1 <- sum(Idents(seu) == ident.1)
  n2 <- sum(Idents(seu) == ident.2)
  if(n1 < 3 || n2 < 3) {
    cat(sprintf("  SKIP %s: too few cells (%d vs %d)\n", outfile, n1, n2))
    return(NULL)
  }
  cat(sprintf("  %s vs %s (%d vs %d cells)...\n", ident.1, ident.2, n1, n2))
  deg <- FindMarkers(seu, ident.1=ident.1, ident.2=ident.2,
      logfc.threshold=0, min.pct=0.1, verbose=FALSE)
  deg$gene <- rownames(deg)
  deg <- deg[, c("gene","p_val","avg_log2FC","pct.1","pct.2","p_val_adj")]
  write.csv(deg, file.path(OUTDIR, outfile), row.names=FALSE)
  cat(sprintf("    %d DEGs saved\n", nrow(deg)))
  return(deg)
}

# DEG 1: Preterm vs NormalLate
run_DEG(seu, 
    ident.1="Preterm Labor / Term Labor",
    ident.2="Normal 3rd trimester / Preeclampsia",
    outfile="deg_Preterm_vs_NormalLate.csv")

# DEG 2: Preterm vs NormalLate (no TL - Preterm Labor only)
# Create temporary identity separating PT from TL
seu$pt_tmp <- as.character(seu$disease_group)
seu$pt_tmp[seu$dataset == "GSE333257" & seu$disease == "PTL/TL"] <- "Preterm_Labor_Only"
Idents(seu) <- seu$pt_tmp
run_DEG(seu,
    ident.1="Preterm_Labor_Only",
    ident.2="Normal 3rd trimester / Preeclampsia",
    outfile="deg_Preterm_vs_NormalLate_noTL.csv")
Idents(seu) <- seu$disease_group  # restore

# DEG 3: PE vs NormalLate (no TL)
run_DEG(seu,
    ident.1="Preeclampsia",
    ident.2="Normal 3rd trimester / Preeclampsia",
    outfile="deg_PE_vs_NormalLate_noTL.csv")

# DEG 4: Miscarriage vs NormalEarly
run_DEG(seu,
    ident.1="Miscarriage / Normal",
    ident.2="Normal 1st trimester",
    outfile="deg_Miscarriage_vs_NormalEarly.csv")

# DEG 5: Infection vs NormalEarly
run_DEG(seu,
    ident.1="Infection",
    ident.2="Normal 1st trimester",
    outfile="deg_Infection_vs_NormalEarly.csv")

# ── 4. module_scores.csv ──
cat("\n=== 4. Generating module_scores.csv ===\n")
# ECM and Immune module scores
ecm_genes <- c("COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL6A1",
    "FN1","FBN1","SPARC","MMP2","MMP9","TIMP1","TIMP2","LOX","LOXL1","POSTN")
immune_genes <- c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","CD74","CIITA",
    "CXCL16","CCL2","CCL3","CCL4","IL1B","TNF","NFKB1","STAT1","IRF1")

ecm_found <- ecm_genes[ecm_genes %in% rownames(seu)]
imm_found <- immune_genes[immune_genes %in% rownames(seu)]
cat(sprintf("  ECM genes: %d/%d, Immune genes: %d/%d\n",
    length(ecm_found), length(ecm_genes), length(imm_found), length(immune_genes)))

seu <- AddModuleScore(seu, features=list(ecm_found), name="ECM")
seu <- AddModuleScore(seu, features=list(imm_found), name="Immune")

mod_df <- data.frame(
    ECM_score = seu$ECM1,
    Immune_score = seu$Immune1,
    subtype = seu$subtype,
    disease_group = seu$disease_group,
    row.names = colnames(seu)
)
write.csv(mod_df, file.path(OUTDIR, "module_scores.csv"))
cat(sprintf("  module_scores.csv: %d rows\n", nrow(mod_df)))

# ── 5. tf_expr.rds ──
cat("\n=== 5. Generating tf_expr.rds ===\n")
# Key TFs for Hofbauer atlas
tf_genes <- c("CEBPA","STAT3","STAT1","NFKB1","RELB","JUN","FOS","IRF1",
    "IRF8","KLF4","ID2","MAFB","SPI1","RUNX1","MAF","MITF","TFEB",
    "PPARG","NR1H3","FLI1","EGR1","EGR3")
tf_found <- tf_genes[tf_genes %in% rownames(seu)]
cat(sprintf("  TF genes found: %d/%d\n", length(tf_found), length(tf_genes)))
tf_expr <- GetAssayData(seu, layer="data")[tf_found, , drop=FALSE]
saveRDS(tf_expr, file.path(OUTDIR, "tf_expr.rds"))
cat(sprintf("  tf_expr.rds: %d TFs x %d cells\n", nrow(tf_expr), ncol(tf_expr)))

# ── 6. subtype_proportions.csv ──
cat("\n=== 6. Generating subtype_proportions.csv ===\n")
prop_tbl <- as.data.frame(prop.table(table(seu$subtype, seu$disease_group), margin=2))
colnames(prop_tbl) <- c("Subtype","Disease","Proportion")
write.csv(prop_tbl, file.path(OUTDIR, "subtype_proportions.csv"), row.names=FALSE)
cat(sprintf("  subtype_proportions.csv: %d rows\n", nrow(prop_tbl)))

cat("\n=== All shiny data generated! ===\n")
cat(sprintf("Output: %s/\n", OUTDIR))
system(paste("ls -lh", OUTDIR))
