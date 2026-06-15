#!/usr/bin/env Rscript
# Assign per-cell disease labels via barcode matching, then regenerate FigB
library(Seurat); library(anndata); library(ggplot2); library(patchwork); library(dplyr)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
seu$UMAP_1 <- Embeddings(seu, "umap")[,1]
seu$UMAP_2 <- Embeddings(seu, "umap")[,2]

# ── Load per-cell disease labels from original h5ad ──
cat("Loading per-cell disease labels...\n")

get_disease_labels <- function(h5ad_path, cell_type_col, label_val, disease_col) {
  ad <- read_h5ad(h5ad_path)
  if(!is.null(cell_type_col)) {
    ad <- ad[ad$obs[[cell_type_col]] == label_val, ]
  }
  labels <- setNames(as.character(ad$obs[[disease_col]]), ad$obs_names)
  cat(sprintf("  %s: %d cells, groups: %s\n", basename(h5ad_path),
      length(labels), paste(names(sort(table(labels), decreasing=T))[1:4], collapse=", ")))
  labels
}

BASE <- "/home/weijuan/文档/胎盘单细胞数据"

gse290578_labels <- get_disease_labels(
  file.path(BASE, "results/phase3_gse290578/gse290578_reclassified.h5ad"),
  "cell_type_new", "Hofbauer", "condition")

gse298602_labels <- get_disease_labels(
  file.path(BASE, "ucsf_integration/results/classification/gse298602_all_hofbauer.h5ad"),
  NULL, NULL, "disease")

hoo_labels <- get_disease_labels(
  file.path(BASE, "processed/hoo_2024_reclassified.h5ad"),
  "cell_type_new", "Hofbauer", "condition")

gse333257_labels <- get_disease_labels(
  file.path(BASE, "processed/my_cohort_processed.h5ad"),
  "cell_type_fine", "Hofbauer", "condition")

# ── Assign labels by constructing full named vector ──
all_cells <- colnames(seu)
disease_vec <- rep(NA_character_, length(all_cells))
names(disease_vec) <- all_cells

assign_labels_vec <- function(vec, seu, ds_name, labels) {
  cells <- colnames(seu)[seu$dataset == ds_name]
  matched <- intersect(cells, names(labels))
  cat(sprintf("%s: %d/%d cells matched\n", ds_name, length(matched), length(cells)))
  vec[matched] <- labels[matched]
  vec
}

disease_vec <- assign_labels_vec(disease_vec, seu, "GSE290578", gse290578_labels)
disease_vec <- assign_labels_vec(disease_vec, seu, "GSE298602", gse298602_labels)
disease_vec <- assign_labels_vec(disease_vec, seu, "E-MTAB-12795", hoo_labels)
disease_vec <- assign_labels_vec(disease_vec, seu, "GSE333257", gse333257_labels)

# Fill remaining with dataset name
na_cells <- is.na(disease_vec)
disease_vec[na_cells] <- as.character(seu$dataset[na_cells])

# Add to Seurat (simple named vector assignment)
seu$disease_detail <- disease_vec[colnames(seu)]

# ── Map to clean disease groups ──
# Standardize labels
label_map <- c(
  # GSE290578
  "Normal"="Normal", "PE"="PE",
  # GSE298602
  "PreE_SF"="PE", "Control"="Normal", "gHTN"="PE",
  # Hoo
  "normal"="Infection", "toxoplasmosis"="Infection",
  "listeriosis"="Infection", "Plasmodium malariae malaria"="Infection",
  # GSE333257
  "PTL"="Preterm", "TL"="Term", "PTNL"="Preterm",
  # Other datasets
  "E-MTAB-12421"="Normal", "E-MTAB-6701"="Normal",
  "UCSF Li 2026"="Normal",
  "GSE173193"="PE", "GSE298119"="PE",
  "GSE214607"="Miscarriage"
)

seu$disease_final <- factor(label_map[seu$disease_detail],
  levels=c("Normal","PE","Miscarriage","Infection","Preterm","Term"))

# Colors
disease_cols <- c("Normal"="#4575B4","PE"="#FC8D59",
  "Miscarriage"="#D73027","Infection"="#FDB462",
  "Preterm"="#E41A1C","Term"="#2D8B57")

cat("\nFinal disease groups:\n")
print(table(seu$disease_final))

# ── Generate FigB ──
set.seed(42); idx <- sample(ncol(seu), 20000)
tpub <- theme_bw(base_size=11) +
  theme(panel.grid=element_line(color="grey92",linewidth=0.2),
        panel.border=element_rect(color="black",fill=NA,linewidth=0.5),
        axis.text=element_blank(), axis.ticks=element_blank(),
        legend.position="right", legend.title=element_blank(),
        legend.text=element_text(size=9), legend.key=element_blank(),
        plot.title=element_text(size=13,face="bold"))

pB <- ggplot(seu@meta.data[idx,], aes(UMAP_1, UMAP_2, color=disease_final)) +
  geom_point(size=0.12, alpha=0.8) + scale_color_manual(values=disease_cols) +
  labs(title="By Disease Group (per-cell)", x="UMAP 1", y="UMAP 2") + tpub +
  guides(color=guide_legend(override.aes=list(size=3), ncol=1))
ggsave(file.path(FIGDIR,"FigB_disease_percell.png"), pB, w=10, h=7, dpi=300, bg="white")
cat("Saved: FigB_disease_percell.png\n")
