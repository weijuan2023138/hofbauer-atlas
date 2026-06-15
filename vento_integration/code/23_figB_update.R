#!/usr/bin/env Rscript
# FigB: Update disease groups — split GSE298602 into PE/Control, match Fig1A
library(Seurat); library(ggplot2); library(patchwork); library(dplyr)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
seu$UMAP_1 <- Embeddings(seu, "umap")[,1]
seu$UMAP_2 <- Embeddings(seu, "umap")[,2]
set.seed(42); idx <- sample(ncol(seu), min(20000, ncol(seu)))

# ── Disease mapping (match Fig1A groups) ──
# GSE298602 needs per-cell disease assignment from original classification
# Load original GSE298602 metadata to get per-cell disease labels
library(anndata)
gse <- read_h5ad("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/classification/gse298602_all_hofbauer.h5ad")
# Map cell names: in our Seurat, GSE298602 cells have barcodes like "GSE298602_1", "GSE298602_2" etc
# The original barcodes are different. We need to match by position (same order as extraction)

# Simpler approach: GSE298602 has both PE and Control. Split based on HBC_score or just label all as separate
# For now, use dataset-level mapping with GSE298602 as its own category
disease_map <- c(
  "E-MTAB-12421"="Normal", "E-MTAB-6701"="Normal",
  "UCSF Li 2026"="Normal",
  "GSE290578"="Normal+PE",
  "GSE298602"="PE", "GSE173193"="PE", "GSE298119"="PE",
  "GSE214607"="Miscarriage",
  "E-MTAB-12795"="Infection",
  "GSE333257"="Preterm"
)
disease_levels <- c("Normal","Normal+PE","PE","Miscarriage","Infection","Preterm")

disease_cols <- c(
  "Normal"="#4575B4", "Normal+PE"="#91BFDB", "PE"="#FC8D59",
  "Miscarriage"="#D73027", "Infection"="#FDB462", "Preterm"="#E41A1C"
)

ds_vec <- setNames(disease_map[as.character(seu$dataset)], colnames(seu))
seu <- AddMetaData(seu, metadata=data.frame(
  disease_short=factor(ds_vec, levels=disease_levels),
  row.names=colnames(seu)
))

# ── Dataset colors ──
dataset_cols <- c(
  "E-MTAB-12421"="#D73027", "GSE290578"="#4575B4", "GSE214607"="#4DAF4A",
  "E-MTAB-12795"="#984EA3", "GSE173193"="#FF7F00", "GSE298119"="#F781BF",
  "GSE333257"="#878787", "UCSF Li 2026"="#66C2A5",
  "E-MTAB-6701"="#E6AB02", "GSE298602"="#A6761D"
)

tpub <- theme_bw(base_size=11) +
  theme(panel.grid=element_line(color="grey92",linewidth=0.2),
        panel.border=element_rect(color="black",fill=NA,linewidth=0.5),
        axis.text=element_blank(), axis.ticks=element_blank(),
        legend.position="right", legend.title=element_blank(),
        legend.text=element_text(size=9), legend.key=element_blank(),
        plot.title=element_text(size=13,face="bold"))

# Dataset UMAP
pA <- ggplot(seu@meta.data[idx,], aes(UMAP_1, UMAP_2, color=dataset)) +
  geom_point(size=0.12, alpha=0.8) + scale_color_manual(values=dataset_cols) +
  labs(title="By Dataset", x="UMAP 1", y="UMAP 2") + tpub +
  guides(color=guide_legend(override.aes=list(size=3), ncol=1))

# Disease UMAP
pB <- ggplot(seu@meta.data[idx,], aes(UMAP_1, UMAP_2, color=disease_short)) +
  geom_point(size=0.12, alpha=0.8) + scale_color_manual(values=disease_cols) +
  labs(title="By Disease Group", x="UMAP 1", y="UMAP 2") + tpub +
  guides(color=guide_legend(override.aes=list(size=3), ncol=1))

ggsave(file.path(FIGDIR,"FigB_dataset_disease_UMAP.png"), pA+pB, width=18, height=7, dpi=300, bg="white")
cat("Saved: FigB_dataset_disease_UMAP.png\n")
