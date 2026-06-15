#!/usr/bin/env Rscript
# FigB: Dataset + Disease UMAP for labeled 10-dataset Atlas
library(Seurat); library(ggplot2); library(patchwork); library(dplyr)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
OUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures/Fig1"

seu <- readRDS(INPUT)
seu$umap_1 <- Embeddings(seu, "umap")[,1]
seu$umap_2 <- Embeddings(seu, "umap")[,2]
set.seed(42); idx <- sample(1:ncol(seu), min(20000, ncol(seu)))

# ── Dataset colors ──
dataset_cols <- c(
  "E-MTAB-12421"="#D73027", "GSE290578"="#4575B4", "GSE214607"="#4DAF4A",
  "E-MTAB-12795"="#984EA3", "GSE173193"="#FF7F00", "GSE298119"="#F781BF",
  "GSE333257"="#878787", "UCSF Li 2026"="#66C2A5",
  "E-MTAB-6701"="#E6AB02", "GSE298602"="#A6761D"
)

# ── Disease group mapping ──
# Determine disease group for each dataset
disease_map <- c(
  "E-MTAB-12421"="Normal 1st", "E-MTAB-6701"="Normal 1st",
  "UCSF Li 2026"="Normal", "GSE290578"="Normal+PE",
  "GSE214607"="Miscarriage", "E-MTAB-12795"="Infection",
  "GSE173193"="PE", "GSE298119"="PE", "GSE298602"="PE/Control",
  "GSE333257"="Preterm"
)
# Create disease short names
ds_vec <- setNames(disease_map[as.character(seu$dataset)], colnames(seu))
seu <- AddMetaData(seu, metadata=data.frame(
  disease_short=factor(ds_vec, levels=c("Normal 1st","Normal","Normal+PE","PE","PE/Control",
                                         "Miscarriage","Infection","Preterm")),
  row.names=colnames(seu)
))

disease_cols <- c(
  "Normal 1st"="#4575B4", "Normal"="#66C2A5",
  "Normal+PE"="#91BFDB", "PE"="#FC8D59", "PE/Control"="#FC8D59",
  "Miscarriage"="#D73027", "Infection"="#FDB462", "Preterm"="#E41A1C"
)

# ── Theme ──
tpub <- theme_bw(base_size=11) +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
        axis.text=element_blank(), axis.ticks=element_blank(),
        legend.position="right", legend.title=element_blank(),
        legend.text=element_text(size=8), legend.key=element_blank(),
        plot.title=element_text(size=13, face="bold"))

# Panel A: By Dataset
pA <- ggplot(seu@meta.data[idx,], aes(umap_1, umap_2, color=dataset)) +
  geom_point(size=0.12, alpha=0.8) + scale_color_manual(values=dataset_cols) +
  labs(title="By Dataset", x="UMAP 1", y="UMAP 2") + tpub +
  guides(color=guide_legend(override.aes=list(size=3), ncol=1))

# Panel B: By Disease Group  
pB <- ggplot(seu@meta.data[idx,], aes(umap_1, umap_2, color=disease_short)) +
  geom_point(size=0.12, alpha=0.8) + scale_color_manual(values=disease_cols) +
  labs(title="By Disease Group", x="UMAP 1", y="UMAP 2") + tpub +
  guides(color=guide_legend(override.aes=list(size=3), ncol=1))

ggsave(file.path(OUT, "补充图FigB_dataset_disease.png"), pA+pB, width=18, height=7, dpi=300, bg="white")
cat("Saved: 补充图FigB_dataset_disease.png\n")
