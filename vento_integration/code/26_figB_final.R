#!/usr/bin/env Rscript
# FigB: Dataset + Disease UMAP with per-cell labels and Fig1A colors
library(Seurat); library(ggplot2); library(patchwork)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
seu$disease_final <- factor(labels$disease_final[1:ncol(seu)],
  levels=c("Normal","PE","Miscarriage","Infection","Preterm"))

seu$UMAP_1 <- Embeddings(seu, "umap")[,1]
seu$UMAP_2 <- Embeddings(seu, "umap")[,2]
set.seed(42); idx <- sample(ncol(seu), 20000)

# Fig1A colors
disease_cols <- c("Normal"="#4575B4","PE"="#FC8D59","Miscarriage"="#D73027",
                  "Infection"="#FDB462","Preterm"="#E41A1C")

dataset_cols <- c(
  "E-MTAB-12421"="#D73027","GSE290578"="#4575B4","GSE214607"="#4DAF4A",
  "E-MTAB-12795"="#984EA3","GSE173193"="#FF7F00","GSE298119"="#F781BF",
  "GSE333257"="#878787","UCSF Li 2026"="#66C2A5",
  "E-MTAB-6701"="#E6AB02","GSE298602"="#A6761D"
)

tpub <- theme_bw(base_size=11) +
  theme(panel.grid=element_line(color="grey92",linewidth=0.2),
        panel.border=element_rect(color="black",fill=NA,linewidth=0.5),
        axis.text=element_blank(), axis.ticks=element_blank(),
        legend.position="right", legend.title=element_blank(),
        legend.text=element_text(size=9), legend.key=element_blank(),
        plot.title=element_text(size=13,face="bold"))

pA <- ggplot(seu@meta.data[idx,], aes(UMAP_1, UMAP_2, color=dataset)) +
  geom_point(size=0.12, alpha=0.8) + scale_color_manual(values=dataset_cols) +
  labs(title="By Dataset", x="UMAP 1", y="UMAP 2") + tpub +
  guides(color=guide_legend(override.aes=list(size=3), ncol=1))

pB <- ggplot(seu@meta.data[idx,], aes(UMAP_1, UMAP_2, color=disease_final)) +
  geom_point(size=0.12, alpha=0.8) + scale_color_manual(values=disease_cols) +
  labs(title="By Disease Group", x="UMAP 1", y="UMAP 2") + tpub +
  guides(color=guide_legend(override.aes=list(size=3), ncol=1))

ggsave(file.path(FIGDIR,"FigB_dataset_disease_UMAP.png"), pA+pB, width=18, height=7, dpi=300, bg="white")
cat("Saved: FigB_dataset_disease_UMAP.png\n")
