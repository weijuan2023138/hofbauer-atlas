#!/usr/bin/env Rscript
# Replicate Fig1 (A/B/C) for new 10-dataset labeled Atlas — identical format
library(Seurat); library(ggplot2); library(patchwork); library(dplyr)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"
dir.create(FIGDIR, showWarnings=FALSE, recursive=TRUE)

seu <- readRDS(INPUT)
# Use UMAP from embeddings, add to meta.data
seu$UMAP_1 <- Embeddings(seu, "umap")[,1]
seu$UMAP_2 <- Embeddings(seu, "umap")[,2]
# Use label transfer subtype
seu$subtype <- seu$subtype_pred

cat(ncol(seu),"cells,",length(unique(seu$subtype)),"subtypes\n")

# ── Colors (exact from original) ──
subtype_colors <- c(
  "Pro-inflammatory"="#C62828", "MHCII+ Antigen-presenting"="#E65100",
  "Homeostatic"="#1565C0", "PRKN+ Autophagy"="#6A1B9A",
  "Vascular remodeling"="#2E7D32", "MKI67+ Proliferating"="#455A64"
)

# Dataset colors — update for new names
dataset_colors <- c(
  "E-MTAB-12421"="#D73027", "GSE290578"="#4575B4", "GSE214607"="#4DAF4A",
  "E-MTAB-12795"="#984EA3", "GSE173193"="#FF7F00", "GSE298119"="#F781BF",
  "GSE333257"="#878787", "UCSF Li 2026"="#66C2A5",
  "E-MTAB-6701"="#E6AB02", "GSE298602"="#A6761D"
)

disease_colors <- c(
  "Normal 1st"="#4DBBD5", "Normal"="#00A087",
  "Normal+PE"="#7E6148", "PE"="#C62828", "PE/Control"="#C62828",
  "Miscarriage"="#F39B7F", "Infection"="#DC0000", "Preterm"="#E18727"
)

# ── Theme (exact from original) ──
tpub <- theme_bw(base_size=11) +
  theme(panel.grid=element_line(color="grey92",linewidth=0.2),
        panel.border=element_rect(color="black",fill=NA,linewidth=0.5),
        axis.text=element_blank(), axis.ticks=element_blank(),
        legend.position="right", legend.title=element_blank(),
        legend.text=element_text(size=9), legend.key=element_blank(),
        plot.title=element_text(size=13,face="bold"))

set.seed(42)
idx <- sample(ncol(seu))

# ═══ FIGURE A: Subtype UMAP ═══
p1 <- ggplot(seu@meta.data[idx,], aes(x=UMAP_1, y=UMAP_2, color=subtype)) +
  geom_point(size=0.15, alpha=0.8) +
  scale_color_manual(values=subtype_colors) +
  labs(title="Hofbauer Atlas — 6 Subtypes", x="UMAP 1", y="UMAP 2") +
  tpub + guides(color=guide_legend(override.aes=list(size=3)))
ggsave(file.path(FIGDIR,"FigA_subtype_UMAP.png"), p1, width=9, height=7, dpi=300, bg="white")
cat("Saved FigA_subtype\n")

# ═══ FIGURE B: Dataset + Disease UMAP ═══
# Disease mapping
disease_map <- c(
  "E-MTAB-12421"="Normal 1st", "E-MTAB-6701"="Normal 1st",
  "UCSF Li 2026"="Normal", "GSE290578"="Normal+PE",
  "GSE214607"="Miscarriage", "E-MTAB-12795"="Infection",
  "GSE173193"="PE", "GSE298119"="PE", "GSE298602"="PE/Control",
  "GSE333257"="Preterm"
)
disease_levels <- c("Normal 1st","Normal","Normal+PE","PE","PE/Control",
                    "Miscarriage","Infection","Preterm")

ds_vec <- setNames(disease_map[as.character(seu$dataset)], colnames(seu))
seu <- AddMetaData(seu, metadata=data.frame(
  disease_short=factor(ds_vec, levels=disease_levels),
  row.names=colnames(seu)
))

p2a <- ggplot(seu@meta.data[idx,], aes(UMAP_1, UMAP_2, color=dataset)) +
  geom_point(size=0.12, alpha=0.8) + scale_color_manual(values=dataset_colors) +
  labs(title="By Dataset", x="UMAP 1", y="UMAP 2") +
  tpub + guides(color=guide_legend(override.aes=list(size=3), ncol=1))

p2b <- ggplot(seu@meta.data[idx,], aes(UMAP_1, UMAP_2, color=disease_short)) +
  geom_point(size=0.12, alpha=0.8) + scale_color_manual(values=disease_colors) +
  labs(title="By Disease Group", x="UMAP 1", y="UMAP 2") +
  tpub + guides(color=guide_legend(override.aes=list(size=3), ncol=1))

ggsave(file.path(FIGDIR,"FigB_dataset_disease_UMAP.png"), p2a+p2b, width=18, height=7, dpi=300, bg="white")
cat("Saved FigB_dataset_disease\n")

# ═══ FIGURE C: DotPlot ═══
DefaultAssay(seu) <- "RNA"
Idents(seu) <- "subtype"

classifier_markers <- c(
  "LYVE1","CD36","DAB2","VSIG4","F13A1","SPP1",
  "FSCN1","COLEC12","ABCG2","CD5L",
  "HLA-DRB5","HLA-DRB1","HLA-DPA1","CST7","HLA-DRA",
  "FOLR2","CD163","MRC1","TREM2","C1QA","C1QB",
  "IL1B","TNF","CXCL8","CCL13","PRKN","C9",
  "MKI67","BUB1B","FN1","PAPPA","FLT1"
)
present <- classifier_markers[classifier_markers %in% rownames(seu)]
cat("Markers present:", length(present), "/", length(classifier_markers), "\n")

p3 <- DotPlot(seu, features=present, assay="RNA") + RotatedAxis() +
  scale_color_gradientn(colors=c("lightgrey","#1B6B93","#C62828")) +
  labs(title="Key Marker Genes by Subtype") + theme_bw(11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
        panel.grid=element_line(color="grey92",linewidth=0.2), legend.position="right")
ggsave(file.path(FIGDIR,"FigC_dotplot_markers.png"), p3, width=18, height=6, dpi=300, bg="white")
cat("Saved FigC_dotplot\n")

cat("\n========================================\n")
cat("ALL FIGURES in", FIGDIR, "\n")
cat("  FigA_subtype_UMAP.png\n")
cat("  FigB_dataset_disease_UMAP.png\n")
cat("  FigC_dotplot_markers.png\n")
cat("========================================\n")
