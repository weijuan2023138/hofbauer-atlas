#!/usr/bin/env Rscript
# Final comprehensive figures: classifier markers, UMAPs, DotPlot
library(Seurat); library(ggplot2); library(patchwork); library(dplyr)

seu <- readRDS("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/Hofbauer_Atlas_Final.rds")
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures"

cat(ncol(seu),"cells,",length(unique(seu$cluster_final)),"clusters\n")

# ============================================================
# Better color scheme
# ============================================================
subtype_colors <- c(
  "Pro-inflammatory" = "#C62828",
  "MHCII+/CCL13+C1Q+" = "#E65100",
  "Homeostatic/SPP1+" = "#1565C0",
  "PRKN+ Autophagy" = "#6A1B9A",
  "Vascular remodeling" = "#2E7D32",
  "MKI67+ Proliferating" = "#455A64"
)

dataset_colors <- c(
  "Arutyunyan" = "#D73027",
  "GSE290578" = "#4575B4",
  "gse214607" = "#4DAF4A",
  "hoo_2024" = "#984EA3",
  "gse173193" = "#FF7F00",
  "gse183338" = "#A65628",
  "gse298119" = "#F781BF",
  "my_preterm_cohort" = "#878787",
  "UCSF_Li_2026" = "#66C2A5"
)

disease_colors <- c(
  "Normal 1st" = "#4DBBD5",
  "Normal 1st/2nd" = "#00A087",
  "Late pregnancy" = "#7E6148",
  "PE" = "#C62828",
  "PTNL" = "#BCAAA4",
  "PTL" = "#E18727",
  "TL" = "#3C5488",
  "Miscarriage" = "#F39B7F",
  "Infection" = "#DC0000"
)

# ============================================================
# Theme
# ============================================================
tpub <- theme_bw(base_size=11) +
  theme(panel.grid=element_line(color="grey92",linewidth=0.2),
        panel.border=element_rect(color="black",fill=NA,linewidth=0.5),
        axis.text=element_blank(), axis.ticks=element_blank(),
        legend.position="right", legend.title=element_blank(),
        legend.text=element_text(size=9), legend.key=element_blank(),
        plot.title=element_text(size=13,face="bold"))

# ============================================================
# FIGURE 1: Subtype UMAP
# ============================================================
set.seed(42)
idx <- sample(ncol(seu))

p1 <- ggplot(seu@meta.data[idx,], aes(x=UMAP_1, y=UMAP_2, color=subtype)) +
  geom_point(size=0.15, alpha=0.8) +
  scale_color_manual(values=subtype_colors) +
  labs(title="Hofbauer Atlas — 6 Subtypes", x="UMAP 1", y="UMAP 2") +
  tpub + guides(color=guide_legend(override.aes=list(size=3)))

ggsave(file.path(FIGDIR,"FigA_subtype_UMAP.png"), p1, width=9, height=7, dpi=300, bg="white")
ggsave(file.path(FIGDIR,"FigA_subtype_UMAP.pdf"), p1, width=9, height=7, bg="white")
cat("Saved FigA_subtype\n")

# ============================================================
# FIGURE 2: Dataset UMAP (original names)
# ============================================================
seu$dataset_short <- factor(seu$dataset,
  levels=c("Arutyunyan","GSE290578","gse214607","hoo_2024",
           "gse173193","gse183338","gse298119","my_preterm_cohort","UCSF_Li_2026"))

p2a <- ggplot(seu@meta.data[idx,], aes(UMAP_1, UMAP_2, color=dataset_short)) +
  geom_point(size=0.12, alpha=0.8) + scale_color_manual(values=dataset_colors) +
  labs(title="By Dataset", x="UMAP 1", y="UMAP 2") +
  tpub + guides(color=guide_legend(override.aes=list(size=3), ncol=1))

# Disease group UMAP
short <- c("Normal 1st trimester"="Normal 1st","Normal 1st/2nd/Term"="Normal 1st/2nd",
           "Normal 3rd trimester / Preeclampsia"="Late pregnancy","Preeclampsia"="PE",
           "Preterm No Labor"="PTNL","Preterm Labor"="PTL","Term Labor"="TL",
           "Miscarriage / Normal"="Miscarriage","Infection"="Infection")
seu$disease_short <- factor(short[seu$disease_group],
  levels=c("Normal 1st","Normal 1st/2nd","Late pregnancy","PE","PTNL","PTL","TL","Miscarriage","Infection"))

p2b <- ggplot(seu@meta.data[idx,], aes(UMAP_1, UMAP_2, color=disease_short)) +
  geom_point(size=0.12, alpha=0.8) + scale_color_manual(values=disease_colors) +
  labs(title="By Disease Group", x="UMAP 1", y="UMAP 2") +
  tpub + guides(color=guide_legend(override.aes=list(size=3), ncol=1))

ggsave(file.path(FIGDIR,"FigB_dataset_disease_UMAP.png"), p2a+p2b, width=18, height=7, dpi=300, bg="white")
ggsave(file.path(FIGDIR,"FigB_dataset_disease_UMAP.pdf"), p2a+p2b, width=18, height=7, bg="white")
cat("Saved FigB_dataset_disease\n")

# ============================================================
# FIGURE 3: DotPlot of key marker genes
# ============================================================
DefaultAssay(seu) <- "RNA"
Idents(seu) <- "subtype"

classifier_markers <- c(
  # HB-up (classifier top genes)
  "LYVE1","CD36","DAB2","VSIG4","F13A1","SPP1",
  "FSCN1","COLEC12","ABCG2","CD5L",
  # MAT-up
  "HLA-DRB5","HLA-DRB1","HLA-DPA1","CST7","HLA-DRA",
  # Hofbauer identity
  "FOLR2","CD163","MRC1","TREM2","C1QA","C1QB",
  # Functional
  "IL1B","TNF","CXCL8","CCL13","PRKN","C9",
  "MKI67","BUB1B","FN1","PAPPA","FLT1"
)

present <- classifier_markers[classifier_markers %in% rownames(seu)]
cat("Markers present:", length(present), "/", length(classifier_markers), "\n")

p3 <- DotPlot(seu, features=present, assay="RNA") +
  RotatedAxis() +
  scale_color_gradientn(colors=c("lightgrey","#1B6B93","#C62828")) +
  labs(title="Key Marker Genes by Subtype") +
  theme_bw(11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
        panel.grid=element_line(color="grey92",linewidth=0.2),
        legend.position="right")

ggsave(file.path(FIGDIR,"FigC_dotplot_markers.png"), p3, width=18, height=6, dpi=300, bg="white")
ggsave(file.path(FIGDIR,"FigC_dotplot_markers.pdf"), p3, width=18, height=6, bg="white")
cat("Saved FigC_dotplot\n")

# ============================================================
# FIGURE 4: Classifier genes table
# ============================================================
classifier_genes <- read.csv(
  "/home/weijuan/文档/胎盘单细胞数据/results/phase1_classifier/classifier_genes.csv")
write.csv(classifier_genes,
  file.path(FIGDIR,"../results/classifier_genes_table.csv"), row.names=FALSE)
cat("Saved classifier_genes_table.csv (", nrow(classifier_genes), "genes)\n")

# ============================================================
# Print summary
# ============================================================
cat("\n========================================\n")
cat("ALL FIGURES:\n")
cat("  FigA_subtype_UMAP.png/pdf\n")
cat("  FigB_dataset_disease_UMAP.png/pdf\n")
cat("  FigC_dotplot_markers.png/pdf\n")
cat("  classifier_genes_table.csv\n")
cat("========================================\n")
