#!/usr/bin/env Rscript
# Final Atlas: UMAP visualization + marker export + beautified figures
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

INPUT <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/Hofbauer_Atlas_Final.rds'
OUTDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results'
FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'
dir.create(FIGDIR, showWarnings=FALSE)

seu <- readRDS(INPUT)
cat("Loaded:", ncol(seu), "cells\n")

# ============================================================
# Color scheme
# ============================================================
subtype_colors <- c(
  "Pro-inflammatory" = "#C62828",
  "MHCII+/CCL13+C1Q+" = "#BF4E1A",
  "Homeostatic/SPP1+" = "#1B6B93",
  "PRKN+ Autophagy" = "#7B4FA0",
  "Vascular remodeling" = "#2D8B57",
  "MKI67+ Proliferating" = "#37474F"
)

dataset_colors <- c(
  "D1" = "#E41A1C", "D2" = "#377EB8", "D3" = "#4DAF4A",
  "D4" = "#984EA3", "D5" = "#FF7F00", "D6" = "#A65628",
  "D7" = "#F781BF", "D8" = "#999999", "D9" = "#66C2A5"
)

# ============================================================
# Rename clusters
# ============================================================
cluster_names <- c(
  "0" = "Pro-inflammatory",
  "1" = "MHCII+/CCL13+C1Q+",
  "2" = "Homeostatic/SPP1+",
  "3" = "PRKN+ Autophagy",
  "4" = "Vascular remodeling",
  "5" = "MKI67+ Proliferating"
)
subtype_vec <- cluster_names[as.character(seu$cluster_final)]
names(subtype_vec) <- colnames(seu)
seu$subtype <- factor(subtype_vec, levels=unname(cluster_names))

# ============================================================
# Publication theme
# ============================================================
theme_pub <- theme_minimal(base_size=11) +
  theme(
    panel.grid=element_blank(),
    axis.line=element_line(color="black", linewidth=0.3),
    axis.ticks=element_line(color="black", linewidth=0.3),
    axis.text=element_blank(),
    axis.title=element_text(size=12),
    legend.position="bottom",
    legend.title=element_blank(),
    legend.text=element_text(size=9),
    legend.key.size=unit(3,"mm"),
    plot.title=element_text(size=13, face="bold"),
    plot.margin=margin(5,5,5,5)
  )

# ============================================================
# UMAP data
# ============================================================
umap <- Embeddings(seu, "umap")
seu$UMAP_1 <- umap[,1]
seu$UMAP_2 <- umap[,2]

# ============================================================
# FIGURE 1: UMAP by subtype
# ============================================================
set.seed(42)
idx <- sample(ncol(seu))

p1 <- ggplot(seu@meta.data[idx,], aes(x=UMAP_1, y=UMAP_2, color=subtype)) +
  geom_point(size=0.3, alpha=0.7) +
  scale_color_manual(values=subtype_colors) +
  labs(title="Hofbauer Atlas — 6 Subtypes (n=17,896)", x="UMAP 1", y="UMAP 2") +
  theme_pub +
  guides(color=guide_legend(nrow=2, override.aes=list(size=2)))

ggsave(file.path(FIGDIR, "Fig1_UMAP_subtype.png"), p1, width=8, height=7, dpi=300, bg="white")
ggsave(file.path(FIGDIR, "Fig1_UMAP_subtype.pdf"), p1, width=8, height=7, bg="white")
cat("Saved Fig1\n")

# ============================================================
# FIGURE 2: UMAP by dataset (numbered)
# ============================================================
dataset_map <- c("Arutyunyan"="D1","GSE290578"="D2","gse214607"="D3",
                 "hoo_2024"="D4","gse173193"="D5","gse183338"="D6",
                 "gse298119"="D7","my_preterm_cohort"="D8","UCSF_Li_2026"="D9")
dl <- dataset_map[seu$dataset]; names(dl) <- colnames(seu)
seu$dataset_label <- factor(dl, levels=paste0("D",1:9))

p2 <- ggplot(seu@meta.data[idx,], aes(x=UMAP_1, y=UMAP_2, color=dataset_label)) +
  geom_point(size=0.3, alpha=0.7) +
  scale_color_manual(values=dataset_colors) +
  labs(title="Hofbauer Atlas — by Dataset", x="UMAP 1", y="UMAP 2") +
  theme_pub +
  guides(color=guide_legend(nrow=2, override.aes=list(size=2)))
ggsave(file.path(FIGDIR, "Fig2_UMAP_dataset.png"), p2, width=8, height=7, dpi=300, bg="white")
ggsave(file.path(FIGDIR, "Fig2_UMAP_dataset.pdf"), p2, width=8, height=7, bg="white")
cat("Saved Fig2\n")

# Split by dataset
p2s <- ggplot(seu@meta.data[idx,], aes(x=UMAP_1, y=UMAP_2, color=subtype)) +
  geom_point(size=0.1, alpha=0.6) +
  scale_color_manual(values=subtype_colors, guide="none") +
  facet_wrap(~dataset_label, ncol=3) +
  labs(x="UMAP 1", y="UMAP 2") +
  theme_pub + theme(strip.text=element_text(size=10, face="bold"))
ggsave(file.path(FIGDIR, "Fig2_split_dataset.png"), p2s, width=14, height=10, dpi=300, bg="white")
cat("Saved Fig2_split\n")

# ============================================================
# FIGURE 3: Barplot by disease
# ============================================================
short <- c("Normal 1st trimester"="Normal\n1st","Normal 1st/2nd/Term"="Normal\n1st/2nd",
           "Normal 3rd trimester / Preeclampsia"="Late\npreg","Preeclampsia"="PE",
           "Preterm No Labor"="PTNL","Preterm Labor"="PTL",
           "Term Labor"="TL","Miscarriage / Normal"="Miscarr.","Infection"="Infection")
prop <- seu@meta.data %>% count(disease_group, subtype) %>%
  group_by(disease_group) %>% mutate(pct=n/sum(n)*100, ds=short[disease_group])
prop$ds <- factor(prop$ds, levels=unname(short))

p3 <- ggplot(prop, aes(x=ds, y=pct, fill=subtype)) +
  geom_bar(stat="identity", width=0.7) +
  scale_fill_manual(values=subtype_colors) +
  labs(x="", y="Proportion (%)", title="Subtype Distribution by Disease") +
  theme_minimal(11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        legend.position="bottom", legend.title=element_blank(),
        legend.text=element_text(size=8), panel.grid.major.x=element_blank())
ggsave(file.path(FIGDIR, "Fig3_barplot_disease.png"), p3, width=12, height=5, dpi=300, bg="white")
cat("Saved Fig3\n")

# ============================================================
# FIGURE 4: Key markers on UMAP
# ============================================================
DefaultAssay(seu) <- "RNA"
markers <- c("IL1B","TNF","HLA-DRA","CCL13","FOLR2","CD163","PRKN","C9","MKI67","SPP1","DAB2")
pl <- list()
for(i in seq_along(markers)) {
  g <- markers[i]
  if(g %in% rownames(seu)) {
    pl[[i]] <- FeaturePlot(seu, features=g, pt.size=0.3, order=TRUE) +
      scale_color_gradientn(colors=c("grey90","#fee0d2","#de2d26")) +
      labs(title=g) + theme_pub + NoLegend()
  }
}
p4 <- wrap_plots(pl, ncol=4)
ggsave(file.path(FIGDIR, "Fig4_marker_UMAP.png"), p4, width=16, height=10, dpi=300, bg="white")
cat("Saved Fig4\n")

# ============================================================
# Export markers with names
# ============================================================
markers <- read.csv(file.path(OUTDIR, 'Hofbauer_Atlas_Final_all_markers.csv'))
markers$subtype_name <- cluster_names[as.character(markers$cluster)]
write.csv(markers, file.path(OUTDIR, 'Hofbauer_Atlas_Final_markers_named.csv'), row.names=FALSE)

top10 <- markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC) %>% arrange(cluster, desc(avg_log2FC))
top10$subtype_name <- cluster_names[as.character(top10$cluster)]
write.csv(top10, file.path(OUTDIR, 'Hofbauer_Atlas_Final_top10_named.csv'), row.names=FALSE)
cat("Saved markers\n")

# Save final RDS
saveRDS(seu, file.path(OUTDIR, 'Hofbauer_Atlas_Final.rds'))
cat(sprintf("\nFinal RDS: %d cells, %d subtypes\n", ncol(seu), length(unique(seu$subtype))))

cat("\n========================================\n")
cat("ALL FILES:\n")
cat("  Hofbauer_Atlas_Final.rds\n")
cat("  Hofbauer_Atlas_Final_markers_named.csv\n")
cat("  Hofbauer_Atlas_Final_top10_named.csv\n")
cat("  Fig1_UMAP_subtype.png/pdf\n")
cat("  Fig2_UMAP_dataset.png/pdf\n")
cat("  Fig2_split_dataset.png\n")
cat("  Fig3_barplot_disease.png/pdf\n")
cat("  Fig4_marker_UMAP.png\n")
cat("========================================\n")
