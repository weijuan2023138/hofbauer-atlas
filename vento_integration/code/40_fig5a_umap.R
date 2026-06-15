#!/usr/bin/env Rscript
# Fig5a: UMAP 2×3 — Hofbauer subtypes across 6 disease groups (new Atlas)
library(Seurat); library(ggplot2); library(patchwork); library(cowplot); library(dplyr)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
seu$disease_final <- factor(labels$disease_final[1:ncol(seu)], 
  levels=c("Normal","PE","Miscarriage","Infection","Preterm"))

# Assign trimester
tri_map <- c("E-MTAB-12421"="Early","E-MTAB-6701"="Early","E-MTAB-12795"="Early",
             "GSE214607"="Early","UCSF Li 2026"="Mid","GSE290578"="Late",
             "GSE333257"="Late","GSE298602"="Late","GSE298119"="Late","GSE173193"="Late")
seu$trimester <- setNames(tri_map[as.character(seu$dataset)], colnames(seu))

# Build disease_clean
seu$disease_clean <- as.character(seu$disease_final)
seu$disease_clean[seu$disease_final=="Normal" & seu$trimester=="Early"] <- "Normal_Early"
seu$disease_clean[seu$disease_final=="Normal" & seu$trimester=="Late"] <- "Normal_Late"
keep <- seu$disease_clean %in% c("Normal_Early","Normal_Late","Miscarriage","PE","Infection","Preterm")
seu <- subset(seu, cells=colnames(seu)[keep])
seu$disease_clean <- factor(seu$disease_clean, 
  levels=c("Normal_Early","Miscarriage","Infection","Normal_Late","PE","Preterm"))

st_cols <- c("Pro-inflammatory"="#C62828","MHCII+ Antigen-presenting"="#E65100",
             "Homeostatic"="#1565C0","PRKN+ Autophagy"="#6A1B9A",
             "Vascular remodeling"="#2E7D32","MKI67+ Proliferating"="#455A64")

labels <- c("Normal_Early"="Normal (Early)  GW4.5-12","Miscarriage"="Miscarriage  GW6-12",
            "Infection"="Infection  GW4-8.5","Normal_Late"="Normal (Late)  GW29-40",
            "PE"="PE  GW29-40","Preterm"="Preterm  GW32-39")

umap_coords <- Embeddings(seu, "umap")
xlim <- range(umap_coords[,1]) + c(-0.5,0.5); ylim <- range(umap_coords[,2]) + c(-0.5,0.5)

plots <- list()
for(d in levels(seu$disease_clean)) {
  mask <- seu$disease_clean == d
  df <- data.frame(x=umap_coords[mask,1], y=umap_coords[mask,2], subtype=seu$subtype_pred[mask])
  df <- df[sample(nrow(df)),]
  p <- ggplot(df, aes(x=x, y=y, color=subtype)) + geom_point(size=0.15, alpha=0.7) +
    scale_color_manual(values=st_cols) + coord_cartesian(xlim=xlim, ylim=ylim) +
    labs(title=labels[d]) + theme_bw(base_size=9) +
    theme(panel.grid.major=element_line(color="grey92", linewidth=0.1),
          panel.grid.minor=element_blank(), panel.border=element_blank(),
          axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(),
          plot.title=element_text(face="bold",size=11,hjust=0.5),
          plot.margin=margin(0,0,0,0), legend.position="none")
  plots[[d]] <- p
}

dummy_df <- data.frame(subtype=factor(names(st_cols),levels=names(st_cols)), x=1:6, y=1:6)
pleg <- ggplot(dummy_df, aes(x=x, y=y, color=subtype)) + geom_point(size=3.5) +
  scale_color_manual(values=st_cols) + theme_void() +
  theme(legend.position="right", legend.title=element_blank(), legend.text=element_text(size=8),
        legend.key.size=unit(0.4,"cm"), legend.margin=margin(0,8,0,0))
leg <- get_legend(pleg)
p_grid <- plot_grid(plotlist=plots, nrow=2, ncol=3)
fig <- plot_grid(p_grid, leg, ncol=2, rel_widths=c(0.86,0.14))
title <- ggdraw() + draw_label("Hofbauer cell subtypes across pregnancy conditions", fontface="bold", size=13)
fig <- plot_grid(title, fig, ncol=1, rel_heights=c(0.05,0.95))
ggsave(file.path(FIGDIR,"Fig5a_UMAP_disease.png"), fig, w=12, h=7, dpi=300, bg="white")
cat("Fig5a done\n")
