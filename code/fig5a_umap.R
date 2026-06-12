#!/usr/bin/env Rscript
# Fig5a: UMAP 2x3 — Hofbauer subtypes across 6 disease groups
library(Seurat); library(ggplot2); library(patchwork); library(cowplot); library(dplyr)

seu <- readRDS("results/Hofbauer_Atlas_Final.rds")
mask_gse <- seu$dataset == "GSE290578"; bcs <- colnames(seu)
is_norm_gse <- grepl("_Norm_", bcs) & mask_gse; is_pt_gse <- grepl("_Pt_", bcs) & mask_gse
seu$disease_clean <- NA
seu$disease_clean[seu$disease == "Normal_1st"] <- "Normal_Early"
seu$disease_clean[seu$disease == "Normal"] <- "Normal_Late"
seu$disease_clean[seu$disease == "RM/NC"] <- "Miscarriage"
seu$disease_clean[seu$disease == "Normal/Listeria/Toxoplasma/Malaria"] <- "Infection"
seu$disease_clean[is_norm_gse] <- "Normal_Late"; seu$disease_clean[is_pt_gse] <- "PE"
seu$disease_clean[seu$disease %in% c("PTL","PTNL")] <- "Preterm"
seu$disease_clean[seu$disease == "TL"] <- "Normal_Late"
keep <- seu$disease_clean %in% c("Normal_Early","Normal_Late","Miscarriage","PE","Infection","Preterm")
seu_sub <- subset(seu, cells=colnames(seu)[keep])

st_cols <- c("Pro-inflammatory"="#C62828","MHCII+ Antigen-presenting"="#E65100",
             "Homeostatic"="#1565C0","PRKN+ Autophagy"="#6A1B9A",
             "Vascular remodeling"="#2E7D32","MKI67+ Proliferating"="#455A64")

row1 <- c("Normal_Early","Miscarriage","Infection")
row2 <- c("Normal_Late","PE","Preterm")
labels <- c("Normal_Early"="Normal (Early)  GW4.5-10","Miscarriage"="Miscarriage  GW6-12",
            "Infection"="Infection  GW11-24","Normal_Late"="Normal (Late)  GW28-38",
            "PE"="PE  GW30-38","Preterm"="Preterm  GW30-36")
umap_coords <- Embeddings(seu_sub, "umap")
xlim <- range(umap_coords[,1]) + c(-0.5,0.5); ylim <- range(umap_coords[,2]) + c(-0.5,0.5)

plots <- list()
for(d in c(row1, row2)) {
  mask <- seu_sub$disease_clean == d
  df <- data.frame(x=umap_coords[mask,1], y=umap_coords[mask,2], subtype=seu_sub$subtype[mask])
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

ggsave("figures/Fig5/Fig5a_UMAP_disease.png", fig, w=12, h=7, dpi=300, bg="white")
cat("Fig5a done\n")
