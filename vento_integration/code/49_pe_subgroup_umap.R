#!/usr/bin/env Rscript
# PE subgroup: UMAP 2×3 + proportion bar chart
library(Seurat); library(ggplot2); library(patchwork); library(cowplot); library(dplyr)
set.seed(42)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
detail <- labels$disease_detail[1:ncol(seu)]

seu$pe_group <- NA_character_
seu$pe_group[seu$dataset=="GSE290578" & detail=="PE"] <- "Early PE\n(GW29-34)"
seu$pe_group[seu$dataset=="GSE298602" & detail=="PreE_SF"] <- "Late Severe\n(GW37-40)"
seu$pe_group[seu$dataset=="GSE298602" & detail=="gHTN"] <- "Late Mild\n(gHTN)"
seu$pe_group[seu$dataset=="GSE298602" & detail=="Control"] <- "GSE298602\nControl"
seu$pe_group[detail %in% c("GSE173193","GSE298119")] <- "Late PE\n(term)"
seu$pe_group[detail=="Normal"] <- "Normal\nLate"

target <- c("Early PE\n(GW29-34)","Late Severe\n(GW37-40)","Late Mild\n(gHTN)",
            "Late PE\n(term)","Normal\nLate","GSE298602\nControl")
seu <- subset(seu, pe_group %in% target)
seu$pe_group <- factor(seu$pe_group, levels=target)

st_cols <- c("Pro-inflammatory"="#C62828","MHCII+ Antigen-presenting"="#E65100",
  "Homeostatic"="#1565C0","PRKN+ Autophagy"="#6A1B9A",
  "Vascular remodeling"="#2E7D32","MKI67+ Proliferating"="#455A64")

umap_coords <- Embeddings(seu, "umap")
xlim <- range(umap_coords[,1]) + c(-0.5,0.5); ylim <- range(umap_coords[,2]) + c(-0.5,0.5)

# ── UMAP 2×3 ──
# Row1: PE groups (early, late severe, late mild)
# Row2: controls (late PE term, normal late, GSE298602 control)
row1 <- c("Early PE\n(GW29-34)","Late Severe\n(GW37-40)","Late Mild\n(gHTN)")
row2 <- c("Late PE\n(term)","Normal\nLate","GSE298602\nControl")

plots <- list()
for(d in c(row1, row2)) {
  mask <- seu$pe_group == d
  df <- data.frame(x=umap_coords[mask,1], y=umap_coords[mask,2], subtype=seu$subtype_pred[mask])
  df <- df[sample(nrow(df)),]
  p <- ggplot(df, aes(x=x, y=y, color=subtype)) + geom_point(size=0.15, alpha=0.7) +
    scale_color_manual(values=st_cols) + coord_cartesian(xlim=xlim, ylim=ylim) +
    labs(title=d) + theme_bw(base_size=9) +
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
title <- ggdraw() + draw_label("Hofbauer subtypes across PE subgroups", fontface="bold", size=13)
fig <- plot_grid(title, fig, ncol=1, rel_heights=c(0.05,0.95))
ggsave(file.path(FIGDIR,"PE_subgroup_UMAP.png"), fig, w=12, h=7, dpi=300, bg="white")

# ── Proportion bar chart ──
prop <- prop.table(table(seu$subtype_pred, seu$pe_group), margin=2)*100
prop_df <- as.data.frame(prop); colnames(prop_df) <- c("Subtype","Group","Proportion")

pbar <- ggplot(prop_df, aes(x=Group, y=Proportion, fill=Subtype)) +
  geom_col(width=0.65) + scale_fill_manual(values=st_cols) +
  labs(title="Subtype proportions: PE subgroups", y="Proportion (%)", x="") +
  theme_bw(11) + theme(panel.grid=element_blank(),
    plot.title=element_text(face="bold",size=13,hjust=0.5,color="black"),
    axis.text.x=element_text(angle=35, hjust=1, size=9, color="black"),
    axis.text.y=element_text(size=9, color="black"), axis.title.y=element_text(size=10,color="black"),
    legend.position="right", legend.text=element_text(size=9), legend.title=element_text(size=10,face="bold"))
ggsave(file.path(FIGDIR,"PE_subgroup_proportions.png"), pbar, w=8, h=5, dpi=300, bg="white")

cat("Saved: PE_subgroup_UMAP.png + PE_subgroup_proportions.png\n")

# Print table
print(round(prop, 1))
