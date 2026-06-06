#!/usr/bin/env Rscript
# Chord diagrams — FB/fEC/vEC → HB with Fig3 spatial colors
library(CellChat); library(ComplexHeatmap); library(grid)

cellchat <- readRDS("results/cellchat_ucsf_mid.rds")

# Match Fig3 spatial plot colors (32_spatial_plot.py)
all_cols <- c(
  "CD14_M"="#66C2A5","CD16_M"="#DA9870","dNK"="#BE979C","EVT"="#AB98C8",
  "FB"="#99D594","fEC"="#FC8D59","HB"="#D73027","SCT"="#4575B4",
  "VCT"="#91BFDB","vEC"="#B3B3B3"
)

for(src in c("FB","fEC","vEC")) {
  fname <- sprintf("figures/Fig4/Fig4_cellchat_chord_HB_%s.png", src)
  png(fname, w=16, h=14, units="in", res=300)
  par(mar=c(1, 1, 3, 4))
  
  netVisual_chord_gene(cellchat, sources.use=src, targets.use="HB",
    title.name=NULL,
    lab.cex=1.2, small.gap=2, big.gap=15,
    color.use=all_cols,
    show.legend=FALSE)
  title(sprintf("%s → HB", src), cex.main=2.5)
  
  lgd <- Legend(
    at = c(src, "HB"),
    type = "grid",
    legend_gp = gpar(fill = c(all_cols[src], all_cols["HB"])),
    title = "Cell State",
    title_gp = gpar(fontsize = 18, fontface = "bold"),
    labels_gp = gpar(fontsize = 16),
    grid_width = unit(8, "mm"),
    grid_height = unit(8, "mm")
  )
  draw(lgd, x = unit(1, "npc") - unit(5, "mm"), 
       y = unit(0.55, "npc"), just = c("right", "center"))
  
  dev.off()
  cat(sprintf("Saved %s\n", fname))
}
cat("Done\n")
