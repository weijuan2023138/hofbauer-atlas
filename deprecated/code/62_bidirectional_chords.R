#!/usr/bin/env Rscript
# Bidirectional chord diagrams: HB ↔ vEC / fEC / FB
library(CellChat); library(ComplexHeatmap); library(grid)

cellchat <- readRDS("results/cellchat_ucsf_mid.rds")
all_cols <- c(
  "CD14_M"="#66C2A5","dNK"="#BE979C","FB"="#99D594",
  "fEC"="#FC8D59","HB"="#D73027","SCT"="#4575B4",
  "VCT"="#91BFDB","vEC"="#B3B3B3"
)

bidirectional_chord <- function(cell_a, cell_b, fname, title_text) {
  png(fname, w=16, h=14, units="in", res=300)
  par(mar=c(1, 1, 3, 4))
  
  netVisual_chord_gene(cellchat,
    sources.use=c(cell_a, cell_b),
    targets.use=c(cell_b, cell_a),
    title.name=title_text,
    lab.cex=1.2, small.gap=2, big.gap=15,
    color.use=all_cols, show.legend=FALSE)
  
  lgd <- Legend(at=c(cell_a, cell_b), type="grid",
    legend_gp=gpar(fill=c(all_cols[cell_a], all_cols[cell_b])),
    title="Cell State", title_gp=gpar(fontsize=18, fontface="bold"),
    labels_gp=gpar(fontsize=16),
    grid_width=unit(8,"mm"), grid_height=unit(8,"mm"))
  draw(lgd, x=unit(1,"npc")-unit(5,"mm"), y=unit(0.55,"npc"),
       just=c("right","center"))
  dev.off()
  cat(sprintf("Saved: %s\n", fname))
}

bidirectional_chord("HB", "vEC", "figures/Fig4/补充图Fig4C_cellchat_chord_HB_vEC_bi.png",
                    "HB \u2194 vEC")
bidirectional_chord("HB", "fEC", "figures/Fig4/补充图Fig4C_cellchat_chord_HB_fEC_bi.png",
                    "HB \u2194 fEC")
bidirectional_chord("HB", "FB",  "figures/Fig4/补充图Fig4C_cellchat_chord_HB_FB_bi.png",
                    "HB \u2194 FB")

cat("\nDone: 3 bidirectional chords\n")
