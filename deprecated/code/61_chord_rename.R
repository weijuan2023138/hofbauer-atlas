#!/usr/bin/env Rscript
# Regenerate chords with clean naming (no title, chord itself labels cell types)
library(CellChat); library(ComplexHeatmap); library(grid)

cellchat <- readRDS("results/cellchat_ucsf_mid.rds")

all_cols <- c(
  "CD14_M"="#66C2A5","dNK"="#BE979C","FB"="#99D594",
  "fEC"="#FC8D59","HB"="#D73027","SCT"="#4575B4",
  "VCT"="#91BFDB","vEC"="#B3B3B3"
)

gen_chord <- function(src, tgt, fname, title_text) {
  png(fname, w=16, h=14, units="in", res=300)
  par(mar=c(1, 1, 3, 4))
  netVisual_chord_gene(cellchat, sources.use=src, targets.use=tgt,
    title.name=NULL,
    lab.cex=1.2, small.gap=2, big.gap=15,
    color.use=all_cols, show.legend=FALSE)
  title(title_text, cex.main=2.5)
  lgd <- Legend(at=c(src, tgt), type="grid",
    legend_gp=gpar(fill=c(all_cols[src], all_cols[tgt])),
    title="Cell State", title_gp=gpar(fontsize=18, fontface="bold"),
    labels_gp=gpar(fontsize=16),
    grid_width=unit(8,"mm"), grid_height=unit(8,"mm"))
  draw(lgd, x=unit(1,"npc")-unit(5,"mm"), y=unit(0.55,"npc"),
       just=c("right","center"))
  dev.off()
  cat(sprintf("Saved %s\n", fname))
}

# Incoming
gen_chord("vEC", "HB", "figures/Fig4/Fig4_cellchat_chord_vEC_to_HB.png",
          "vEC \u2192 HB")
gen_chord("FB",  "HB", "figures/Fig4/Fig4_cellchat_chord_FB_to_HB.png",
          "FB \u2192 HB")
gen_chord("fEC", "HB", "figures/Fig4/Fig4_cellchat_chord_fEC_to_HB.png",
          "fEC \u2192 HB")

# Outgoing
gen_chord("HB", "vEC",    "figures/Fig4/Fig4_cellchat_chord_HB_to_vEC.png",
          "HB \u2192 vEC")
gen_chord("HB", "CD14_M", "figures/Fig4/Fig4_cellchat_chord_HB_to_CD14_M.png",
          "HB \u2192 CD14_M")
gen_chord("HB", "dNK",    "figures/Fig4/Fig4_cellchat_chord_HB_to_dNK.png",
          "HB \u2192 dNK")

cat("\nDone. 3 incoming + 3 outgoing chords with clean names.\n")
