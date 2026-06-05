#!/usr/bin/env Rscript
# Chord diagrams — original working version, just bigger title
library(CellChat); library(dplyr)

cellchat <- readRDS("results/cellchat_ucsf_mid.rds")

for(src in c("FB","fEC","vEC")) {
  png(sprintf("figures/Fig4/Fig4_cellchat_chord_HB_%s.png", src),
    w=14, h=14, units="in", res=300)
  
  netVisual_chord_gene(cellchat, sources.use=src, targets.use="HB",
    title.name=NULL,
    lab.cex=1.2, small.gap=2, big.gap=15,
    legend.pos.x=5, legend.pos.y=200)

  title(sprintf("%s → HB", src), cex.main=2.5)
  dev.off()
  cat(sprintf("Saved %s\n", src))
}
cat("Done\n")
