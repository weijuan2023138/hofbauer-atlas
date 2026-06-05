#!/usr/bin/env Rscript
# Fig4: TGF-β superfamily chord — all TGFb/BMP/GDF/ACTIVIN interactions
library(CellChat); library(dplyr)

cellchat <- readRDS("results/cellchat_ucsf_mid.rds")

# Combine TGFb + BMP + GDF pathways
tgfb_all <- bind_rows(
  subsetCommunication(cellchat, signaling="TGFb"),
  subsetCommunication(cellchat, signaling="BMP"),
  subsetCommunication(cellchat, signaling="GDF")
)
cat(sprintf("Total TGF-β superfamily interactions: %d\n", nrow(tgfb_all)))
cat("Sources:", unique(tgfb_all$source), "\n")
cat("Targets:", unique(tgfb_all$target), "\n")

# Show all
for(i in 1:nrow(tgfb_all)) {
  cat(sprintf("  %s -> %s: %s -> %s [%s] %.2e\n",
    tgfb_all$source[i], tgfb_all$target[i],
    tgfb_all$ligand[i], tgfb_all$receptor[i],
    tgfb_all$pathway_name[i], tgfb_all$prob[i]))
}

# Generate chord diagram with all TGFb/BMP/GDF
png("figures/Fig4/Fig4_cellchat_TGFb_all.png", w=10, h=10, units="in", res=300)
# Use netVisual_chord_gene but for multiple pathways we need to aggregate
# Instead, create a custom network from the combined data
netVisual_chord_cell(cellchat, signaling=c("TGFb","BMP","GDF"),
  title.name="TGF-β superfamily signaling")
dev.off()

cat("Saved Fig4_cellchat_TGFb_all.png\n")
