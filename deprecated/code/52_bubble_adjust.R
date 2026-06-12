#!/usr/bin/env Rscript
# Bubble plots — narrower, bold font
library(CellChat); library(ggplot2)

cellchat <- readRDS("results/cellchat_ucsf_mid.rds")

# Incoming: who talks to HB
png("figures/Fig4/Fig4A_cellchat_bubble_incoming.png", w=6, h=7, units="in", res=300)
netVisual_bubble(cellchat,
  sources.use=c("FB","fEC","vEC","VCT","SCT","dNK","CD14_M"),
  targets.use="HB", remove.isolate=FALSE,
  title.name="Signaling to Hofbauer cells",
  color.text=TRUE,
  font.size=10, font.size.title=14)
dev.off()

# Outgoing: who HB talks to
png("figures/Fig4/Fig4E_cellchat_bubble_outgoing.png", w=6, h=7, units="in", res=300)
netVisual_bubble(cellchat,
  sources.use="HB",
  targets.use=c("FB","fEC","vEC","VCT","SCT","dNK","CD14_M"),
  remove.isolate=FALSE,
  title.name="Hofbauer signaling to neighbors",
  color.text=TRUE,
  font.size=10, font.size.title=14)
dev.off()

cat("Done\n")
