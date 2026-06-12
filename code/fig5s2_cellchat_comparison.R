#!/usr/bin/env Rscript
# Fig5S2/S3: CellChat infection vs control comparison
library(CellChat)

cc_ctrl <- readRDS("results/cellchat_hoo2024_control.rds")
cc_inf  <- readRDS("results/cellchat_hoo2024_infected.rds")
merged  <- mergeCellChat(list(Control=cc_ctrl, Infected=cc_inf), add.names=c("Control","Infected"))

# S2: LR pair communication probability
png("figures/Fig5/补充图Fig5B_LR_probability.png", w=2400, h=1800, res=300, bg="white")
netVisual_bubble(merged, comparison=c(1,2), angle=0,
  title.name="LR pair communication probability",
  font.size=10, font.size.title=14)
dev.off()

# S3: Pathway communication strength change
png("figures/Fig5/补充图Fig5C_pathway_change.png", w=2400, h=1500, res=300, bg="white")
rankNet(merged, mode="comparison", stacked=F, do.stat=T, font.size=10)
dev.off()

cat("Fig5S2/S3 done\n")
