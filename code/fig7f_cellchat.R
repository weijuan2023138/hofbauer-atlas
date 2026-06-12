#!/usr/bin/env Rscript
# Fig7f: Intra-HB CellChat — 6 subtype communication network
setwd("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration")
suppressMessages(library(Seurat)); library(CellChat); library(patchwork); library(ggplot2)

so <- readRDS("results/Hofbauer_Atlas_Final.rds")
DefaultAssay(so) <- "RNA"

subtypes <- c("Vascular remodeling","MHCII+ Antigen-presenting","Pro-inflammatory",
              "Homeostatic","PRKN+ Autophagy","MKI67+ Proliferating")
so$subtype <- factor(so$subtype, levels=subtypes)

# ===== CellChat: HB subtypes as groups =====
Idents(so) <- so$subtype
data_input <- GetAssayData(so, assay="RNA", layer="data")

# Downsample to 500 per subtype for speed
set.seed(42)
cells_use <- c()
for(st in subtypes) {
  cidx <- which(so$subtype == st)
  if(length(cidx) > 500) cidx <- sample(cidx, 500)
  cells_use <- c(cells_use, cidx)
}
so_sub <- so[, cells_use]
data_input <- GetAssayData(so_sub, assay="RNA", layer="data")

# Create CellChat directly from matrix (bypass Seurat API incompatibility)
data_matrix <- as.matrix(GetAssayData(so_sub, assay="RNA", layer="data"))
meta <- so_sub@meta.data

cc <- createCellChat(object=data_matrix, meta=meta, group.by="subtype")
cc@DB <- CellChatDB.human
cc <- subsetData(cc)
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)
cc <- projectData(cc, PPI.human)
cc <- computeCommunProb(cc, type="triMean", trim=0.1)
cc <- filterCommunication(cc, min.cells=10)
cc <- computeCommunProbPathway(cc)
cc <- aggregateNet(cc)
cc <- netAnalysis_computeCentrality(cc, slot.name="netP")

# ===== Panel 1: Interaction count/strength netVisual =====
dir.create("figures/Fig7", showWarnings=FALSE)
png("figures/Fig7/Fig7f_network.png", w=10, h=5, units="in", res=300, bg="white")
par(mfrow=c(1,2))
netVisual_circle(cc@net$count,  vertex.weight=as.numeric(table(cc@idents)),
  weight.scale=TRUE, label.edge=FALSE, title.name="Number of interactions")
netVisual_circle(cc@net$weight, vertex.weight=as.numeric(table(cc@idents)),
  weight.scale=TRUE, label.edge=FALSE, title.name="Interaction strength")
dev.off()

# ===== Panel 2: Heatmap of significant L-R pairs =====
png("figures/Fig7/Fig7f_heatmap.png", w=14, h=10, units="in", res=300, bg="white")
netAnalysis_signalingRole_heatmap(cc, pattern="outgoing", width=5, height=14, font.size=8)
dev.off()

# ===== Panel 3: Chord diagram of top pathways =====
pathways <- cc@netP$pathways
if(length(pathways) > 0) {
  top_path <- pathways[1:min(3, length(pathways))]
  for(pw in top_path) {
    png(sprintf("figures/Fig7/Fig7f_chord_%s.png", pw), w=8, h=8, units="in", res=300, bg="white")
    try(netVisual_chord_gene(cc, signaling=pw, slot.name="net", lab.cex=0.8, legend.pos.x=8, legend.pos.y=25))
    dev.off()
  }
}

# ===== Summary stats =====
cat("\n=== Subtype communication summary ===\n")
cat("Total interactions:", sum(cc@net$count), "\n")
cat("Total strength:", round(sum(cc@net$weight), 2), "\n")
for(st in subtypes) {
  out_count <- sum(cc@net$count[st,])
  in_count  <- sum(cc@net$count[,st])
  cat(sprintf("  %-30s  out=%d  in=%d\n", st, out_count, in_count))
}

# Save CellChat object
saveRDS(cc, "results/cellchat_hb_subtypes.rds")
message("Fig7f done")
