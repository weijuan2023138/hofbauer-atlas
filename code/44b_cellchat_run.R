#!/usr/bin/env Rscript
# CellChat on UCSF 2nd trimester: HB-stromal L-R analysis
library(Matrix); library(CellChat); library(dplyr)

# ── Load data ──
cat("Loading data...\n")
counts <- readMM("results/ucsf_mid_chat_counts.mtx")
data_mat <- readMM("results/ucsf_mid_chat_data.mtx")
meta <- read.csv("results/ucsf_mid_chat_meta.csv", row.names=1)
genes <- read.csv("results/ucsf_mid_chat_genes.csv")$gene
barcodes <- read.csv("results/ucsf_mid_chat_barcodes.csv")$barcode

rownames(counts) <- barcodes; colnames(counts) <- genes
rownames(data_mat) <- barcodes; colnames(data_mat) <- genes

# CellChat expects genes x cells
data_input <- t(data_mat)

cat(sprintf("Loaded: %d cells, %d genes\n", ncol(data_input), nrow(data_input)))
cat("Cell types:\n")
print(table(meta$celltype_shortname))

# ── Run CellChat ──
cat("\nRunning CellChat...\n")
cellchat <- createCellChat(object=data_input, meta=meta, group.by="celltype_shortname")
cellchat@DB <- CellChatDB.human

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type="triMean", population.size=TRUE)
cellchat <- filterCommunication(cellchat, min.cells=10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, "results/cellchat_ucsf_mid.rds")
cat("CellChat saved\n")

# ── HB-centric analysis ──
# 1. HB outgoing
hb_out <- subsetCommunication(cellchat, sources.use="HB")
hb_out_top <- hb_out %>% group_by(target, pathway_name) %>%
  summarise(total=sum(prob), .groups="drop") %>%
  arrange(desc(total)) %>% head(20)
cat("\n=== HB → others (top pathways) ===\n")
print(hb_out_top, n=15)

# 2. Incoming to HB (FB/fEC)
hb_fb <- subsetCommunication(cellchat, sources.use="FB", targets.use="HB")
cat("\n=== FB → HB (top 15 L-R pairs) ===\n")
hb_fb %>% arrange(desc(prob)) %>% head(15) %>%
  select(ligand, receptor, pathway_name, prob) %>% print()

hb_fec <- subsetCommunication(cellchat, sources.use="fEC", targets.use="HB")
cat("\n=== fEC → HB (top 15 L-R pairs) ===\n")
hb_fec %>% arrange(desc(prob)) %>% head(15) %>%
  select(ligand, receptor, pathway_name, prob) %>% print()

hb_vec <- subsetCommunication(cellchat, sources.use="vEC", targets.use="HB")
cat("\n=== vEC → HB (top 10 L-R pairs) ===\n")
hb_vec %>% arrange(desc(prob)) %>% head(10) %>%
  select(ligand, receptor, pathway_name, prob) %>% print()

# 3. Key pathways
cat("\n=== Key HB-related pathways ===\n")
hb_pathways <- hb_out %>% group_by(pathway_name) %>%
  summarise(outgoing=sum(prob), .groups="drop") %>% arrange(desc(outgoing)) %>% head(10)
print(hb_pathways)

cat("\nCellChat complete!\n")
