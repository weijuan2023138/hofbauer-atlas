#!/usr/bin/env Rscript
# Extract incoming and outgoing pathways for HB
library(CellChat)

cellchat <- readRDS("results/cellchat_ucsf_mid.rds")

# Incoming pathways (to HB)
inc_raw <- subsetCommunication(cellchat, targets.use="HB")
inc_pathways <- unique(inc_raw$pathway_name)

# Outgoing pathways (from HB)
out_raw <- subsetCommunication(cellchat, sources.use="HB")
out_pathways <- unique(out_raw$pathway_name)

# Write to files
writeLines(inc_pathways, "results/incoming_pathways.txt")
writeLines(out_pathways, "results/outgoing_pathways.txt")

cat("Incoming pathways:", length(inc_pathways), "\n")
cat("Outgoing pathways:", length(out_pathways), "\n")