#!/usr/bin/env Rscript
# Extract all pathway names from CellChat object for classification
library(CellChat)

cellchat <- readRDS("results/cellchat_ucsf_mid.rds")
db <- cellchat@DB$interaction

# Get unique pathway names
pathways <- unique(db$pathway_name)
pathways <- pathways[!is.na(pathways)]

# Write to file
writeLines(pathways, "results/pathway_list.txt")

cat("Extracted", length(pathways), "pathways to results/pathway_list.txt\n")