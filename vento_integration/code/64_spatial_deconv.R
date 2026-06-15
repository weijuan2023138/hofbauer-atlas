#!/usr/bin/env Rscript
# Spatial subtype mapping via AddModuleScore (fast, no SPOTlight dependency)
library(Seurat); library(ggplot2); library(dplyr)
set.seed(42)

REF <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
SPDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/spatial"

ref <- readRDS(REF)
Idents(ref) <- "subtype_pred"

# Marker genes per subtype (top 20)
markers_list <- list()
for(st in levels(Idents(ref))) {
  deg <- FindMarkers(ref, ident.1=st, only.pos=TRUE, logfc.threshold=0.5, max.cells.per.ident=100)
  markers_list[[st]] <- head(rownames(deg), 20)
  cat(sprintf("%-30s %d markers\n", st, length(markers_list[[st]])))
}

SAMPLES <- c("Pla_Camb9518737","Pla_HDBR9518710","WS_PLA_S9101764",
  "WS_PLA_S9101765","WS_PLA_S9101766","WS_PLA_S9101767",
  "WS_PLA_S9101769","WS_PLA_S9101770")

for(sid in SAMPLES) {
  sp_file <- file.path(SPDIR, paste0("spatial_", sid, ".rds"))
  if(!file.exists(sp_file)) { cat(sprintf("SKIP: %s\n", sid)); next }
  sp <- readRDS(sp_file)
  
  # AddModuleScore for each subtype
  for(st in names(markers_list)) {
    genes <- intersect(markers_list[[st]], rownames(sp))
    if(length(genes) < 5) next
    score_name <- make.names(st)
    sp <- AddModuleScore(sp, features=list(genes), name=score_name, ctrl=min(50,nrow(sp)))
  }
  
  # Rename module columns
  for(st in names(markers_list)) {
    old <- paste0(make.names(st), "1")
    if(old %in% colnames(sp@meta.data)) {
      colnames(sp@meta.data)[colnames(sp@meta.data)==old] <- paste0("Spatial_", st)
    }
  }
  
  # Plot
  score_cols <- grep("^Spatial_", colnames(sp@meta.data), value=TRUE)
  if(length(score_cols) > 0) {
    for(sc in score_cols) {
      p <- SpatialFeaturePlot(sp, features=sc, alpha=c(0.1,1), max.cutoff="q95") +
        labs(title=paste(sid, "-", gsub("Spatial_","",sc)))
      ggsave(file.path(SPDIR, sprintf("module_%s_%s.png", sid, gsub("Spatial_","",sc))),
             p, w=5, h=5, dpi=300, bg="white")
    }
  }
  
  saveRDS(sp, file.path(SPDIR, paste0("spatial_modules_", sid, ".rds")))
  cat(sprintf("DONE: %s\n", sid))
}
cat("\nSpatial module mapping complete\n")
