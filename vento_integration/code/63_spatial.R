#!/usr/bin/env Rscript
# Spatial transcriptomics: Hofbauer markers across 8 early-gestation Visium samples
library(Seurat); library(ggplot2); library(patchwork)

SAMPLES <- c("Pla_Camb9518737","Pla_HDBR9518710",
  "WS_PLA_S9101764","WS_PLA_S9101765","WS_PLA_S9101766",
  "WS_PLA_S9101767","WS_PLA_S9101769","WS_PLA_S9101770")

RAW <- "/home/weijuan/文档/胎盘单细胞数据/raw_data/Arutyunyan_2023/E-MTAB-12698_visium"
TMP <- "/tmp/visium_extract"
OUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/spatial"
WD <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/spatial"

dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
dir.create(TMP, showWarnings=FALSE)

markers <- c("FOLR2","CD163","DAB2","FN1","HLA-DRA","FLT1","SPP1","COL1A2")

for(sid in SAMPLES) {
  tar_file <- file.path(RAW, paste0(sid, "_spaceranger_output.tar.gz"))
  if(!file.exists(tar_file)) { cat(sprintf("MISSING: %s\n", sid)); next }
  
  # Extract
  ext_dir <- file.path(TMP, sid)
  dir.create(ext_dir, showWarnings=FALSE, recursive=TRUE)
  system(sprintf("tar -xzf %s -C %s 2>/dev/null", tar_file, ext_dir))
  
  # Find the inner output directory (contains both filtered_feature_bc_matrix and spatial)
  inner_dir <- list.dirs(ext_dir, recursive=FALSE)[1]
  if(is.na(inner_dir)) { cat(sprintf("NO INNER: %s\n", sid)); next }
  sp_dir <- list.dirs(inner_dir, recursive=FALSE)
  sp_dir <- sp_dir[grepl("filtered_feature_bc_matrix|spatial", sp_dir)]
  mat_dir <- inner_dir  # Seurat expects the parent directory
  
  # Load with Seurat
  sp <- Load10X_Spatial(mat_dir, slice=sid)
  sp <- SCTransform(sp, assay="Spatial", verbose=FALSE)
  
  # Plot markers
  avail <- intersect(markers, rownames(sp))
  if(length(avail)==0) { cat(sprintf("NO GENES: %s\n", sid)); next }
  
  pdf(file.path(WD, sprintf("spatial_%s.pdf", sid)), w=14, h=7)
  print(SpatialFeaturePlot(sp, features=avail[1:min(4,length(avail))], ncol=4, alpha=c(0.1,1)))
  if(length(avail)>4) print(SpatialFeaturePlot(sp, features=avail[5:min(8,length(avail))], ncol=4, alpha=c(0.1,1)))
  dev.off()
  
  # Also PNG for easy viewing
  p <- SpatialFeaturePlot(sp, features=avail[1:min(4,length(avail))], ncol=4, alpha=c(0.1,1))
  ggsave(file.path(WD, sprintf("spatial_%s.png", sid)), p, w=16, h=4, dpi=300, bg="white")
  
  saveRDS(sp, file.path(WD, sprintf("spatial_%s.rds", sid)))
  cat(sprintf("DONE: %s (%d genes)\n", sid, length(avail)))
}

cat("\nAll spatial samples processed\n")
