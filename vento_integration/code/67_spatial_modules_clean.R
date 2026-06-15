#!/usr/bin/env Rscript
# Fix spatial module plots — smaller legend, no overlap
library(Seurat); library(ggplot2)

SPDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/spatial"
files <- list.files(SPDIR, pattern="spatial_modules_.*rds$", full=TRUE)

for(rf in files) {
  sp <- readRDS(rf)
  sid <- gsub("spatial_modules_|\\.rds", "", basename(rf))
  score_cols <- grep("^Spatial_", colnames(sp@meta.data), value=TRUE)
  
  for(sc in score_cols) {
    st_name <- gsub("Spatial_", "", sc)
    p <- SpatialFeaturePlot(sp, features=sc, alpha=c(0.1,1), max.cutoff="q95", ncol=1, combine=FALSE)
    p[[1]] <- p[[1]] + labs(title=paste(st_name, "-", sid)) +
      theme(legend.key.size=unit(0.3, "cm"), legend.text=element_text(size=6),
        legend.title=element_text(size=7))
    ggsave(file.path(SPDIR, sprintf("module_%s_%s.png", sid, st_name)),
           p[[1]], w=6.5, h=6, dpi=300, bg="white")
  }
  cat(sprintf("%s ", sid))
}
cat("\nDone\n")
