#!/usr/bin/env Rscript
# Coverage tracks with actual fragment data (single sample)
library(Signac); library(Seurat); library(ggplot2)

FRAG_DIR <- '/home/weijuan/文档/胎盘单细胞数据/raw_data/UCSF_Li_2026/fragments'
FIGDIR   <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'

# Load barcodes
hb <- read.csv('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/hb_barcodes.csv')

# Pick one mid + one term sample
mid_sid <- "ZY020"; term_sid <- "ZY011"

build_obj <- function(sid) {
  cells <- hb$barcode[hb$sample==sid]
  cells_bare <- sub("^[^_]+_", "", cells)
  f <- file.path(FRAG_DIR, paste0(sid, "_fragments.tsv.bgz"))
  frag <- CreateFragmentObject(f, cells=cells_bare)
  obj <- CreateSeuratObject(CreateChromatinAssay(
    counts=FeatureMatrix(fragments=frag, features=CallPeaks(frag, macs2.path="/home/weijuan/.local/bin/macs3"), cells=cells_bare),
    fragments=frag
  ), assay="ATAC")
  obj
}

cat("Building mid...\n"); mid <- build_obj(mid_sid)
cat(sprintf("Mid: %d cells, %d peaks\n", ncol(mid), nrow(mid)))
cat("Building term...\n"); term <- build_obj(term_sid)
cat(sprintf("Term: %d cells, %d peaks\n", ncol(term), nrow(term)))

# Merge on shared peaks
shared <- intersect(rownames(mid), rownames(term))
cat(sprintf("Shared peaks: %d\n", length(shared)))
combined <- merge(mid[shared,], term[shared,])
combined$group <- c(rep("Mid",ncol(mid)), rep("Term",ncol(term)))
combined <- RunTFIDF(combined); combined <- FindTopFeatures(combined, min.cutoff='q5')
combined <- RunSVD(combined)

# Coverage plots
for(g in c("NFKB1","RELB","FCGR3A","C1QA","CD36","TREM2","SOD2")) {
  png(file.path(FIGDIR,paste0("Fig3f_track_",g,".png")), w=10, h=4.5, units="in", res=300)
  try(print(CoveragePlot(combined, region=g, group.by="group", extend.upstream=3000, extend.downstream=3000) +
    scale_fill_manual(values=c("Mid"="#FDAE61","Term"="#D73027")) +
    ggtitle(g) + theme(plot.title=element_text(face="bold.italic",size=14))))
  dev.off()
  cat(g, "\n")
}
cat("Done\n")
