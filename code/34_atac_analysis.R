#!/usr/bin/env Rscript
# ATAC v3: per-sample peaks -> merge GRanges -> combined matrix
library(Signac); library(Seurat); library(ggplot2); library(dplyr); library(GenomicRanges)

FRAG_DIR <- '/home/weijuan/文档/胎盘单细胞数据/raw_data/UCSF_Li_2026/fragments'
OUTDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results'
FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'

hb_meta <- read.csv(file.path(OUTDIR, "hb_barcodes.csv"))

# ---- Per-sample peak calling ----
all_peaks <- GRanges()
all_counts_list <- list()
all_frags <- list()

for(sid in unique(hb_meta$sample)) {
  f <- file.path(FRAG_DIR, paste0(sid, "_fragments.tsv.bgz"))
  if(!file.exists(f)) next
  cells <- sub("^[^_]+_", "", hb_meta$barcode[hb_meta$sample==sid])
  cat(sprintf("%s: %d cells...\n", sid, length(cells)))
  
  frag <- CreateFragmentObject(f, cells=cells)
  peaks <- CallPeaks(frag, effective.genome.size=2.7e9, macs2.path="/home/weijuan/.local/bin/macs3")
  cat(sprintf("  peaks: %d\n", length(peaks)))
  
  counts <- FeatureMatrix(fragments=frag, features=peaks, cells=cells)
  all_counts_list[[sid]] <- counts
  all_frags[[sid]] <- frag
  all_peaks <- c(all_peaks, peaks)
}

# Reduce to non-overlapping consensus
all_peaks <- reduce(all_peaks)
cat(sprintf("Consensus peaks: %d\n", length(all_peaks)))

# Re-count all samples against consensus peaks
all_counts <- NULL
all_cells <- c()
cell_groups <- c()
cell_samples <- c()

for(sid in names(all_frags)) {
  cells <- sub("^[^_]+_", "", hb_meta$barcode[hb_meta$sample==sid])
  grp <- hb_meta$group[hb_meta$sample==sid][1]
  counts <- FeatureMatrix(fragments=all_frags[[sid]], features=all_peaks, cells=cells)
  colnames(counts) <- paste0(sid, "_", colnames(counts))  # unique cell names
  
  if(is.null(all_counts)) {
    all_counts <- counts
  } else {
    all_counts <- cbind(all_counts, counts)
  }
  all_cells <- c(all_cells, colnames(counts))
  cell_groups <- c(cell_groups, rep(grp, ncol(counts)))
  cell_samples <- c(cell_samples, rep(sid, ncol(counts)))
}
names(cell_groups) <- all_cells
names(cell_samples) <- all_cells

cat(sprintf("Matrix: %d peaks x %d cells\n", nrow(all_counts), ncol(all_counts)))

# ---- Seurat object ----
combined_frags <- Reduce(c, all_frags)
chrom <- CreateChromatinAssay(counts=all_counts, sep=c(":","-"))
obj <- CreateSeuratObject(chrom, assay="ATAC")
obj$group <- ifelse(cell_groups[colnames(obj)] == "2nd_trimester", "Mid", "Term")
obj$sample <- cell_samples[colnames(obj)]

obj <- RunTFIDF(obj); obj <- FindTopFeatures(obj, min.cutoff='q5')
obj <- RunSVD(obj); obj <- RunUMAP(obj, dims=2:30, reduction='lsi')

p1 <- DimPlot(obj, group.by='group', cols=c("Mid"="#FDAE61","Term"="#D73027"), pt.size=1) +
  ggtitle("Hofbauer ATAC-seq") + theme_bw() + theme(panel.grid=element_blank())
ggsave(file.path(FIGDIR,"Fig3c_atac_umap.png"), p1, w=6, h=5, dpi=300)

# DA
Idents(obj) <- "group"
da <- FindMarkers(obj, ident.1="Term", ident.2="Mid", min.pct=0.05, test.use="LR", latent.vars="nCount_ATAC")
da$sig <- ifelse(!is.na(da$p_val_adj) & da$p_val_adj<0.05 & abs(da$avg_log2FC)>0.25, "Sig", "NS")
cat(sprintf("DA: Term-up=%d Mid-up=%d\n",
  sum(da$sig=="Sig" & da$avg_log2FC>0), sum(da$sig=="Sig" & da$avg_log2FC<0)))

p2 <- ggplot(da, aes(x=avg_log2FC, y=-log10(p_val_adj), color=sig)) +
  geom_point(size=0.5, alpha=0.6) + scale_color_manual(values=c("Sig"="red","NS"="grey80")) +
  labs(title="Differential Accessibility: Term vs Mid") + theme_bw() + theme(legend.position="none")
ggsave(file.path(FIGDIR,"Fig3d_atac_volcano.png"), p2, w=6, h=5, dpi=300)

saveRDS(obj, file.path(OUTDIR, "Hofbauer_ATAC_mid_term.rds"))
write.csv(da, file.path(OUTDIR, "atac_differential_peaks.csv"))
cat("Done\n")
