#!/usr/bin/env Rscript
# Fig6c Motif Enrichment: Find ATAC peaks near communication genes
setwd("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration")
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)

comm_genes <- c("SPP1","FN1","COL1A2","COL1A1","TGFB1","IGF1","PTPRM","CD44","ITGAV","ITGB1",
  "ITGA1","ITGA2","ITGB3","ITGB5","CD47","CXCL8","CCL2","CCL3","CCL4","CXCL2",
  "TNF","IL1B","IL6","IL10","CSF1","VEGFA","HGF","BMP2","MMP9","THBS1")

# Gene coordinates
edb <- EnsDb.Hsapiens.v86
gene_ranges <- genes(edb, filter=GeneNameFilter(comm_genes))
gene_ranges <- gene_ranges[seqnames(gene_ranges) %in% c(1:22,"X")]
seqlevels(gene_ranges) <- paste0("chr", seqlevels(gene_ranges))
gene_50k <- flank(gene_ranges, width=50000, both=TRUE)
gene_50k <- gene_50k[!is.na(gene_50k)]

# Load differential peaks
peaks <- read.csv("results/atac_differential_peaks.csv", row.names=1)
gr <- GRanges(
  seqnames = gsub("-.*","", rownames(peaks)),
  ranges = IRanges(
    start = as.integer(gsub(".*-(\\d+)-.*","\\1", rownames(peaks))),
    end   = as.integer(gsub(".*-\\d+-(\\d+)","\\1", rownames(peaks)))
  )
)

# Overlap with comm gene regions
ov <- findOverlaps(gr, gene_50k)
comm_peaks <- gr[unique(queryHits(ov))]
hit_genes <- gene_ranges$gene_name[unique(subjectHits(ov))]

cat(sprintf("Total differential peaks: %d\n", length(gr)))
cat(sprintf("Peaks overlapping comm genes (±50kb): %d\n", length(comm_peaks)))
cat(sprintf("Comm genes hit: %d/%d\n", length(hit_genes), length(comm_genes)))
cat("Genes:", paste(sort(hit_genes), collapse=", "), "\n")

# Write comm-gene peaks to FASTA
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, comm_peaks)
names(seqs) <- paste0(seqnames(comm_peaks),":",start(comm_peaks),"-",end(comm_peaks))
writeXStringSet(seqs, "results/comm_gene_diff_peaks.fa")

cat("Done: results/comm_gene_diff_peaks.fa\n")
