#!/usr/bin/env Rscript
# Fig3e: ATAC-RNA joint scatter (simple nearest-gene annotation)
library(ggplot2); library(dplyr); library(ggrepel)

OUTDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results'
FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'

# Load ATAC DA
da <- read.csv(file.path(OUTDIR, "atac_differential_peaks.csv"))
da <- da[!is.na(da$p_val_adj),]

suppressMessages(library(Signac))
suppressMessages(library(EnsDb.Hsapiens.v86))
annot <- genes(EnsDb.Hsapiens.v86)
annot <- keepStandardChromosomes(annot, pruning.mode="coarse")
seqlevels(annot) <- paste0("chr", seqlevels(annot))
parts <- strsplit(da$X, "-")
peaks_gr <- GRanges(
  seqnames = sapply(parts, `[`, 1),
  ranges = IRanges(start=as.integer(sapply(parts, `[`, 2)),
                   end=as.integer(sapply(parts, `[`, 3)))
)
nearest <- distanceToNearest(peaks_gr, annot)
gene_ids <- rep(NA_character_, length(peaks_gr))
gene_ids[queryHits(nearest)] <- names(annot)[subjectHits(nearest)]
da$gene <- gene_ids
da$gene_symbol <- NA_character_
valid_genes <- da$gene[!is.na(da$gene)]
valid_symbols <- mapIds(EnsDb.Hsapiens.v86, keys=valid_genes, column="SYMBOL", keytype="GENEID")
da$gene_symbol[!is.na(da$gene)] <- valid_symbols

cat(sprintf("Peaks annotated: %d/%d\n", sum(!is.na(da$gene_symbol)), nrow(da)))

rna <- read.csv(file.path(OUTDIR, "dev_trimester_DEGs.csv"))
rna_late <- rna[rna$trimester=="Late", c("gene","avg_log2FC","p_val_adj")]
merged <- merge(da, rna_late, by.x="gene_symbol", by.y="gene", suffixes=c("_atac","_rna"))
cat(sprintf("ATAC-RNA matched: %d genes\n", nrow(merged)))

# Classify concordance: same direction + ATAC sig
merged$sig <- "NS"
merged$sig[merged$p_val_adj_atac<0.05 & merged$avg_log2FC_atac>0.25 & merged$avg_log2FC_rna>0] <- "Concordant"
merged$sig[merged$p_val_adj_atac<0.05 & merged$avg_log2FC_atac< -0.25 & merged$avg_log2FC_rna<0] <- "Concordant"

key <- c("NFKB1","RELB","NR4A3","NFKB2","CCL8","CTSS","HLA-DRA","FCGR3A","IL1B","C1QA","SOD2","CD74","CXCL8","IRF1","STAT1","JUNB","FOS","EGR1")
merged$label <- ifelse(merged$gene_symbol %in% key & merged$sig=="Concordant", merged$gene_symbol, "")

p <- ggplot(merged, aes(x=avg_log2FC_rna, y=avg_log2FC_atac, color=sig)) +
  geom_hline(yintercept=0, color="grey80", linewidth=0.3) +
  geom_vline(xintercept=0, color="grey80", linewidth=0.3) +
  geom_point(size=0.6, alpha=0.4) +
  geom_point(data=subset(merged, sig=="Concordant"), size=1, alpha=0.8) +
  geom_text_repel(aes(label=label), size=3, max.overlaps=25, color="black", fontface="italic") +
  scale_color_manual(values=c("Concordant"="#D73027","NS"="grey80")) +
  labs(x="RNA log2FC (Late vs Early)", y="ATAC log2FC (Term vs Mid)",
       title="ATAC-RNA Concordance",
       subtitle=sprintf("%d genes, %.0f%% same direction", nrow(merged),
         100*sum(sign(merged$avg_log2FC_atac)==sign(merged$avg_log2FC_rna))/nrow(merged))) +
  theme_bw() + theme(panel.grid=element_blank(), panel.border=element_rect(color="black",linewidth=0.4),
    legend.position="none")

ggsave(file.path(FIGDIR,"Fig3e_atac_rna_scatter.png"), p, w=6.5, h=5.5, dpi=300)
ggsave("figures/Fig2/Fig2c_atac_rna_scatter.png", p, w=6.5, h=5.5, dpi=300)
cat("\nSaved Fig3e_atac_rna_scatter.png\n")
cat(sprintf("Concordant: %d genes\n", sum(merged$sig=="Concordant")))
