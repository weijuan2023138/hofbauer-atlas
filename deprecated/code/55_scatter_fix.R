#!/usr/bin/env Rscript
# Scatter: match MASTER — solid points, percentages, LHX2/SPATS2L/CRYL1
library(ggplot2); library(dplyr); library(ggrepel)
suppressMessages(library(Signac)); suppressMessages(library(EnsDb.Hsapiens.v86))

da <- read.csv("results/atac_differential_peaks.csv"); da <- da[!is.na(da$p_val_adj),]
annot <- genes(EnsDb.Hsapiens.v86); annot <- keepStandardChromosomes(annot, pruning.mode="coarse")
seqlevels(annot) <- paste0("chr", seqlevels(annot))
parts <- strsplit(da$X, "-")
peaks_gr <- GRanges(seqnames=sapply(parts,`[`,1), ranges=IRanges(start=as.integer(sapply(parts,`[`,2)), end=as.integer(sapply(parts,`[`,3))))
nearest <- distanceToNearest(peaks_gr, annot)
gene_ids <- rep(NA_character_, length(peaks_gr))
gene_ids[queryHits(nearest)] <- names(annot)[subjectHits(nearest)]
da$gene_symbol <- NA_character_
valid <- !is.na(gene_ids)
da$gene_symbol[valid] <- mapIds(EnsDb.Hsapiens.v86, keys=gene_ids[valid], column="SYMBOL", keytype="GENEID")

rna <- read.csv("results/dev_trimester_DEGs.csv")
rna_late <- rna[rna$trimester=="Late", c("gene","avg_log2FC","p_val_adj")]
merged <- merge(da, rna_late, by.x="gene_symbol", by.y="gene", suffixes=c("_atac","_rna"))

same_dir <- sign(merged$avg_log2FC_atac) == sign(merged$avg_log2FC_rna)
merged$sig <- ifelse(same_dir, "Concordant", "Discordant")
n_conc <- sum(same_dir); n_disc <- sum(!same_dir)
pct_conc <- round(100*n_conc/nrow(merged))
pct_disc <- round(100*n_disc/nrow(merged))

# Labels — include LHX2, SPATS2L, CRYL1
label_genes <- c("VNN2","TAPT1","MAP3K5","DEXI","KLF3","LHX2","SPATS2L","CRYL1","SPP1","EMB","KCNQ1")
label_genes <- intersect(label_genes, merged$gene_symbol)
merged$label <- ""
for(g in label_genes) {
  idx <- which(merged$gene_symbol == g)
  if(length(idx) > 0) idx <- idx[which.max(abs(merged$avg_log2FC_atac[idx]))]
  merged$label[idx] <- g
}

cat(sprintf("Concordant=%d Discordant=%d Labels=%d\n", n_conc, n_disc, sum(merged$label!="")))

p <- ggplot(merged, aes(x=avg_log2FC_rna, y=avg_log2FC_atac, color=sig)) +
  geom_hline(yintercept=0, color="grey80", linewidth=0.3) +
  geom_vline(xintercept=0, color="grey80", linewidth=0.3) +
  geom_point(size=0.6, alpha=0.4) +
  geom_text_repel(aes(label=label), size=3, color="black", fontface="italic",
    segment.size=0, box.padding=0.05, force=3, max.overlaps=Inf) +
  scale_color_manual(values=c("Concordant"="#D73027","Discordant"="#4575B4"),
    labels=c(sprintf("Concordant (%d%%)",pct_conc), sprintf("Discordant (%d%%)",pct_disc))) +
  labs(x="RNA log2FC (Late vs Early)", y="ATAC log2FC (Term vs Mid)",
       title="ATAC-RNA Concordance") +
  theme_bw() + theme(panel.grid=element_blank(),
    legend.position="right", legend.text=element_text(face="bold"), legend.title=element_blank(),
    plot.title=element_text(face="bold", hjust=0.5))

ggsave("figures/Fig2/Fig2c_atac_rna_scatter.png", p, w=8, h=5.5, dpi=300)
ggsave("figures/Fig3/Fig3e_atac_rna_scatter.png", p, w=8, h=5.5, dpi=300)
cat("Saved\n")
