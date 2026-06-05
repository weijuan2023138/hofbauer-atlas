#!/usr/bin/env Rscript
# Fix scatter plot: Concordant + Discordant + 92%
library(ggplot2); library(dplyr); library(ggrepel)
suppressMessages(library(Signac)); suppressMessages(library(EnsDb.Hsapiens.v86))

# ── Load & annotate ATAC peaks ──
da <- read.csv("results/atac_differential_peaks.csv")
da <- da[!is.na(da$p_val_adj),]
annot <- genes(EnsDb.Hsapiens.v86)
annot <- keepStandardChromosomes(annot, pruning.mode="coarse")
seqlevels(annot) <- paste0("chr", seqlevels(annot))
parts <- strsplit(da$X, "-")
peaks_gr <- GRanges(
  seqnames = sapply(parts, `[`, 1),
  ranges   = IRanges(start=as.integer(sapply(parts, `[`, 2)),
                      end=as.integer(sapply(parts, `[`, 3))))
nearest <- distanceToNearest(peaks_gr, annot)
gene_ids <- rep(NA_character_, length(peaks_gr))
gene_ids[queryHits(nearest)] <- names(annot)[subjectHits(nearest)]
da$gene_symbol <- NA_character_
valid <- !is.na(gene_ids)
da$gene_symbol[valid] <- mapIds(EnsDb.Hsapiens.v86, keys=gene_ids[valid],
                                 column="SYMBOL", keytype="GENEID")

# ── Merge with RNA DEGs ──
rna <- read.csv("results/dev_trimester_DEGs.csv")
rna_late <- rna[rna$trimester == "Late", c("gene", "avg_log2FC", "p_val_adj")]
merged <- merge(da, rna_late, by.x="gene_symbol", by.y="gene",
                suffixes=c("_atac", "_rna"))

# ── Classify ──
same_dir <- sign(merged$avg_log2FC_atac) == sign(merged$avg_log2FC_rna)
pct_same  <- round(100 * sum(same_dir) / nrow(merged))

atac_sig  <- merged$p_val_adj_atac < 0.05 & abs(merged$avg_log2FC_atac) > 0.25
merged$sig <- "NS"
merged$sig[atac_sig & same_dir]  <- "Concordant"
merged$sig[atac_sig & !same_dir] <- "Discordant"
n_conc <- sum(merged$sig == "Concordant")
n_disc <- sum(merged$sig == "Discordant")

# ── Labels ──
key <- c("NFKB1","RELB","NR4A3","CCL8","CTSS","HLA-DRA","FCGR3A",
         "C1QA","SOD2","CD74","CXCL8")
merged$label <- ifelse(merged$gene_symbol %in% key & merged$sig == "Concordant",
                       merged$gene_symbol, "")

# ── Plot ──
p <- ggplot(merged, aes(x=avg_log2FC_rna, y=avg_log2FC_atac, color=sig)) +
  geom_hline(yintercept=0, color="grey80", linewidth=0.3) +
  geom_vline(xintercept=0, color="grey80", linewidth=0.3) +
  geom_point(size=0.6, alpha=0.4) +
  geom_text_repel(aes(label=label), size=3, max.overlaps=25,
                  color="black", fontface="italic") +
  scale_color_manual(
    values = c("Concordant"="#D73027", "Discordant"="#4575B4", "NS"="grey80"),
    labels = c(sprintf("Concordant (%d)", n_conc),
               sprintf("Discordant (%d)", n_disc))) +
  labs(x="RNA log2FC (Late vs Early)",
       y="ATAC log2FC (Term vs Mid)",
       title="ATAC-RNA Concordance",
       subtitle=sprintf("%d genes, %d%% same direction", nrow(merged), pct_same)) +
  theme_bw() +
  theme(panel.grid    = element_blank(),
        legend.position  = "right",
        legend.text      = element_text(face="bold"),
        legend.title     = element_blank(),
        plot.title       = element_text(face="bold", hjust=0.5),
        plot.subtitle    = element_text(hjust=0.5))

ggsave("figures/Fig2/Fig2c_atac_rna_scatter.png", p, w=7, h=5.5, dpi=300)
ggsave("figures/Fig3/Fig3e_atac_rna_scatter.png", p, w=7, h=5.5, dpi=300)

cat(sprintf("Concordant=%d  Discordant=%d  NS=%d  SameDir=%d%%\n",
            n_conc, n_disc, sum(merged$sig=="NS"), pct_same))
