#!/usr/bin/env Rscript
# Fig6b: TSS-centered ATAC mirror plot with gene models (bp units)
library(EnsDb.Hsapiens.v86); library(GenomicRanges); library(ggplot2); library(patchwork); library(dplyr)

atac <- readRDS("results/Hofbauer_ATAC_mid_term.rds")
annot <- genes(EnsDb.Hsapiens.v86); annot <- keepStandardChromosomes(annot, pruning.mode="coarse")
seqlevels(annot) <- paste0("chr", seqlevels(annot))

frag_dir <- "/home/weijuan/文档/胎盘单细胞数据/raw_data/UCSF_Li_2026/fragments"
mid_samples <- c("ZY012","ZY014","ZY020")
term_samples <- c("ZY011","ZY019")

read_cov <- function(gene_symbol, flank=5000, bin=30) {
  gene <- annot[annot$symbol==gene_symbol]
  tss <- ifelse(as.character(strand(gene)[1])=="-", end(gene)[1], start(gene)[1])
  rs <- max(1, tss-flank); re <- tss+flank
  chr <- as.character(seqnames(gene)[1])
  region <- sprintf("%s:%d-%d", chr, rs, re)
  bins <- seq(rs, re, by=bin)
  get_one <- function(samples, label) {
    cov <- rep(0, length(bins)-1)
    for(s in samples) {
      path <- file.path(frag_dir, paste0(s, "_fragments.tsv.bgz"))
      cmd <- sprintf("tabix %s %s 2>/dev/null", path, region)
      lines <- system(cmd, intern=TRUE)
      if(length(lines)==0) next
      parts <- strsplit(lines, "\t")
      starts <- as.integer(sapply(parts, `[`, 2))
      ends <- as.integer(sapply(parts, `[`, 3))
      for(i in seq_along(starts)) {
        si <- max(1, ceiling((starts[i]-rs)/bin))
        ei <- min(length(cov), floor((ends[i]-rs)/bin))
        if(si <= ei) cov[si:ei] <- cov[si:ei] + 1
      }
    }
    cov
  }
  # Normalize by cell count: Mid=1230, Term=925
  mid_cov <- get_one(mid_samples) / 1230
  term_cov <- get_one(term_samples) / 925
  list(mid=mid_cov, term=term_cov, bins=bins, rs=rs, re=re)
}

make_gene_track <- function(gene_symbol, rs, re) {
  gene <- annot[annot$symbol==gene_symbol]
  ex <- exons(EnsDb.Hsapiens.v86, columns=c("exon_id","gene_name"),
    filter=GeneNameFilter(gene_symbol))
  ex <- keepStandardChromosomes(ex, pruning.mode="coarse")
  seqlevels(ex) <- paste0("chr", seqlevels(ex))
  if(length(ex)==0) return(ggplot())
  ex_df <- as.data.frame(ex) %>%
    filter(seqnames==as.character(seqnames(gene)[1]), end>=rs, start<=re)
  strand_dir <- ifelse(as.character(strand(gene)[1])=="-", -1, 1)
  ggplot(ex_df) +
    geom_rect(aes(xmin=start, xmax=end, ymin=-0.2, ymax=0.2), fill="grey30", color=NA) +
    annotate("segment", x=min(ex_df$start), xend=max(ex_df$end), y=0, yend=0,
      linewidth=0.4, color="grey50") +
    annotate("segment", x=start(gene)[1], xend=start(gene)[1]+strand_dir*1500,
      y=0, yend=0, arrow=arrow(length=unit(0.06,"inch"),type="closed"),
      linewidth=0.6, color="grey30") +
    xlim(rs, re) + ylim(-0.4, 0.4) + theme_void()
}

build_panel <- function(gene_symbol, flank=5000, bin=30) {
  cov <- read_cov(gene_symbol, flank, bin)
  pos_bp <- (cov$bins[-1] + cov$bins[-length(cov$bins)]) / 2
  df <- data.frame(pos=pos_bp, Mid=cov$mid, Term=-cov$term)
  ymax <- max(abs(c(df$Mid, df$Term)), na.rm=TRUE) * 1.1
  
  # Format x-axis: show kb from TSS
  tss <- cov$rs + flank  # TSS position in bp
  x_breaks <- seq(cov$rs, cov$re, length.out=5)
  x_labels <- sprintf("%.0f kb", (x_breaks - tss) / 1000)
  
  p_main <- ggplot(df) +
    geom_ribbon(aes(x=pos, ymin=0, ymax=Mid), fill="#FDAE61", alpha=0.85) +
    geom_ribbon(aes(x=pos, ymin=0, ymax=Term), fill="#D73027", alpha=0.85) +
    geom_hline(yintercept=0, linewidth=0.4, color="grey30") +
    coord_cartesian(xlim=c(cov$rs, cov$re), expand=FALSE) +
    annotate("text", x=cov$re - (cov$re-cov$rs)*0.05, y=ymax*0.7, label="Mid", color="#FDAE61",
      fontface="bold", size=5, hjust=1) +
    annotate("text", x=cov$re - (cov$re-cov$rs)*0.05, y=-ymax*0.7, label="Term", color="#D73027",
      fontface="bold", size=5, hjust=1) +
    scale_x_continuous(breaks=x_breaks, labels=x_labels) +
    scale_y_continuous(limits=c(-ymax, ymax)) +
    labs(title=gene_symbol, x="Distance from TSS (kb)", y="Norm. frag./cell") +
    theme_bw(9) +
    theme(plot.title=element_text(face="bold.italic", size=14, hjust=0.5, color="black"),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=9),
          panel.grid=element_blank())
  
  p_gene <- make_gene_track(gene_symbol, cov$rs, cov$re)
  p_main / p_gene + plot_layout(heights=c(1, 0.12))
}

p1 <- build_panel("TGFB1", flank=10000, bin=50)
p2 <- build_panel("CD47", flank=10000, bin=50)
p3 <- build_panel("PTPRM", flank=10000, bin=50)
p4 <- build_panel("SPP1", flank=10000, bin=50)

combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title="ATAC-seq coverage at communication gene promoters",
    theme=theme(plot.title=element_text(face="bold", size=16, hjust=0.5)))
ggsave("figures/Fig6/Fig6b_ATAC_tracks.png", combined, w=16, h=11, dpi=300, bg="white")
cat("Fig6b done\n")
