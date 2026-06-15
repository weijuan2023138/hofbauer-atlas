#!/usr/bin/env Rscript
# 补充FigA: Exact format match — 5×2 classifier marker UMAP grid
library(Seurat); library(ggplot2); library(patchwork)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
seu$UMAP_1 <- Embeddings(seu, "umap")[,1]
seu$UMAP_2 <- Embeddings(seu, "umap")[,2]
set.seed(42); idx <- sample(ncol(seu), min(20000, ncol(seu)))

# Theme matching original: no grid, no axis text, compact
tpub <- theme_bw(base_size=10) +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(color="black",fill=NA,linewidth=0.3),
        axis.text=element_blank(), axis.ticks=element_blank(),
        axis.title=element_text(size=8),
        legend.position="right", legend.title=element_blank(),
        legend.text=element_text(size=7), legend.key.size=unit(0.25,"cm"),
        plot.title=element_text(size=11,face="bold",hjust=0.5),
        plot.margin=margin(3,3,3,3))

# Top row: Hofbauer-up genes
row1 <- c("FOLR2","CD163","DAB2","MAF","F13A1")
# Bottom row: Maternal macrophage-up + other key genes
row2 <- c("HLA-DRA","HLA-DRB1","CD74","FCGR3A","TREM2")

make_umap <- function(gene) {
  if(!gene %in% rownames(seu)) return(NULL)
  df <- cbind(seu@meta.data[idx,], FetchData(seu, vars=gene, cells=colnames(seu)[idx]))
  ggplot(df, aes(UMAP_1, UMAP_2, color=.data[[gene]])) +
    geom_point(size=0.12, alpha=0.7) +
    scale_color_gradientn(colors=c("grey90","#BDD7E7","#2171B5","#08306B")) +
    labs(title=gene, x="UMAP 1", y="UMAP 2") + tpub +
    theme(plot.title=element_text(size=14, face="bold"))
}

plots1 <- lapply(row1, make_umap)
plots2 <- lapply(row2, make_umap)

# Stack: top row + bottom row
p_all <- wrap_plots(c(plots1, plots2), ncol=5)
ggsave(file.path(FIGDIR,"补充FigA_markers_classifier_UMAP.png"), p_all, w=22, h=9, dpi=300, bg="white")
cat("Saved: 补充FigA_markers_classifier_UMAP.png\n")
