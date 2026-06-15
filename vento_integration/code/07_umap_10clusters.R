#!/usr/bin/env Rscript
# UMAP for 10 Hofbauer clusters to assess merging
library(Seurat); library(ggplot2); library(patchwork)

INPUT  <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_harmony_10datasets.rds"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)

# Keep only Hofbauer clusters
hb_clusters <- c("0","1","2","4","5","6","7","8","9","16")
seu <- subset(seu, RNA_snn_res.0.15 %in% hb_clusters)
seu$umap_1 <- Embeddings(seu, "umap")[,1]
seu$umap_2 <- Embeddings(seu, "umap")[,2]
cat(sprintf("Hofbauer: %d cells\n", ncol(seu)))

set.seed(42)
idx <- sample(1:ncol(seu), min(30000, ncol(seu)))

tpub <- theme_bw(base_size=11) +
  theme(panel.grid=element_line(color="grey92", linewidth=0.2),
        panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
        axis.text=element_blank(), axis.ticks=element_blank(),
        legend.position="right", legend.title=element_blank(),
        legend.text=element_text(size=8), legend.key=element_blank(),
        plot.title=element_text(size=12,face="bold"))

# ── UMAP by cluster ──
p1 <- ggplot(seu@meta.data[idx,], aes(x=umap_1, y=umap_2, color=RNA_snn_res.0.15)) +
  geom_point(size=0.15, alpha=0.6) + tpub +
  ggtitle("10 Hofbauer Clusters") +
  guides(color=guide_legend(override.aes=list(size=3), ncol=2))
ggsave(file.path(OUTDIR, "UMAP_10clusters.png"), p1, w=9, h=7, dpi=300, bg="white")

# ── Key marker genes on UMAP ──
markers <- c("FOLR2","CD163","IL1B","CCL3","MKI67","TOP2A",
             "SPP1","PAPPA","HSPA1A","HSP90AA1","HLA-DRA","HLA-DQB1",
             "PRKN","TREM2","APOE","C1QB","FN1","MMP9")

for(g in markers) {
  if(!g %in% rownames(seu)) next
  expr <- FetchData(seu, vars=g, cells=colnames(seu)[idx])
  df <- cbind(seu@meta.data[idx,], expr)
  p <- ggplot(df, aes(x=umap_1, y=umap_2, color=.data[[g]])) +
    geom_point(size=0.1, alpha=0.5) +
    scale_color_gradientn(colors=c("grey90","#BDD7E7","#2171B5","#08306B"), name=g) +
    tpub + ggtitle(g)
  ggsave(file.path(OUTDIR, paste0("UMAP_", g, ".png")), p, w=7, h=6, dpi=300, bg="white")
}
cat(sprintf("Saved %d gene UMAPs to %s\n", length(markers), OUTDIR))
