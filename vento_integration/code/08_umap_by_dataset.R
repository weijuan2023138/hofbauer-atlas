#!/usr/bin/env Rscript
# UMAP by dataset to assess Harmony integration quality
library(Seurat); library(ggplot2)

INPUT  <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_harmony_10datasets.rds"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
hb_clusters <- c("0","1","2","4","5","6","7","8","9","16")
seu <- subset(seu, RNA_snn_res.0.15 %in% hb_clusters)
seu$umap_1 <- Embeddings(seu, "umap")[,1]
seu$umap_2 <- Embeddings(seu, "umap")[,2]

# Short names
seu$ds_short <- dplyr::recode(seu$dataset,
  "Arutyunyan"="E-MTAB-12421", "Vento_Tormo_2018"="E-MTAB-6701",
  "hoo_2024"="E-MTAB-12795", "gse214607"="GSE214607",
  "my_preterm_cohort"="GSE333257", "UCSF_Li_2026"="UCSF Li",
  "GSE290578"="GSE290578", "gse298602"="GSE298602",
  "gse298119"="GSE298119", "gse329173"="GSE329173",
  "gse173193"="GSE173193"
)

set.seed(42)
idx <- sample(1:ncol(seu), 30000)

tpub <- theme_bw(base_size=10) +
  theme(panel.grid=element_blank(), panel.border=element_rect(color="black",fill=NA,linewidth=0.5),
        axis.text=element_blank(), axis.ticks=element_blank(),
        legend.position="none", plot.title=element_text(size=11,face="bold"))

# One panel per dataset
datasets <- sort(unique(seu$ds_short))

for(ds in datasets) {
  sub_idx <- which(seu$ds_short[idx] == ds)
  n <- sum(seu$ds_short == ds)
  p <- ggplot(seu@meta.data[idx,], aes(x=umap_1, y=umap_2)) +
    geom_point(size=0.1, alpha=0.15, color="grey90") +
    geom_point(data=seu@meta.data[idx[sub_idx],], aes(x=umap_1, y=umap_2),
               size=0.2, alpha=0.7, color="#D73027") +
    tpub + ggtitle(sprintf("%s (%d cells)", ds, n))
  ggsave(file.path(OUTDIR, paste0("UMAP_ds_", gsub("/","_",ds), ".png")),
         p, w=5, h=5, dpi=300, bg="white")
}

# Combined facet
p_all <- ggplot(seu@meta.data[idx,], aes(x=umap_1, y=umap_2)) +
  geom_point(size=0.08, alpha=0.4, color="grey70") +
  facet_wrap(~ds_short, ncol=4) +
  geom_point(size=0.08, alpha=0.6, color="#D73027") +
  tpub + ggtitle("By Dataset")
ggsave(file.path(OUTDIR, "UMAP_by_dataset_facet.png"), p_all, w=16, h=12, dpi=300, bg="white")
cat("Done!\n")
