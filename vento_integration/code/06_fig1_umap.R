#!/usr/bin/env Rscript
# 10-dataset Hofbauer Atlas: clustering → markers → subtype annotation → Fig1 UMAP
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)

INPUT  <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_harmony_10datasets.rds"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

seu <- readRDS(INPUT)
cat(sprintf("Loaded: %d cells, %d clusters (res=0.15)\n", ncol(seu), length(unique(seu$RNA_snn_res.0.15))))

# ── 1. Add disease group labels ──
seu$disease_group <- dplyr::recode(seu$disease_group,
  "Normal 1st trimester" = "Normal 1st",
  "Normal 1st/2nd/Term"  = "Normal",
  "Normal 3rd trimester / Preeclampsia" = "Normal+PE",
  "Preeclampsia" = "PE",
  "Preeclampsia / Control" = "PE/Control",
  "Severe Preeclampsia" = "Severe PE",
  "Preterm Labor / Term Labor" = "Preterm",
  "Miscarriage / Normal" = "Miscarriage",
  "Infection" = "Infection"
)

# Also add dataset_short
seu$dataset_short <- dplyr::recode(seu$dataset,
  "Arutyunyan" = "E-MTAB-12421",
  "Vento_Tormo_2018" = "E-MTAB-6701",
  "hoo_2024" = "E-MTAB-12795",
  "gse214607" = "GSE214607",
  "my_preterm_cohort" = "GSE333257",
  "UCSF_Li_2026" = "UCSF Li 2026",
  "GSE290578" = "GSE290578",
  "gse298602" = "GSE298602",
  "gse298119" = "GSE298119",
  "gse329173" = "GSE329173",
  "gse173193" = "GSE173193"
)

# ── 2. Disease group color palette (from Shiny app.R) ──
disease_cols <- c(
  "Normal 1st" = "#4575B4", "Normal" = "#4575B4",
  "Normal+PE" = "#91BFDB", "PE" = "#FC8D59",
  "PE/Control" = "#FC8D59", "Severe PE" = "#D73027",
  "Miscarriage" = "#D73027", "Infection" = "#FDB462",
  "Preterm" = "#E41A1C"
)

# ── 3. Find cluster markers ──
Idents(seu) <- "RNA_snn_res.0.15"
seu <- JoinLayers(seu)
cat("\nFinding cluster markers...\n")
markers <- FindAllMarkers(seu, only.pos=TRUE, min.pct=0.3, logfc.threshold=0.5, test.use='t')
markers <- markers[!grepl("^(RPL|RPS|MT-|MALAT1|LINC|AC0|AL|AP0|RP11|RP1-|CTD-)", markers$gene), ]

# ── 4. Check for non-Hofbauer clusters (maternal macrophage markers) ──
cat("\nChecking for maternal macrophage contamination...\n")
mat_markers <- c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "CD74", "FCGR3A")
hb_markers <- c("FOLR2", "CD163", "DAB2", "MAF", "TREM2")

for(cl in sort(unique(Idents(seu)))) {
  cells <- WhichCells(seu, idents=cl)
  mat_expr <- mean(FetchData(seu, vars=intersect(mat_markers, rownames(seu)), cells=cells)[,1])
  hb_expr <- mean(FetchData(seu, vars=intersect(hb_markers, rownames(seu)), cells=cells)[,1])
  top3 <- head(markers[markers$cluster==cl, "gene"], 3)
  cat(sprintf("  C%2d: %5d cells | HLA-DRA=%.2f FOLR2=%.2f | top: %s\n",
      cl, length(cells), mat_expr, hb_expr, paste(top3, collapse=", ")))
}

# ── 5. Generate UMAP plots ──
set.seed(42)
idx <- sample(1:ncol(seu), min(30000, ncol(seu)))

tpub <- theme_bw(base_size=11) +
  theme(panel.grid=element_line(color="grey92", linewidth=0.2),
        panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
        axis.text=element_blank(), axis.ticks=element_blank(),
        legend.position="right", legend.title=element_blank(),
        legend.text=element_text(size=8), legend.key=element_blank(),
        plot.title=element_text(size=12, face="bold"))

# UMAP by cluster
p1 <- ggplot(seu@meta.data[idx,], aes(x=umap_1, y=umap_2, color=RNA_snn_res.0.15)) +
  geom_point(size=0.15, alpha=0.6) + tpub +
  ggtitle("Hofbauer Atlas — 17 Clusters (res=0.15)") +
  guides(color=guide_legend(override.aes=list(size=3), ncol=2))
ggsave(file.path(OUTDIR, "UMAP_clusters_10datasets.png"), p1, w=9, h=7, dpi=300, bg="white")

# UMAP by dataset
p2 <- ggplot(seu@meta.data[idx,], aes(x=umap_1, y=umap_2, color=dataset_short)) +
  geom_point(size=0.15, alpha=0.6) + tpub +
  ggtitle("By Dataset") +
  guides(color=guide_legend(override.aes=list(size=3)))
ggsave(file.path(OUTDIR, "UMAP_dataset_10datasets.png"), p2, w=10, h=7, dpi=300, bg="white")

# UMAP by disease group
p3 <- ggplot(seu@meta.data[idx,], aes(x=umap_1, y=umap_2, color=disease_group)) +
  geom_point(size=0.15, alpha=0.6) + tpub +
  scale_color_manual(values=disease_cols, na.value="grey70") +
  ggtitle("By Disease Group") +
  guides(color=guide_legend(override.aes=list(size=3)))
ggsave(file.path(OUTDIR, "UMAP_disease_10datasets.png"), p3, w=10, h=7, dpi=300, bg="white")

cat("\nDone! Figures saved to:", OUTDIR, "\n")
cat("  UMAP_clusters_10datasets.png\n")
cat("  UMAP_dataset_10datasets.png\n")
cat("  UMAP_disease_10datasets.png\n")
