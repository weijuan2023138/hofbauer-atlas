#!/usr/bin/env Rscript
# Regenerate ALL Fig1 UMAP figures for new 10-dataset Atlas
library(Seurat); library(ggplot2); library(patchwork); library(dplyr)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"
dir.create(FIGDIR, showWarnings=FALSE, recursive=TRUE)

seu <- readRDS(INPUT)
seu$UMAP_1 <- Embeddings(seu, "umap")[,1]
seu$UMAP_2 <- Embeddings(seu, "umap")[,2]
seu$subtype <- seu$subtype_pred

# ── Colors ──
subtype_colors <- c(
  "Pro-inflammatory"="#C62828", "MHCII+ Antigen-presenting"="#E65100",
  "Homeostatic"="#1565C0", "PRKN+ Autophagy"="#6A1B9A",
  "Vascular remodeling"="#2E7D32", "MKI67+ Proliferating"="#455A64"
)
tri_cols <- c("Early"="#4575B4", "Mid"="#FDAE61", "Late"="#D73027")

# ── Theme ──
tpub <- theme_bw(base_size=11) +
  theme(panel.grid=element_line(color="grey92",linewidth=0.2),
        panel.border=element_rect(color="black",fill=NA,linewidth=0.5),
        axis.text=element_blank(), axis.ticks=element_blank(),
        legend.position="right", legend.title=element_blank(),
        legend.text=element_text(size=9), legend.key=element_blank(),
        plot.title=element_text(size=13,face="bold"))

set.seed(42); idx <- sample(ncol(seu))

# ── Trimester mapping ──
trimester_map <- c(
  "E-MTAB-12421"="Early", "E-MTAB-6701"="Early",
  "E-MTAB-12795"="Early", "GSE214607"="Early",
  "UCSF Li 2026"="Mid",
  "GSE290578"="Late", "GSE333257"="Late", "GSE298602"="Late",
  "GSE298119"="Late", "GSE173193"="Late"
)
ds_vec <- setNames(trimester_map[as.character(seu$dataset)], colnames(seu))
seu <- AddMetaData(seu, metadata=data.frame(
  trimester=factor(ds_vec, levels=c("Early","Mid","Late")),
  row.names=colnames(seu)
))

# ═══ Fig1B: Subtype UMAP ═══
p1 <- ggplot(seu@meta.data[idx,], aes(UMAP_1, UMAP_2, color=subtype)) +
  geom_point(size=0.15, alpha=0.8) + scale_color_manual(values=subtype_colors) +
  labs(title="Hofbauer Atlas — 6 Subtypes", x="UMAP 1", y="UMAP 2") + tpub +
  guides(color=guide_legend(override.aes=list(size=3)))
ggsave(file.path(FIGDIR,"Fig1B_subtype.png"), p1, w=9, h=7, dpi=300, bg="white")
cat("Fig1B_subtype ✓\n")

# ═══ Fig1C: Trimester UMAP ═══
p2 <- ggplot(seu@meta.data[idx,], aes(UMAP_1, UMAP_2, color=trimester)) +
  geom_point(size=0.15, alpha=0.8) + scale_color_manual(values=tri_cols) +
  labs(title="Hofbauer Atlas — Trimesters", x="UMAP 1", y="UMAP 2") + tpub +
  guides(color=guide_legend(override.aes=list(size=3)))
ggsave(file.path(FIGDIR,"Fig1C_trimesters_UMAP.png"), p2, w=9, h=7, dpi=300, bg="white")
cat("Fig1C_trimesters ✓\n")

# ═══ 补充FigA: Classifier marker UMAPs ═══
DefaultAssay(seu) <- "RNA"
classifier_genes <- c("FOLR2","CD163","DAB2","MAF","TREM2","F13A1",
                      "HLA-DRA","HLA-DRB1","CD74",
                      "IL1B","CCL3","CCL4","EGR3",
                      "SPP1","PAPPA","FN1","MMP9","NOTUM",
                      "MKI67","TOP2A","BIRC5",
                      "PRKN","C9","SQSTM1","HSPA1A")

p3_list <- list()
for(g in classifier_genes) {
  if(!g %in% rownames(seu)) next
  df <- cbind(seu@meta.data[idx,], FetchData(seu, vars=g, cells=colnames(seu)[idx]))
  p <- ggplot(df, aes(UMAP_1, UMAP_2, color=.data[[g]])) +
    geom_point(size=0.1, alpha=0.5) +
    scale_color_gradientn(colors=c("grey90","#BDD7E7","#2171B5","#08306B")) +
    labs(title=g, x="UMAP 1", y="UMAP 2") + tpub +
    theme(legend.title=element_text(size=8))
  ggsave(file.path(FIGDIR, paste0("补充FigA_", g, "_UMAP.png")), p, w=6, h=5.5, dpi=300, bg="white")
}
cat(sprintf("补充FigA: %d gene UMAPs ✓\n", length(classifier_genes)))

# ═══ Fig1A Overview: dataset gestational timeline ═══
datasets <- data.frame(
  name = c("E-MTAB-6701","E-MTAB-12421","E-MTAB-12795","GSE214607",
           "UCSF Li 2026","GSE290578","GSE298602","GSE333257",
           "GSE298119","GSE173193"),
  gw_start = c(6,4.5,4,6,11,29,37,32,37,37),
  gw_end   = c(12,12.5,8.5,8,39,40,40,39,40,40),
  type = c("Normal","Normal","Infection","Miscarriage","Normal",
           "Normal+PE","PE/Control","Preterm","PE","PE"),
  stringsAsFactors=FALSE
)
datasets$type <- factor(datasets$type,
  levels=c("Normal","Normal+PE","PE/Control","PE","Miscarriage","Infection","Preterm"))
type_cols <- c("Normal"="#4575B4","Normal+PE"="#91BFDB","PE/Control"="#FC8D59",
               "PE"="#FC8D59","Miscarriage"="#D73027","Infection"="#FDB462","Preterm"="#E41A1C")

datasets <- datasets[order(datasets$gw_start),]
datasets$y_pos <- nrow(datasets):1

pA <- ggplot(datasets) +
  geom_segment(aes(x=gw_start, xend=gw_end, y=y_pos, yend=y_pos, color=type),
               linewidth=4.5, alpha=0.88) +
  geom_text(aes(x=3.5, y=y_pos, label=name), hjust=1, size=3, lineheight=0.9) +
  scale_color_manual(values=type_cols, name="Condition") +
  scale_x_continuous(breaks=seq(5,40,5), limits=c(0,42), labels=paste0("GW",seq(5,40,5))) +
  scale_y_continuous(limits=c(0.5, nrow(datasets)+0.5)) +
  labs(x="Gestational Age", y="") +
  theme_minimal() + theme(text=element_text(size=11), axis.text.y=element_blank(),
    axis.ticks.y=element_blank(), axis.text.x=element_text(color="black",size=9),
    axis.title.x=element_text(size=10), legend.position="top",
    legend.title=element_text(size=9), legend.text=element_text(size=8),
    legend.key.size=unit(0.3,"cm"), panel.grid.major.y=element_blank(),
    panel.grid.minor=element_blank(),
    panel.grid.major.x=element_line(color="grey90",linewidth=0.3))
ggsave(file.path(FIGDIR,"Fig1A_data_overview.png"), pA, w=10, h=5.5, dpi=300, bg="white")
cat("Fig1A_data_overview ✓\n")

cat("\nAll figures saved to:", FIGDIR, "\n")
