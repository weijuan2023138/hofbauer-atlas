#!/usr/bin/env Rscript
# PE subgroup UMAP visualizations
library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
set.seed(42)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
detail <- labels$disease_detail[1:ncol(seu)]

# PE subgroup assignment
seu$pe_group <- "Other"
seu$pe_group[seu$dataset=="GSE290578" & detail=="PE"] <- "Early_PE\n(GW29-34)"
seu$pe_group[seu$dataset=="GSE298602" & detail=="PreE_SF"] <- "Late_Severe\n(GW37-40)"
seu$pe_group[seu$dataset=="GSE298602" & detail=="gHTN"] <- "Late_Mild\n(gHTN)"
seu$pe_group[seu$dataset=="GSE298602" & detail=="Control"] <- "Control\n(GSE298602)"
seu$pe_group[detail %in% c("GSE173193","GSE298119")] <- "Late_PE\n(term)"
seu$pe_group[detail=="Normal"] <- "Normal_Late"

target <- c("Early_PE\n(GW29-34)","Late_Severe\n(GW37-40)","Late_Mild\n(gHTN)",
            "Late_PE\n(term)","Normal_Late","Control\n(GSE298602)")
seu <- subset(seu, pe_group %in% target)
seu$pe_group <- factor(seu$pe_group, levels=target)

seu$UMAP_1 <- Embeddings(seu, "umap")[,1]
seu$UMAP_2 <- Embeddings(seu, "umap")[,2]

# Colors
pe_cols <- c("Early_PE\n(GW29-34)"="#D73027","Late_Severe\n(GW37-40)"="#FC8D59",
             "Late_Mild\n(gHTN)"="#FDAE61","Late_PE\n(term)"="#FDB462",
             "Normal_Late"="#4575B4","Control\n(GSE298602)"="#66C2A5")
st_cols <- c("Pro-inflammatory"="#C62828","MHCII+ Antigen-presenting"="#E65100",
             "Homeostatic"="#1565C0","PRKN+ Autophagy"="#6A1B9A",
             "Vascular remodeling"="#2E7D32","MKI67+ Proliferating"="#455A64")
tpub <- theme_bw(11) + theme(panel.grid=element_line(color="grey92",linewidth=0.2),
  panel.border=element_rect(color="black",fill=NA,linewidth=0.5),
  axis.text=element_blank(), axis.ticks=element_blank(), legend.position="right",
  legend.title=element_blank(), legend.text=element_text(size=8), legend.key=element_blank(),
  plot.title=element_text(size=12,face="bold",hjust=0.5))

set.seed(42); idx <- sample(ncol(seu), min(20000, ncol(seu)))

# Panel A: UMAP by PE group
pA <- ggplot(seu@meta.data[idx,], aes(UMAP_1, UMAP_2, color=pe_group)) +
  geom_point(size=0.15, alpha=0.8) + scale_color_manual(values=pe_cols) +
  labs(title="PE Subgroups: UMAP by Group", x="UMAP 1", y="UMAP 2") + tpub +
  guides(color=guide_legend(override.aes=list(size=3), ncol=1))

# Panel B: UMAP by subtype (all PE cells)
pe_cells <- seu$pe_group != "Normal_Late" & seu$pe_group != "Control\n(GSE298602)"
pB <- ggplot(seu@meta.data[idx & pe_cells[idx],], aes(UMAP_1, UMAP_2, color=subtype_pred)) +
  geom_point(size=0.15, alpha=0.8) + scale_color_manual(values=st_cols) +
  labs(title="PE Cells Only: UMAP by Subtype", x="UMAP 1", y="UMAP 2") + tpub +
  guides(color=guide_legend(override.aes=list(size=3)))

# Panel C: Key gene feature UMAPs
gene_umaps <- function(g) {
  if(!g %in% rownames(seu)) return(NULL)
  df <- cbind(seu@meta.data[idx,], FetchData(seu, vars=g, cells=colnames(seu)[idx]))
  ggplot(df, aes(UMAP_1, UMAP_2, color=.data[[g]])) +
    geom_point(size=0.12, alpha=0.7) +
    scale_color_gradientn(colors=c("grey90","#BDD7E7","#2171B5","#08306B")) +
    labs(title=g, x="UMAP 1", y="UMAP 2") + tpub + theme(legend.title=element_text(size=8))
}

p_FLT1 <- gene_umaps("FLT1"); p_FN1 <- gene_umaps("FN1")
p_HLA <- gene_umaps("HLA-DRA"); p_CEBPA <- gene_umaps("CEBPA")
p_PAPPA <- gene_umaps("PAPPA"); p_SPP1 <- gene_umaps("SPP1")

# Save
ggsave(file.path(FIGDIR,"PE_UMAP_groups.png"), pA, w=9, h=7, dpi=300, bg="white")
ggsave(file.path(FIGDIR,"PE_UMAP_subtypes.png"), pB, w=9, h=7, dpi=300, bg="white")
if(!is.null(p_FLT1)) ggsave(file.path(FIGDIR,"PE_UMAP_FLT1.png"), p_FLT1, w=7, h=6, dpi=300, bg="white")
if(!is.null(p_FN1)) ggsave(file.path(FIGDIR,"PE_UMAP_FN1.png"), p_FN1, w=7, h=6, dpi=300, bg="white")
if(!is.null(p_HLA)) ggsave(file.path(FIGDIR,"PE_UMAP_HLA-DRA.png"), p_HLA, w=7, h=6, dpi=300, bg="white")
if(!is.null(p_CEBPA)) ggsave(file.path(FIGDIR,"PE_UMAP_CEBPA.png"), p_CEBPA, w=7, h=6, dpi=300, bg="white")
if(!is.null(p_PAPPA)) ggsave(file.path(FIGDIR,"PE_UMAP_PAPPA.png"), p_PAPPA, w=7, h=6, dpi=300, bg="white")
if(!is.null(p_SPP1)) ggsave(file.path(FIGDIR,"PE_UMAP_SPP1.png"), p_SPP1, w=7, h=6, dpi=300, bg="white")

cat("All UMAPs saved\n")
