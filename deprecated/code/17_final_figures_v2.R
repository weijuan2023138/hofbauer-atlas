#!/usr/bin/env Rscript
# Final figures: subtype UMAP, dataset UMAP, DotPlot, classifier genes
library(Seurat); library(ggplot2); library(patchwork); library(dplyr)

seu <- readRDS("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/Hofbauer_Atlas_Final.rds")
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures"
cat(ncol(seu),"cells\n")

# ============================================================
# Colors
# ============================================================
st_cols <- c("Pro-inflammatory"="#C62828","MHCII+/CCL13+C1Q+"="#E65100",
  "Homeostatic/SPP1+"="#1565C0","PRKN+ Autophagy"="#6A1B9A",
  "Vascular remodeling"="#2E7D32","MKI67+ Proliferating"="#455A64")

ds_cols <- c("Arutyunyan"="#D73027","GSE290578"="#4575B4","gse214607"="#4DAF4A",
  "hoo_2024"="#984EA3","gse173193"="#FF7F00","gse183338"="#A65628",
  "gse298119"="#F781BF","my_preterm_cohort"="#878787","UCSF_Li_2026"="#66C2A5")

dis_cols <- c("Normal 1st"="#4DBBD5","Normal 1st/2nd"="#00A087","Late pregnancy"="#7E6148",
  "PE"="#C62828","PTNL"="#BCAAA4","PTL"="#E18727","TL"="#3C5488",
  "Miscarriage"="#F39B7F","Infection"="#DC0000")

# ============================================================
# Theme
# ============================================================
tpub <- theme_bw(11) +
  theme(panel.grid=element_line(color="grey92",linewidth=0.2),
        panel.border=element_rect(color="black",fill=NA,linewidth=0.5),
        axis.text=element_blank(), axis.ticks=element_blank(),
        legend.position="right", legend.title=element_blank(),
        legend.text=element_text(size=9), legend.key=element_blank(),
        plot.title=element_text(size=13,face="bold"))

set.seed(42); idx <- sample(ncol(seu))

# ============================================================
# FIGURE 1: Subtype UMAP
# ============================================================
p1 <- ggplot(seu@meta.data[idx,], aes(UMAP_1,UMAP_2,color=subtype)) +
  geom_point(size=0.15,alpha=0.8) + scale_color_manual(values=st_cols) +
  labs(title="Hofbauer Atlas (n=17,896)",x="UMAP 1",y="UMAP 2") +
  tpub + guides(color=guide_legend(override.aes=list(size=3)))
ggsave(file.path(FIGDIR,"FigA_subtype.png"),p1,w=9,h=7,dpi=300,bg="white")
cat("Saved FigA\n")

# ============================================================
# FIGURE 2: Dataset + Disease UMAP
# ============================================================
# Add metadata safely
ds <- factor(seu$dataset, levels=names(ds_cols))
names(ds) <- colnames(seu); seu$ds <- ds

short <- c("Normal 1st trimester"="Normal 1st","Normal 1st/2nd/Term"="Normal 1st/2nd",
  "Normal 3rd trimester / Preeclampsia"="Late pregnancy","Preeclampsia"="PE",
  "Preterm No Labor"="PTNL","Preterm Labor"="PTL","Term Labor"="TL",
  "Miscarriage / Normal"="Miscarriage","Infection"="Infection")
dis <- factor(short[seu$disease_group], levels=c("Normal 1st","Normal 1st/2nd",
  "Late pregnancy","PE","PTNL","PTL","TL","Miscarriage","Infection"))
names(dis) <- colnames(seu); seu$dis <- dis

p2a <- ggplot(seu@meta.data[idx,], aes(UMAP_1,UMAP_2,color=ds)) +
  geom_point(size=0.12,alpha=0.8) + scale_color_manual(values=ds_cols) +
  labs(title="By Dataset",x="UMAP 1",y="UMAP 2") +
  tpub + guides(color=guide_legend(override.aes=list(size=3),ncol=1))

p2b <- ggplot(seu@meta.data[idx,], aes(UMAP_1,UMAP_2,color=dis)) +
  geom_point(size=0.12,alpha=0.8) + scale_color_manual(values=dis_cols) +
  labs(title="By Disease Group",x="UMAP 1",y="UMAP 2") +
  tpub + guides(color=guide_legend(override.aes=list(size=3),ncol=1))

ggsave(file.path(FIGDIR,"FigB_dataset_disease.png"),p2a+p2b,w=18,h=7,dpi=300,bg="white")
cat("Saved FigB\n")

# ============================================================
# FIGURE 3: DotPlot
# ============================================================
DefaultAssay(seu) <- "RNA"; Idents(seu) <- "subtype"

markers <- c("FOLR2","CD163","MRC1","LYVE1","CD36","TREM2","MAF","DAB2",
  "IL1B","TNF","CXCL8","NFKBIZ","CCL13","AIF1",
  "HLA-DRA","HLA-DRB1","HLA-DPA1","C1QA","C1QB","FCGR3A",
  "PRKN","C9","SOX5","MKI67","BUB1B","KIF4A",
  "SPP1","FN1","PAPPA","FLT1")

present <- markers[markers %in% rownames(seu)]
cat(length(present),"markers present\n")

p3 <- DotPlot(seu, features=present, assay="RNA") + RotatedAxis() +
  scale_color_gradientn(colors=c("lightgrey","#1565C0","#C62828")) +
  labs(title="Key Marker Genes") +
  theme_bw(11) + theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
    panel.grid=element_line(color="grey92",linewidth=0.2), legend.position="right")

ggsave(file.path(FIGDIR,"FigC_dotplot.png"),p3,w=20,h=6,dpi=300,bg="white")
cat("Saved FigC\n")

# ============================================================
# Classifier genes
# ============================================================
clf <- read.csv("/home/weijuan/文档/胎盘单细胞数据/results/phase1_classifier/classifier_genes.csv")
write.csv(clf, file.path(FIGDIR,"../results/classifier_genes.csv"), row.names=FALSE)
cat("Saved classifier_genes (",nrow(clf),"genes)\n")

cat("\nDone: FigA_subtype, FigB_dataset_disease, FigC_dotplot, classifier_genes.csv\n")
