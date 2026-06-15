#!/usr/bin/env Rscript
# 补充FigC dotplot — exact replica of old format + bold black text
library(Seurat); library(ggplot2)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIG   <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures/补充FigC_dotplot.png"

seu <- readRDS(INPUT)
DefaultAssay(seu) <- "RNA"
Idents(seu) <- "subtype_pred"

classifier_markers <- c(
  "LYVE1","CD36","DAB2","VSIG4","F13A1","SPP1",
  "FSCN1","COLEC12","ABCG2","CD5L",
  "HLA-DRB5","HLA-DRB1","HLA-DPA1","CST7","HLA-DRA",
  "FOLR2","CD163","MRC1","TREM2","C1QA","C1QB",
  "IL1B","TNF","CXCL8","CCL13","PRKN","C9",
  "MKI67","BUB1B","FN1","PAPPA","FLT1"
)
present <- classifier_markers[classifier_markers %in% rownames(seu)]

p <- DotPlot(seu, features=present, assay="RNA") +
  RotatedAxis() +
  scale_color_gradientn(colors=c("lightgrey","#1B6B93","#C62828")) +
  labs(title="Key Marker Genes by Subtype") +
  theme_bw(11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8, color="black", face="bold"),
        axis.text.y=element_text(color="black", face="bold"),
        panel.grid=element_line(color="grey92", linewidth=0.2),
        legend.position="right")

ggsave(FIG, p, width=16, height=5, dpi=300, bg="white")
cat("Saved: 补充FigC_dotplot.png\n")
