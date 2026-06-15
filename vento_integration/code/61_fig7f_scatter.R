#!/usr/bin/env Rscript
# Fig7f: ECM vs Immune scatter by PE subtype
library(Seurat); library(ggplot2); library(dplyr)
set.seed(42)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
detail <- labels$disease_detail[1:ncol(seu)]

dc <- rep(NA_character_, ncol(seu))
dc[seu$dataset=="GSE290578" & detail=="PE"] <- "Early PE"
dc[detail %in% c("PreE_SF","gHTN","GSE173193","GSE298119")] <- "Late PE"
dc[(seu$dataset=="GSE290578" & detail=="Normal") | (seu$dataset=="GSE298602" & detail=="Control")] <- "Normal Late"
seu <- subset(seu, cells=colnames(seu)[dc %in% c("Normal Late","Early PE","Late PE")])
seu$group <- factor(dc[dc %in% c("Normal Late","Early PE","Late PE")], levels=c("Normal Late","Early PE","Late PE"))
seu <- JoinLayers(seu)

ecm <- rowMeans(FetchData(seu, vars=intersect(c("FN1","SPP1","COL1A1","MMP9"), rownames(seu))))
imm <- rowMeans(FetchData(seu, vars=intersect(c("IL1B","TNF","CXCL8","CD44","CD47"), rownames(seu))))
st_cols <- c("Pro-inflammatory"="#C62828","MHCII+ Antigen-presenting"="#E65100",
  "Homeostatic"="#1565C0","PRKN+ Autophagy"="#6A1B9A",
  "Vascular remodeling"="#2E7D32","MKI67+ Proliferating"="#455A64")

df <- data.frame(ECM=ecm, Immune=imm, Subtype=seu$subtype_pred, Group=seu$group)
df <- df[df$Group!="Normal Late",]
df_plot <- do.call(rbind, lapply(split(df, list(df$Group,df$Subtype)), function(d) d[sample(nrow(d), min(200, nrow(d))),]))

p <- ggplot(df_plot, aes(x=ECM, y=Immune, color=Subtype)) +
  geom_point(size=1, alpha=0.6) + scale_color_manual(values=st_cols) +
  facet_wrap(~Group, scales="free") +
  labs(title="ECM vs Immune module — PE subtypes", x="ECM module score", y="Immune module score") +
  theme_bw(12) + theme(panel.grid.major=element_line(color="grey92", linewidth=0.3),
    panel.grid.minor=element_blank(), strip.text=element_text(face="bold", size=14, color="black"),
    axis.title=element_text(face="bold", size=12), legend.text=element_text(size=9),
    plot.title=element_text(face="bold", size=15, hjust=0.5))

ggsave(file.path(FIGDIR, "Fig7f_ECM_Immune_scatter.png"), p, w=10, h=5.5, dpi=300, bg="white")
cat("Fig7f saved\n")
