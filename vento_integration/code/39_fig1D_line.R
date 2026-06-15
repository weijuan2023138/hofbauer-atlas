#!/usr/bin/env Rscript
# Fig1D: Subtype composition line chart (Early vs Late, Normal only)
library(Seurat); library(ggplot2); library(dplyr)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
seu$disease_final <- factor(labels$disease_final[1:ncol(seu)], levels=c("Normal","PE","Miscarriage","Infection","Preterm"))
seu <- subset(seu, disease_final=="Normal")

tri_map <- c("E-MTAB-12421"="Early","E-MTAB-6701"="Early","E-MTAB-12795"="Early","GSE214607"="Early",
             "UCSF Li 2026"="Mid","GSE290578"="Late","GSE333257"="Late",
             "GSE298602"="Late","GSE298119"="Late","GSE173193"="Late")
seu$tri <- setNames(tri_map[as.character(seu$dataset)], colnames(seu))
seu <- subset(seu, tri %in% c("Early","Late"))
seu$tri <- droplevels(seu$tri)

prop <- prop.table(table(seu$subtype_pred, seu$tri), margin=2)
df <- as.data.frame(prop); colnames(df) <- c("Subtype","Trimester","Proportion")

sc <- c("Pro-inflammatory"="#C62828","MHCII+ Antigen-presenting"="#E65100",
        "Homeostatic"="#1565C0","PRKN+ Autophagy"="#6A1B9A",
        "Vascular remodeling"="#2E7D32","MKI67+ Proliferating"="#455A64")

p <- ggplot(df, aes(Trimester, Proportion*100, color=Subtype, group=Subtype)) +
  geom_line(linewidth=1.2) + geom_point(size=3) +
  scale_color_manual(values=sc) +
  labs(y="Proportion (%)", x="", title="Subtype Composition Across Gestation") +
  theme_bw(11) + theme(panel.grid=element_line(color="grey92",linewidth=0.2),
    legend.title=element_blank(), legend.text=element_text(size=9),
    axis.text=element_text(size=10), axis.title=element_text(size=11),
    plot.title=element_text(size=12,face="bold",hjust=0.5))

ggsave(file.path(FIGDIR,"Fig1D_subtype_line.png"), p, w=6, h=5, dpi=300, bg="white")
cat("Saved\n")
