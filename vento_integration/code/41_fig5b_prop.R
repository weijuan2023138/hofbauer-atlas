#!/usr/bin/env Rscript
# Fig5b: Subtype proportions bar chart (6 disease groups, new Atlas)
library(Seurat); library(ggplot2); library(dplyr)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
seu$disease_final <- factor(labels$disease_final[1:ncol(seu)], 
  levels=c("Normal","PE","Miscarriage","Infection","Preterm"))

tri_map <- c("E-MTAB-12421"="Early","E-MTAB-6701"="Early","E-MTAB-12795"="Early",
             "GSE214607"="Early","UCSF Li 2026"="Mid","GSE290578"="Late",
             "GSE333257"="Late","GSE298602"="Late","GSE298119"="Late","GSE173193"="Late")
seu$trimester <- setNames(tri_map[as.character(seu$dataset)], colnames(seu))

seu$disease_clean <- as.character(seu$disease_final)
seu$disease_clean[seu$disease_final=="Normal" & seu$trimester=="Early"] <- "Normal_Early"
seu$disease_clean[seu$disease_final=="Normal" & seu$trimester=="Late"] <- "Normal_Late"
keep <- seu$disease_clean %in% c("Normal_Early","Normal_Late","Miscarriage","PE","Infection","Preterm")
seu <- subset(seu, cells=colnames(seu)[keep])
seu$disease_clean <- factor(seu$disease_clean, 
  levels=c("Normal_Early","Miscarriage","Infection","Normal_Late","PE","Preterm"))

st_cols <- c("Pro-inflammatory"="#C62828","MHCII+ Antigen-presenting"="#E65100",
             "Homeostatic"="#1565C0","PRKN+ Autophagy"="#6A1B9A",
             "Vascular remodeling"="#2E7D32","MKI67+ Proliferating"="#455A64")

prop <- prop.table(table(seu$subtype_pred, seu$disease_clean), margin=2)
prop_df <- as.data.frame(prop); colnames(prop_df) <- c("Subtype","Disease","Proportion")

p <- ggplot(prop_df, aes(x=Disease, y=Proportion*100, fill=Subtype)) +
  geom_col(width=0.65) + scale_fill_manual(values=st_cols) +
  labs(title="Subtype proportions", y="Proportion (%)", x="") +
  theme_bw(base_size=11) +
  theme(panel.grid=element_blank(), plot.title=element_text(face="bold",size=13,hjust=0.5,color="black"),
        legend.position="right",
        axis.text.x=element_text(face="italic",size=9,color="black"),
        axis.text.y=element_text(size=9,color="black"),
        axis.title.y=element_text(size=10,color="black"),
        legend.text=element_text(size=9,color="black"),
        legend.title=element_text(size=10,face="bold",color="black"))
ggsave(file.path(FIGDIR,"Fig5b_subtype_proportions.png"), p, w=7, h=5, dpi=300, bg="white")
cat("Fig5b done\n")
