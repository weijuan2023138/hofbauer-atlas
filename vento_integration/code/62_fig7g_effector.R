#!/usr/bin/env Rscript
# Fig7g: Effector gene expression bar chart — PE subtypes
library(Seurat); library(ggplot2); library(dplyr)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
detail <- labels$disease_detail[1:ncol(seu)]

dc <- rep(NA_character_, ncol(seu))
dc[(seu$dataset=="GSE290578" & detail=="Normal") | (seu$dataset=="GSE298602" & detail=="Control")] <- "Normal_Late"
dc[seu$dataset=="GSE290578" & detail=="PE"] <- "Early_PE"
dc[detail %in% c("PreE_SF","gHTN","GSE173193","GSE298119")] <- "Late_PE"
seu <- subset(seu, cells=colnames(seu)[dc %in% c("Normal_Late","Early_PE","Late_PE")])
seu$group <- factor(dc[dc %in% c("Normal_Late","Early_PE","Late_PE")], levels=c("Normal_Late","Early_PE","Late_PE"))
seu <- JoinLayers(seu)

genes <- c("FLT1","FN1","SPP1","COL1A2","HLA-DRA","CD74","STAT3")
genes <- intersect(genes, rownames(seu))
grp_cols <- c(Normal_Late="#66C2A5",Early_PE="#FC8D59",Late_PE="#8DA0CB")

df_plot <- data.frame()
for(g in genes) for(grp in levels(seu$group)) {
  vals <- FetchData(seu, g)[seu$group==grp, 1]
  df_plot <- rbind(df_plot, data.frame(Gene=g, Group=grp, Mean=mean(vals), SE=sd(vals)/sqrt(length(vals))))
}
df_plot$Gene <- factor(df_plot$Gene, levels=genes)
df_plot$Group <- factor(df_plot$Group, levels=c("Normal_Late","Early_PE","Late_PE"))

p <- ggplot(df_plot, aes(x=Gene, y=Mean, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(0.8), width=0.7, color="black", linewidth=0.3) +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), position=position_dodge(0.8), width=0.2, linewidth=0.3) +
  scale_fill_manual(values=grp_cols) +
  labs(title="Effector gene expression — PE subtypes", y="Mean expression", x="") +
  theme_bw(14) + theme(axis.text.x=element_text(size=14, color="black", face="bold.italic"),
    axis.text.y=element_text(size=12, color="black"), axis.title.y=element_text(size=13, face="bold"),
    panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
    plot.title=element_text(hjust=0.5, size=16, face="bold"),
    legend.position="top", legend.title=element_blank(), legend.text=element_text(size=13))

ggsave(file.path(FIGDIR, "Fig7g_effector_bars.png"), p, w=10, h=5.5, dpi=300, bg="white")
cat("Fig7g saved\n")
