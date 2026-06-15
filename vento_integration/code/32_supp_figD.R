#!/usr/bin/env Rscript
# 补充FigD: Extra dotplot (Proliferation / Metabolism / Immune Effector)
library(Seurat); library(ggplot2); library(dplyr)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
seu$disease_final <- factor(labels$disease_final[1:ncol(seu)], 
  levels=c("Normal","PE","Miscarriage","Infection","Preterm"))
seu <- subset(seu, disease_final=="Normal")

tri_map <- c("E-MTAB-12421"="Early","E-MTAB-6701"="Early","E-MTAB-12795"="Early",
             "GSE214607"="Early","UCSF Li 2026"="Mid","GSE290578"="Late",
             "GSE333257"="Late","GSE298602"="Late","GSE298119"="Late","GSE173193"="Late")
seu$tri <- setNames(tri_map[as.character(seu$dataset)], colnames(seu))
seu <- subset(seu, tri %in% c("Early","Late"))
seu$tri <- factor(seu$tri, levels=c("Early","Late"))

row1 <- c("MKI67","TOP2A","CCNA2","PCNA","BIRC5","CDK1","AURKB","CDC20")
row2 <- c("LDHA","HK2","GAPDH","PKM","ALDOA","PFKL","PGAM1","PRDX1","CAT")
row3 <- c("HLA-DRB1","CD74","HLA-DMA","HLA-DMB","C1QA","C1QB","C1QC","C3","C4A","CFB","CFH")

all_genes <- intersect(c(row1, row2, row3), rownames(seu))
expr <- FetchData(seu, vars=c(all_genes, "tri"))
plot_data <- data.frame()
for(gene in all_genes) {
  for(tr in c("Early","Late")) {
    vals <- expr[expr$tri==tr, gene]
    plot_data <- rbind(plot_data, data.frame(Gene=gene, Trimester=tr,
      Mean=mean(vals), SEM=sd(vals)/sqrt(length(vals))))
  }
}
plot_data$Trimester <- factor(plot_data$Trimester, levels=c("Early","Late"))
plot_data$Module <- ifelse(plot_data$Gene %in% row1, "Proliferation",
                    ifelse(plot_data$Gene %in% row2, "Metabolism & Stress", "Immune Effector"))
plot_data$Module <- factor(plot_data$Module, levels=c("Proliferation","Metabolism & Stress","Immune Effector"))
plot_data$Gene <- factor(plot_data$Gene, levels=all_genes)
tri_cols <- c("Early"="#4575B4","Late"="#D73027")

p <- ggplot(plot_data, aes(x=Trimester, y=Mean, color=Trimester)) +
  geom_line(aes(group=Gene), color="grey75", linewidth=0.3) +
  geom_point(size=2.2) + geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=0.12, linewidth=0.4) +
  scale_color_manual(values=tri_cols) +
  facet_grid(Module ~ Gene, scales="free_y", switch="y") +
  labs(y="Mean Expression ± SEM", x="", subtitle="Early (GW4.5–10) vs Late (GW32–38)") +
  theme_bw() + theme(axis.text=element_text(color="black",size=8), axis.text.x=element_text(size=9),
    axis.title.y=element_text(size=10), legend.position="top", legend.title=element_blank(),
    legend.key.size=unit(0.4,"cm"), legend.text=element_text(size=10),
    strip.text.y=element_text(size=9,face="bold"), strip.text.x=element_text(size=8,face="bold.italic"),
    strip.background=element_rect(fill="grey95",color="black",linewidth=0.4),
    panel.grid.major=element_line(color="grey92",linewidth=0.3), panel.grid.minor=element_blank(),
    panel.border=element_rect(color="black",linewidth=0.4), panel.spacing=unit(0.3,"lines"),
    plot.subtitle=element_text(hjust=0.5,size=10))
ggsave(file.path(FIGDIR,"补充FigD_dotplot_extra.png"), p, w=20, h=6.5, dpi=300, bg="white")
cat(sprintf("Saved: %d genes in 3 modules\n", length(all_genes)))
