#!/usr/bin/env Rscript
# Fig1G: Dev dotplot — exact match to old 25_dotplot_genes.R format
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

row1 <- c("TREM2","AXL","CEBPA","ID2","CD5L","NOTCH2")
row2 <- c("TIMP1","VIM","MMP14","ENO1","PGK1","COL1A2","FN1","PDGFB","SOD2")
row3 <- c("CCL8","IL18","IFI30","HLA-DRA","CTSS","FCGR3A")

all_genes <- intersect(c(row1, row2, row3), rownames(seu))
expr <- FetchData(seu, vars=c(all_genes, "tri"))
plot_data <- data.frame()
for(gene in all_genes) {
  for(tr in levels(seu$tri)) {
    vals <- expr[expr$tri==tr, gene]
    plot_data <- rbind(plot_data, data.frame(Gene=gene, Trimester=tr,
      Mean=mean(vals), SEM=sd(vals)/sqrt(length(vals))))
  }
}
plot_data$Trimester <- factor(plot_data$Trimester, levels=c("Early","Late"))
plot_data$Module <- ifelse(plot_data$Gene %in% row1, "Development",
                    ifelse(plot_data$Gene %in% row2, "Remodeling", "Immunity"))
plot_data$Module <- factor(plot_data$Module, levels=c("Development","Remodeling","Immunity"))
plot_data$Gene <- factor(plot_data$Gene, levels=all_genes)
tri_cols <- c("Early"="#4575B4","Late"="#D73027")

p <- ggplot(plot_data, aes(x=Trimester, y=Mean, color=Trimester)) +
  geom_line(aes(group=Gene), color="grey60", linewidth=0.4) +
  geom_point(size=2.8) +
  geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=0.1, linewidth=0.5) +
  scale_color_manual(values=tri_cols) +
  facet_grid(Module ~ Gene, scales="free_y", switch="y") +
  labs(y="Mean Expression ± SEM", x="", subtitle="Early (GW4.5–10) vs Late (GW32–38)") +
  theme_bw() + theme(axis.text=element_text(color="black",size=8),
    axis.text.x=element_text(size=9,color="black"), axis.title.y=element_text(size=10),
    legend.position="top", legend.title=element_blank(),
    legend.key.size=unit(0.4,"cm"), legend.text=element_text(size=10),
    strip.text.y=element_text(size=10,face="bold"),
    strip.text.x=element_text(size=9,face="bold.italic"),
    strip.background=element_rect(fill="grey95",color="black",linewidth=0.4),
    panel.grid.major=element_line(color="grey92",linewidth=0.3), panel.grid.minor=element_blank(),
    panel.border=element_rect(color="black",linewidth=0.4), panel.spacing=unit(0.4,"lines"),
    plot.subtitle=element_text(hjust=0.5,size=10,face="plain"))
w <- max(11, length(all_genes)*0.9)
ggsave(file.path(FIGDIR,"Fig1G_dev_dotplot_genes.png"), p, w=w, h=6, dpi=300, bg="white")

for(gene in all_genes) {
  means <- sapply(levels(seu$tri), function(tr) mean(expr[expr$tri==tr, gene]))
  fc <- means[2]/means[1]
  dir <- if(fc>1.3) "↑↑" else if(fc>1.1) "↑" else if(fc<0.7) "↓↓" else if(fc<0.9) "↓" else "→"
  cat(sprintf("%-10s  E=%.3f  L=%.3f  FC=%.2f  %s\n", gene, means[1], means[2], fc, dir))
}
cat(sprintf("\n%d genes | %d cells\n", length(all_genes), ncol(seu)))
