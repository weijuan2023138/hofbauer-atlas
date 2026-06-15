#!/usr/bin/env Rscript
# Fig1G: Developmental gene dotplot (Early vs Late) — adapted for new Atlas
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

# Gene rows
row1 <- c("TREM2","AXL","CEBPA","ID2","CD5L","NOTCH2")
row2 <- c("TIMP1","VIM","MMP14","ENO1","PGK1","COL1A2","FN1","PDGFB","SOD2")
row3 <- c("CCL8","IL18","IFI30","HLA-DRA","CTSS","FCGR3A")

all_genes <- c(row1, row2, row3)
all_genes <- intersect(all_genes, rownames(seu))

expr <- FetchData(seu, vars=c(all_genes, "tri"))
plot_data <- data.frame()
for(gene in all_genes) {
  for(tr in levels(seu$tri)) {
    vals <- expr[expr$tri==tr, gene]
    plot_data <- rbind(plot_data, data.frame(
      Gene=gene, Trimester=tr,
      Mean=mean(vals), SEM=sd(vals)/sqrt(length(vals))
    ))
  }
}

# Assign category
plot_data$Category <- "Progenitor & Diff"
plot_data$Category[plot_data$Gene %in% row2] <- "ECM & Metabolism"
plot_data$Category[plot_data$Gene %in% row3] <- "Immune Activation"
plot_data$Category <- factor(plot_data$Category, levels=c("Progenitor & Diff","ECM & Metabolism","Immune Activation"))

tri_cols <- c("Early"="#4575B4","Late"="#D73027")

p <- ggplot(plot_data, aes(Gene, Mean, fill=Trimester)) +
  geom_bar(stat="identity", position=position_dodge(0.7), width=0.6) +
  geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), position=position_dodge(0.7), width=0.2, linewidth=0.3) +
  scale_fill_manual(values=tri_cols) +
  facet_wrap(~Category, scales="free_x", nrow=3) +
  labs(y="Mean Expression", x="") +
  theme_bw(11) + theme(panel.grid=element_blank(), legend.position="top",
    legend.title=element_blank(), legend.text=element_text(size=9),
    strip.text=element_text(size=10,face="bold"),
    strip.background=element_rect(fill="white"),
    axis.text.x=element_text(angle=45, hjust=1, size=9))
ggsave(file.path(FIGDIR,"Fig1G_dev_dotplot_genes.png"), p, w=10, h=8, dpi=300, bg="white")
cat("Saved Fig1G_dev_dotplot_genes.png\n")
