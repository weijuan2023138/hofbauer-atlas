#!/usr/bin/env Rscript
# Supplementary dot plot: 3 merged modules, no overlap with main figure
library(Seurat); library(ggplot2); library(dplyr)

seu <- readRDS('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/Hofbauer_Atlas_Final.rds')
cond <- read.csv("/tmp/gse290578_conditions.csv", stringsAsFactors=FALSE)
rownames(cond) <- cond$barcode; seu$gse <- NA
for(bc in colnames(seu)) if(bc %in% cond$barcode) seu$gse[bc] <- cond[bc,"condition"]
seu$tri <- NA; seu$tri[seu$dataset=="Arutyunyan"] <- "Early"
seu$tri[seu$dataset=="GSE290578" & seu$gse=="Normal"] <- "Late"
seu_normal <- subset(seu, cells=colnames(seu)[!is.na(seu$tri)])

row1 <- c("MKI67","TOP2A","CCNA2","PCNA","BIRC5","CDK1","AURKB","CDC20")
row2 <- c("LDHA","HK2","GAPDH","PKM","ALDOA","PFKL","PGAM1","PRDX1","CAT")
row3 <- c("HLA-DRB1","CD74","HLA-DMA","HLA-DMB","C1QA","C1QB","C1QC","C3","C4A","CFB","CFH")

all_genes <- c(row1, row2, row3)
all_genes <- intersect(all_genes, rownames(seu_normal))

expr <- FetchData(seu_normal, vars=c(all_genes, "tri"))
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
                    ifelse(plot_data$Gene %in% row2, "Metabolism & Stress",
                           "Immune Effector"))
plot_data$Module <- factor(plot_data$Module, levels=c("Proliferation","Metabolism & Stress","Immune Effector"))
plot_data$Gene <- factor(plot_data$Gene, levels=all_genes)

tri_cols <- c("Early"="#4575B4","Late"="#D73027")

p <- ggplot(plot_data, aes(x=Trimester, y=Mean, color=Trimester)) +
  geom_line(aes(group=Gene), color="grey75", linewidth=0.3) +
  geom_point(size=2.2) +
  geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=0.12, linewidth=0.4) +
  scale_color_manual(values=tri_cols) +
  facet_grid(Module ~ Gene, scales="free_y", switch="y") +
  labs(y="Mean Expression \u00b1 SEM", x="",
       subtitle="Early (GW4.5\u201310) vs Late (GW32\u201338)") +
  theme_bw() +
  theme(
    axis.text=element_text(color="black", size=8),
    axis.text.x=element_text(size=9),
    axis.title.y=element_text(size=10),
    legend.position="top", legend.title=element_blank(),
    legend.key.size=unit(0.4,"cm"), legend.text=element_text(size=10),
    strip.text.y=element_text(size=9, face="bold"),
    strip.text.x=element_text(size=8, face="bold.italic"),
    strip.background=element_rect(fill="grey95", color="black", linewidth=0.4),
    panel.grid.major=element_line(color="grey92", linewidth=0.3),
    panel.grid.minor=element_blank(),
    panel.border=element_rect(color="black", linewidth=0.4),
    panel.spacing=unit(0.3,"lines"),
    plot.subtitle=element_text(hjust=0.5, size=10)
  )

FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'
ggsave(file.path(FIGDIR,"FigS_dotplot_extra.png"), p, w=20, h=6.5, dpi=300, bg="white")
cat(sprintf("Saved: %d genes in 3 modules\n", length(all_genes)))
