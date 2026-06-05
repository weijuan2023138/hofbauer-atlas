#!/usr/bin/env Rscript
# Violin plot: developmental gene expression across trimesters
library(Seurat); library(ggplot2); library(dplyr); library(patchwork)

seu <- readRDS('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/Hofbauer_Atlas_Final.rds')
cond <- read.csv("/tmp/gse290578_conditions.csv", stringsAsFactors=FALSE)
rownames(cond) <- cond$barcode
seu$gse <- NA
for(bc in colnames(seu)) if(bc %in% cond$barcode) seu$gse[bc] <- cond[bc,"condition"]

seu$tri <- NA
seu$tri[seu$dataset=="Arutyunyan"] <- "Early"
seu$tri[seu$dataset=="UCSF_Li_2026"] <- "Mid"
seu$tri[seu$dataset=="GSE290578" & seu$gse=="Normal"] <- "Late"
seu$tri <- factor(seu$tri, levels=c("Early","Mid","Late"))
seu_normal <- subset(seu, cells=colnames(seu)[!is.na(seu$tri)])
Idents(seu_normal) <- "tri"

genes <- c(
  # Differentiation & Development
  "TREM2","AXL","CEBPA","ID2","CD36","TIMP1","FOLR2",
  # Tissue Remodeling & Signaling
  "RHOA","ITGAV","TGFB1","VIM","ANGPT2",
  # Immune Maturation
  "CD74","HLA-DRB1","S100A8","S100A9","CCL2","NFKB1","CXCL8"
)
genes <- intersect(genes, rownames(seu_normal))
cat(sprintf("%d genes\n", length(genes)))

tri_cols <- c("Early"="#4575B4","Mid"="#FDAE61","Late"="#D73027")

# Build individual violin plots with shared legend
plots <- list()
for(gene in genes) {
  expr <- FetchData(seu_normal, vars=c(gene, "tri"))
  colnames(expr)[1] <- "Expression"
  
  p <- ggplot(expr, aes(x=tri, y=Expression, fill=tri)) +
    geom_violin(scale="width", trim=TRUE, alpha=0.8, color=NA, linewidth=0) +
    geom_boxplot(width=0.2, outlier.size=0.15, alpha=0.5, linewidth=0.2,
                 color="grey30", fill="white") +
    scale_fill_manual(values=tri_cols, guide="none") +
    labs(title=gene, x="", y="") +
    theme_classic() +
    theme(
      plot.title=element_text(size=10, face="bold.italic", hjust=0.5),
      axis.text.x=element_text(size=8, color="black"),
      axis.text.y=element_text(size=7, color="black"),
      axis.line=element_line(linewidth=0.3),
      axis.ticks=element_line(linewidth=0.3),
      plot.margin=margin(2,2,2,2)
    )
  plots[[gene]] <- p
}

# Arrange in grid: 3 columns, rows vary by module
mod1 <- wrap_plots(plots[genes[genes %in% c("TREM2","AXL","CEBPA","ID2","CD36","TIMP1","FOLR2")]], ncol=4) +
  plot_annotation(title="Differentiation & Development", 
    theme=theme(plot.title=element_text(size=11, face="bold", hjust=0.5)))
mod2 <- wrap_plots(plots[genes[genes %in% c("RHOA","ITGAV","TGFB1","VIM","ANGPT2")]], ncol=5) +
  plot_annotation(title="Tissue Remodeling & Signaling",
    theme=theme(plot.title=element_text(size=11, face="bold", hjust=0.5)))
mod3 <- wrap_plots(plots[genes[genes %in% c("CD74","HLA-DRB1","S100A8","S100A9","CCL2","NFKB1","CXCL8")]], ncol=4) +
  plot_annotation(title="Immune Maturation",
    theme=theme(plot.title=element_text(size=11, face="bold", hjust=0.5)))

final <- mod1 / mod2 / mod3

FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'
ggsave(file.path(FIGDIR,"Fig_dev_violin.png"), final, w=14, h=10, dpi=300, bg="white")
cat("\nSaved Fig_dev_violin.png\n")
