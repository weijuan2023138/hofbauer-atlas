#!/usr/bin/env Rscript
# Fig1X: Protein marker validation — feature plots of canonical Hofbauer markers
library(Seurat); library(ggplot2); library(patchwork)

seu <- readRDS("results/Hofbauer_Atlas_Final.rds")
cond <- read.csv("/tmp/gse290578_conditions.csv", stringsAsFactors=FALSE)
rownames(cond) <- cond$barcode; seu$gse <- NA
for(bc in colnames(seu)) if(bc %in% cond$barcode) seu$gse[bc] <- cond[bc,"condition"]
seu$tri <- NA; seu$tri[seu$dataset=="Arutyunyan"] <- "Early"
seu$tri[seu$dataset=="UCSF_Li_2026"] <- "Mid"
seu$tri[seu$dataset=="GSE290578" & seu$gse=="Normal"] <- "Late"
seu_normal <- subset(seu, cells=colnames(seu)[!is.na(seu$tri)])
DefaultAssay(seu_normal) <- "RNA"

markers <- c("FOLR2","CD163","MRC1","CD68","FCGR3A","HLA-DRA","TREM2","C1QA")
markers <- intersect(markers, rownames(seu_normal))

# Feature plots — 2 rows of 4
fp_list <- list()
for(g in markers) {
  fp_list[[g]] <- FeaturePlot(seu_normal, features=g, pt.size=0.3, max.cutoff="q99",
    order=TRUE, cols=c("grey92","#440154","#FDE725")) +
    ggtitle(g) + NoLegend() + NoAxes() +
    theme(plot.title=element_text(size=11, face="bold", hjust=0.5))
}

png("figures/Fig1/Fig1X_protein_validation.png", w=16, h=4.5, units="in", res=300)
wrap_plots(fp_list, nrow=1) +
  plot_annotation(title="Hofbauer Cell Canonical Protein Markers",
    theme=theme(plot.title=element_text(hjust=0.5, face="bold", size=13)))
dev.off()
cat("Protein validation saved\n")
