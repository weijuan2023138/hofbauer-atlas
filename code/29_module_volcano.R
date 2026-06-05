#!/usr/bin/env Rscript
# Supplementary: Original 8 modules, genes added (no module renaming/merging)
library(Seurat); library(ggplot2); library(dplyr); library(patchwork)

seu <- readRDS('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/Hofbauer_Atlas_Final.rds')
cond <- read.csv("/tmp/gse290578_conditions.csv", stringsAsFactors=FALSE)
rownames(cond) <- cond$barcode; seu$gse <- NA
for(bc in colnames(seu)) if(bc %in% cond$barcode) seu$gse[bc] <- cond[bc,"condition"]
seu$tri <- NA; seu$tri[seu$dataset=="Arutyunyan"] <- "Early"
seu$tri[seu$dataset=="GSE290578" & seu$gse=="Normal"] <- "Late"
seu_normal <- subset(seu, cells=colnames(seu)[!is.na(seu$tri)])

# Original 8 modules + expanded genes
mods <- list(
  "Progenitor & Diff" = c("TREM2","AXL","CEBPA","ID2","CD5L","NOTCH2","SOX4","MAFB","HMGA2"),
  "Proliferation"    = c("MKI67","TOP2A","CCNA2","PCNA","BIRC5","CDK1","AURKB","CDC20","CCNB1"),
  "Glycolysis"       = c("ENO1","PGK1","LDHA","HK2","GAPDH","PKM","ALDOA","TPI1"),
  "ECM Remodeling"   = c("FN1","PDGFB","COL1A2","MMP14","TIMP1","VIM","COL4A1","COL6A1","SPARC"),
  "Oxidative Stress" = c("SOD2","PRDX1","CAT"),
  "Inflammation"     = c("CCL8","IL18","NLRP3","TNF","IL1B","NFKB1","RELB","CXCL8","JUNB"),
  "Antigen Presentation" = c("HLA-DRA","HLA-DRB1","CTSS","IFI30","CD74","HLA-DMA","HLA-DMB","FCER1G"),
  "Complement"       = c("C1QA","C1QB","C1QC","C3","C4A","CFB","CFH","SERPING1")
)

# Simplified internal names for AddModuleScore
short_names <- c("Prog","Prolif","Gly","ECM","OxS","Inflam","AgPres","Compl")

for(i in seq_along(mods)) {
  mods[[i]] <- intersect(mods[[i]], rownames(seu_normal))
  cat(sprintf("%-22s %d genes\n", names(mods)[i], length(mods[[i]])))
  seu_normal <- AddModuleScore(seu_normal, features=mods[i], name=short_names[i], ctrl=min(30, length(mods[[i]])))
}

meta <- seu_normal@meta.data
for(i in seq_along(short_names)) {
  old <- grep(paste0("^",short_names[i],"1$"), colnames(meta), value=TRUE)
  if(length(old)==1) colnames(meta)[colnames(meta)==old] <- short_names[i]
}

tri_cols <- c("Early"="#4575B4","Late"="#D73027")
module_names <- names(mods)

plots <- list()
for(i in seq_along(module_names)) {
  nm <- short_names[i]; disp <- module_names[i]
  wt <- wilcox.test(meta[meta$tri=="Late",nm], meta[meta$tri=="Early",nm])
  lbl <- ifelse(wt$p.value<0.001,"***",ifelse(wt$p.value<0.01,"**",ifelse(wt$p.value<0.05,"*","ns")))
  
  p <- ggplot(meta, aes_string(x="tri", y=nm, fill="tri")) +
    geom_violin(alpha=0.75, color=NA, scale="width") +
    geom_boxplot(width=0.15, outlier.size=0.15, alpha=0.5, fill="white", linewidth=0.25) +
    scale_fill_manual(values=tri_cols, guide="none") +
    labs(title=disp, y="", x="") +
    annotate("text", x=1.5, y=Inf, label=lbl, vjust=1.5, size=4.5, fontface="bold") +
    theme_bw() + theme(panel.grid=element_blank(), panel.border=element_rect(color="black",linewidth=0.35),
      axis.text=element_text(color="black",size=7), axis.text.x=element_text(size=8),
      plot.title=element_text(hjust=0.5,face="bold",size=8.5,lineheight=0.9),
      plot.margin=margin(2,3,2,3))
  plots[[i]] <- p
}

final <- wrap_plots(plots, ncol=4) +
  plot_annotation(title="Functional Module Scores: Early vs Late (Extended)",
    theme=theme(plot.title=element_text(hjust=0.5,face="bold",size=14)))

FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'
ggsave(file.path(FIGDIR,"FigS_module_scores.png"), final, w=12, h=7, dpi=300, bg="white")
cat("\nSaved FigS_module_scores.png\n")
