#!/usr/bin/env Rscript
# Fig1F: Module scores (Early vs Late) — adapted for new 10-dataset Atlas
library(Seurat); library(ggplot2); library(dplyr); library(patchwork)

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
cat(sprintf("Normal cells (Early+Late): %d\n", ncol(seu)))

seu <- JoinLayers(seu)

# Module definitions
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

short_names <- c("Prog","Prolif","Gly","ECM","OxS","Inflam","AgPres","Compl")

for(i in seq_along(mods)) {
  mods[[i]] <- intersect(mods[[i]], rownames(seu))
  cat(sprintf("%-22s %d genes\n", names(mods)[i], length(mods[[i]])))
  seu <- AddModuleScore(seu, features=mods[i], name=short_names[i], ctrl=min(30, length(mods[[i]])))
}

# Rename module score columns
meta <- seu@meta.data
for(i in seq_along(short_names)) {
  old <- grep(paste0("^",short_names[i],"1$"), colnames(meta), value=TRUE)
  if(length(old)==1) colnames(meta)[colnames(meta)==old] <- short_names[i]
}

t.test(meta$Prog[meta$tri=="Early"], meta$Prog[meta$tri=="Late"])$p.value
tri_cols <- c("Early"="#4575B4","Late"="#D73027")

plots <- list()
for(i in seq_along(names(mods))) {
  mod_name <- names(mods)[i]
  sn <- short_names[i]
  df <- meta[, c(sn, "tri")]
  colnames(df)[1] <- "Score"
  
  p <- ggplot(df, aes(tri, Score, fill=tri)) +
    geom_violin(scale="width", alpha=0.7, draw_quantiles=0.5, linewidth=0.3) +
    scale_fill_manual(values=tri_cols) +
    labs(title=mod_name, x="", y="Module Score") +
    theme_bw(9) + theme(legend.position="none", panel.grid=element_blank(),
      plot.title=element_text(size=10,face="bold"))
  plots[[i]] <- p
}

p_all <- wrap_plots(plots, ncol=4)
ggsave(file.path(FIGDIR,"Fig1F_module_scores.png"), p_all, w=16, h=8, dpi=300, bg="white")
cat("Saved Fig1F_module_scores.png\n")
