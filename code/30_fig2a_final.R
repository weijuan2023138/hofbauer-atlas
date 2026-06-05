#!/usr/bin/env Rscript
# Fig2a: 5-module scores (matches dot plot + proliferation + complement)
library(Seurat); library(ggplot2); library(dplyr); library(patchwork)

seu <- readRDS('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/Hofbauer_Atlas_Final.rds')
cond <- read.csv("/tmp/gse290578_conditions.csv", stringsAsFactors=FALSE)
rownames(cond) <- cond$barcode; seu$gse <- NA
for(bc in colnames(seu)) if(bc %in% cond$barcode) seu$gse[bc] <- cond[bc,"condition"]
seu$tri <- NA; seu$tri[seu$dataset=="Arutyunyan"] <- "Early"
seu$tri[seu$dataset=="GSE290578" & seu$gse=="Normal"] <- "Late"
seu_normal <- subset(seu, cells=colnames(seu)[!is.na(seu$tri)])

mod_names <- c("Progenitor","Proliferation","Remodeling","Immunity","Complement")
display <- c("Progenitor"="Progenitor","Proliferation"="Proliferation",
  "Remodeling"="Remodeling","Immunity"="Immunity","Complement"="Complement")

modules <- list(
  c("TREM2","AXL","CEBPA","ID2","CD5L","NOTCH2"),
  c("MKI67","TOP2A","CCNA2","PCNA","BIRC5","CDK1"),
  c("TIMP1","VIM","MMP14","ENO1","PGK1","COL1A2","FN1","PDGFB","SOD2"),
  c("CCL8","IL18","IFI30","HLA-DRA","CTSS","FCGR3A"),
  c("C1QA","C1QB","C1QC","C3")
)
names(modules) <- mod_names

for(i in seq_along(modules)) {
  modules[[i]] <- intersect(modules[[i]], rownames(seu_normal))
  seu_normal <- AddModuleScore(seu_normal, features=modules[i], name=mod_names[i], ctrl=30)
}
meta <- seu_normal@meta.data
for(nm in mod_names) {
  old <- grep(paste0("^",nm,"1$"), colnames(meta), value=TRUE)
  if(length(old)==1) colnames(meta)[colnames(meta)==old] <- nm
}

tri_cols <- c("Early"="#4575B4","Late"="#D73027")
stats <- data.frame()
for(nm in mod_names) {
  wt <- wilcox.test(meta[meta$tri=="Late",nm], meta[meta$tri=="Early",nm])
  stats <- rbind(stats, data.frame(Module=nm, delta=mean(meta[meta$tri=="Late",nm])-mean(meta[meta$tri=="Early",nm]), pval=wt$p.value))
}
stats$p_adj <- p.adjust(stats$pval, "BH")
stats$label <- ifelse(stats$p_adj<0.001,"***",ifelse(stats$p_adj<0.01,"**",ifelse(stats$p_adj<0.05,"*","ns")))

plots <- list()
for(i in seq_along(mod_names)) {
  nm <- mod_names[i]; disp <- display[nm]; lbl <- stats$label[i]
  p <- ggplot(meta, aes_string(x="tri", y=nm, fill="tri")) +
    geom_violin(alpha=0.75, color=NA, scale="width") +
    geom_boxplot(width=0.15, outlier.size=0.15, alpha=0.5, fill="white", linewidth=0.25) +
    scale_fill_manual(values=tri_cols, guide="none") +
    labs(title=disp, y="", x="") +
    annotate("text", x=1.5, y=Inf, label=lbl, vjust=1.5, size=4.5, fontface="bold") +
    theme_bw() + theme(panel.grid=element_blank(), panel.border=element_rect(color="black",linewidth=0.35),
      axis.text=element_text(color="black",size=8), axis.text.x=element_text(size=9),
      plot.title=element_text(hjust=0.5,face="bold",size=10), plot.margin=margin(3,4,3,4))
  plots[[nm]] <- p
}

final <- wrap_plots(plots, ncol=5) +
  plot_annotation(title="Functional Module Scores: Early vs Late",
    theme=theme(plot.title=element_text(hjust=0.5,face="bold",size=14)))

FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'
ggsave(file.path(FIGDIR,"Fig2a_module_violin.png"), final, w=14, h=3.2, dpi=300, bg="white")
cat("Saved Fig2a_module_violin.png (5 modules)\n")
for(i in 1:nrow(stats)) cat(sprintf("%-15s delta=%+.3f %s\n", stats$Module[i], stats$delta[i], stats$label[i]))
