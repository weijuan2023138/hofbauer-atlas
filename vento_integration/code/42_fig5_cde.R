#!/usr/bin/env Rscript
# Fig5c/d/e: GSEA + gene dotplot + TF violins — new Atlas
library(Seurat); library(ggplot2); library(dplyr); library(fgsea); library(tidyr)
set.seed(42)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"
RESDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"
GMT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/ref/h.all.v2023.2.Hs.symbols.gmt"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
seu$disease_final <- factor(labels$disease_final[1:ncol(seu)], levels=c("Normal","PE","Miscarriage","Infection","Preterm"))

tri_map <- c("E-MTAB-12421"="Early","E-MTAB-6701"="Early","E-MTAB-12795"="Early",
             "GSE214607"="Early","UCSF Li 2026"="Mid","GSE290578"="Late",
             "GSE333257"="Late","GSE298602"="Late","GSE298119"="Late","GSE173193"="Late")
seu$trimester <- setNames(tri_map[as.character(seu$dataset)], colnames(seu))

seu$disease_clean <- as.character(seu$disease_final)
seu$disease_clean[seu$disease_final=="Normal" & seu$trimester=="Early"] <- "Normal_Early"
seu$disease_clean[seu$disease_final=="Normal" & seu$trimester=="Late"] <- "Normal_Late"
# Exclude TL from Normal_Late (match old pipeline)
tl_mask <- labels$disease_detail == "TL"
seu$disease_clean[tl_mask] <- "Excluded"
cat(sprintf("Excluded %d TL cells from Normal_Late\n", sum(tl_mask)))
keep <- seu$disease_clean %in% c("Normal_Early","Normal_Late","Miscarriage","PE","Infection","Preterm")
seu <- subset(seu, cells=colnames(seu)[keep])

# ═══ Fig5c: GSEA Hallmark ═══
cat("=== Fig5c: GSEA ===\n")
hallmarks <- gmtPathways(GMT)
run_gsea <- function(disease, control, label) {
  cells <- colnames(seu)[seu$disease_clean %in% c(disease, control)]
  sub <- subset(seu, cells=cells)
  sub$group <- ifelse(sub$disease_clean==disease, "Disease", "Control")
  Idents(sub) <- "group"; sub <- JoinLayers(sub)
  deg <- FindMarkers(sub, ident.1="Disease", ident.2="Control", logfc.threshold=0, min.pct=0.1, test.use="wilcox")
  deg$gene <- rownames(deg)
  ranks <- deg$avg_log2FC; names(ranks) <- deg$gene; ranks <- sort(ranks, decreasing=TRUE)
  gsea <- fgsea(pathways=hallmarks, stats=ranks, minSize=10, maxSize=500)
  gsea$comparison <- label
  list(deg=deg, gsea=gsea)
}

res_mis <- run_gsea("Miscarriage","Normal_Early","Miscarriage")
res_inf <- run_gsea("Infection","Normal_Early","Infection")
res_pe  <- run_gsea("PE","Normal_Late","PE")
res_pt  <- run_gsea("Preterm","Normal_Late","Preterm")

gsea_all <- rbind(res_mis$gsea, res_inf$gsea, res_pe$gsea, res_pt$gsea)
gsea_all$comparison <- factor(gsea_all$comparison, levels=c("Miscarriage","Infection","PE","Preterm"))

top_paths <- gsea_all %>% filter(padj < 0.1) %>%
  group_by(comparison) %>% slice_max(abs(NES), n=12) %>% pull(pathway) %>% unique()
gsea_plot <- gsea_all %>% filter(pathway %in% top_paths)
gsea_plot$pathway <- gsub("HALLMARK_","",gsea_plot$pathway)

pc <- ggplot(gsea_plot, aes(x=comparison, y=pathway, size=-log10(padj), color=NES)) +
  geom_point() + scale_color_gradient2(low="#4575B4",mid="white",high="#D73027",midpoint=0) +
  scale_size_continuous(range=c(1.5,6),name="-log10(FDR)") +
  labs(title="GSEA Hallmark: Disease vs trimester-matched Normal",x="",y="") +
  theme_bw(10) + theme(axis.text.y=element_text(size=8),axis.text.x=element_text(size=10,face="bold"),
    plot.title=element_text(face="bold",size=12,hjust=0.5),panel.grid.major=element_line(color="grey92"))
ggsave(file.path(FIGDIR,"Fig5c_GSEA_dotplot.png"), pc, w=10, h=8, dpi=300, bg="white")
write.csv(res_mis$deg, file.path(RESDIR,"deg_Miscarriage_vs_NormalEarly.csv"))
write.csv(res_inf$deg, file.path(RESDIR,"deg_Infection_vs_NormalEarly.csv"))
write.csv(res_pe$deg,  file.path(RESDIR,"deg_PE_vs_NormalLate.csv"))
write.csv(res_pt$deg,  file.path(RESDIR,"deg_Preterm_vs_NormalLate.csv"))
cat("Fig5c + DEGs saved\n")

# ═══ Fig5d: Gene Z-score dotplot ═══
cat("=== Fig5d: Gene dotplot ===\n")
genes <- c("SQSTM1","BNIP3","PRKN","TREM2","RELB","NFKB1","STAT3","STAT1","ID2","MAFB",
           "CEBPA","IFI27","FCGR2B","HLA-DQA1","HLA-DQB1","NNMT","HLA-G",
           "CFI","AEBP1","HIGD2A","COX14","ATP5MF",
           "PRG2","PSG9","TNF","IL1B","CXCL8","HLA-DRA","FCGR3A","C1QC","C1QB","C1QA","ITGB1",
           "ITGAV","CD44","PTPRM","TGFB1","COL1A2","FN1","SPP1")
genes <- intersect(genes, rownames(seu))

groups <- c("Normal_Early","Miscarriage","Infection","Normal_Late","PE","Preterm")
group_labels <- c("Normal\n(Early)","Miscarriage","Infection","Normal\n(Late)","PE","Preterm")

expr <- FetchData(seu, vars=c(genes, "disease_clean"))
plot_data <- data.frame()
for(g in genes) {
  grp_means <- sapply(groups, function(gr) mean(expr[expr$disease_clean==gr, g]))
  z <- (grp_means - mean(grp_means)) / sd(grp_means)
  for(i in seq_along(groups)) {
    plot_data <- rbind(plot_data, data.frame(Gene=g,Group=groups[i],GroupLabel=group_labels[i],Zscore=z[i],AbsZ=abs(z[i])))
  }
}
plot_data$Gene <- factor(plot_data$Gene, levels=genes)
plot_data$GroupLabel <- factor(plot_data$GroupLabel, levels=group_labels)

ph <- ggplot(plot_data, aes(x=Gene, y=GroupLabel)) +
  geom_point(aes(fill=Zscore, size=AbsZ), shape=21, stroke=0.3, color="grey80") +
  scale_fill_gradient2(low="#313695",mid="#FFFFBF",high="#A50026",midpoint=0,name="Z-score") +
  scale_size_continuous(range=c(2,8),name="|Z-score|",breaks=c(0.5,1,1.5,2)) +
  labs(title="Key gene expression",x="",y="") + theme_minimal(11) +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=9,color="black"),
    axis.text.y=element_text(size=10,color="black"),
    axis.line=element_line(color="black",linewidth=0.4),
    panel.grid.major=element_line(color="grey92",linewidth=0.3),panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold",size=14,hjust=0.5),
    legend.position="right",legend.title=element_text(size=10,face="bold"),legend.text=element_text(size=9))
ggsave(file.path(FIGDIR,"Fig5d_gene_dotplot.png"), ph, w=12, h=4.5, dpi=300, bg="white")
cat("Fig5d done\n")

# ═══ Fig5e: TF violins ═══
cat("=== Fig5e: TF violins ===\n")
fill_cols <- c("#4575B4","#D73027","#FDB462","#66C2A5","#FC8D59","#E41A1C")
all_tfs <- c("STAT1","STAT3","NFKB1","RELB","CEBPA","MAFB","ID2","KLF4","JUN","FOS","IRF1","IRF8")

mis <- read.csv(file.path(RESDIR,"deg_Miscarriage_vs_NormalEarly.csv"), row.names=1)
inf <- read.csv(file.path(RESDIR,"deg_Infection_vs_NormalEarly.csv"), row.names=1)
pe  <- read.csv(file.path(RESDIR,"deg_PE_vs_NormalLate.csv"), row.names=1)
pt  <- read.csv(file.path(RESDIR,"deg_Preterm_vs_NormalLate.csv"), row.names=1)
get_padj <- function(deg, gene) { if(gene %in% rownames(deg) && !is.na(deg[gene,"p_val_adj"])) deg[gene,"p_val_adj"] else 1 }

seu <- JoinLayers(seu)
expr_mat <- GetAssayData(seu, assay="RNA", layer="data"); meta <- seu@meta.data

for(g in all_tfs) {
  if(!g %in% rownames(expr_mat)) next
  df <- data.frame(expr=expr_mat[g,], Group=factor(meta$disease_clean, levels=groups))
  max_y <- max(df$expr) * 1.15; ys <- max_y * c(1.06, 1.12)
  
  p <- ggplot(df, aes(x=Group, y=expr, fill=Group)) +
    geom_violin(scale="width", trim=TRUE, linewidth=0.3, alpha=0.8) +
    geom_boxplot(width=0.04, outlier.size=0.02, alpha=0.5, fill="white", linewidth=0.2, fatten=1.5) +
    scale_fill_manual(values=fill_cols)
  
  comps <- list(list(1,2,1.5,ys[1],get_padj(mis,g)), list(1,3,2.5,ys[2],get_padj(inf,g)),
                list(4,5,4.5,ys[1],get_padj(pe,g)),  list(4,6,5,ys[2],get_padj(pt,g)))
  for(cp in comps) {
    pv <- cp[[5]]
    lbl <- ifelse(pv<0.001,"***",ifelse(pv<0.01,"**",ifelse(pv<0.05,"*","ns")))
    p <- p + annotate("segment",x=cp[[1]],xend=cp[[2]],y=cp[[4]],yend=cp[[4]],color="black",linewidth=0.5) +
      annotate("text",x=cp[[3]],y=cp[[4]]*1.005,label=lbl,size=2.8,fontface="bold",vjust=0)
  }
  
  p <- p + labs(title=g, x="", y="Expression") + theme_classic(10) +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=10,color="black",face="bold"),
      axis.text.y=element_text(size=8,color="black"),
      axis.line=element_line(color="black",linewidth=0.4),
      plot.title=element_text(face="bold",size=13,hjust=0.5,color="black"),
      legend.position="none")
  
  if(g %in% c("STAT1","NFKB1","CEBPA","JUN")) {
    ggsave(file.path(FIGDIR,sprintf("Fig5e_TF_%s.png",g)), p, w=5, h=4, dpi=300, bg="white")
  } else {
    ggsave(file.path(FIGDIR,sprintf("补充Fig5D_TF_%s.png",g)), p, w=5, h=4, dpi=300, bg="white")
  }
}
cat("Fig5e done\n")
