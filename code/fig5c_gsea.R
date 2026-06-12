#!/usr/bin/env Rscript
# Fig5c: GSEA Hallmark dotplot — 4 diseases vs trimester-matched normal
library(Seurat); library(ggplot2); library(dplyr); library(fgsea)
set.seed(42)

seu <- readRDS("results/Hofbauer_Atlas_Final.rds")
mask_gse <- seu$dataset == "GSE290578"; bcs <- colnames(seu)
is_norm_gse <- grepl("_Norm_", bcs) & mask_gse; is_pt_gse <- grepl("_Pt_", bcs) & mask_gse
seu$disease_clean <- NA
seu$disease_clean[seu$disease == "Normal_1st"] <- "Normal_Early"
seu$disease_clean[seu$disease == "Normal"] <- "Normal_Late"
seu$disease_clean[seu$disease == "RM/NC"] <- "Miscarriage"
seu$disease_clean[seu$disease == "Normal/Listeria/Toxoplasma/Malaria"] <- "Infection"
seu$disease_clean[is_norm_gse] <- "Normal_Late"; seu$disease_clean[is_pt_gse] <- "PE"
seu$disease_clean[seu$disease %in% c("PTL","PTNL")] <- "Preterm"
seu$disease_clean[seu$disease == "TL"] <- "Normal_Late"
keep <- seu$disease_clean %in% c("Normal_Early","Normal_Late","Miscarriage","PE","Infection","Preterm")
seu_sub <- subset(seu, cells=colnames(seu)[keep])
seu_sub$disease_clean <- as.character(seu_sub$disease_clean)

hallmarks <- gmtPathways("ref/h.all.v2023.2.Hs.symbols.gmt")

run_gsea <- function(disease, control, label) {
  cells <- colnames(seu_sub)[seu_sub$disease_clean %in% c(disease, control)]
  sub <- subset(seu_sub, cells=cells)
  sub$group <- ifelse(sub$disease_clean == disease, "Disease", "Control")
  Idents(sub) <- "group"
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
  geom_point() + scale_color_gradient2(low="#4575B4", mid="white", high="#D73027", midpoint=0) +
  scale_size_continuous(range=c(1.5, 6), name="-log10(FDR)") +
  labs(title="GSEA Hallmark: Disease vs trimester-matched Normal", x="", y="") +
  theme_bw(base_size=10) +
  theme(axis.text.y=element_text(size=8), axis.text.x=element_text(size=10,face="bold"),
        plot.title=element_text(face="bold",size=12,hjust=0.5),
        panel.grid.major=element_line(color="grey92"))
ggsave("figures/Fig5/Fig5c_GSEA_dotplot.png", pc, w=10, h=8, dpi=300, bg="white")

# Save DEGs
write.csv(res_mis$deg, "results/deg_Miscarriage_vs_NormalEarly.csv")
write.csv(res_inf$deg, "results/deg_Infection_vs_NormalEarly.csv")
write.csv(res_pe$deg,  "results/deg_PE_vs_NormalLate.csv")
write.csv(res_pt$deg,  "results/deg_Preterm_vs_NormalLate.csv")
cat("Fig5c + DEGs done\n")
