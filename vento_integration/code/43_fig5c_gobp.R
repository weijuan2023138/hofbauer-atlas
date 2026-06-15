#!/usr/bin/env Rscript
# Fig5c: GO:BP GSEA dotplot — match old format exactly
library(Seurat); library(clusterProfiler); library(org.Hs.eg.db)
library(ggplot2); library(dplyr); library(stringr)
set.seed(42)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
seu$disease_final <- factor(labels$disease_final[1:ncol(seu)], levels=c("Normal","PE","Miscarriage","Infection","Preterm"))
tri_map <- c("E-MTAB-12421"="Early","E-MTAB-6701"="Early","E-MTAB-12795"="Early","GSE214607"="Early","UCSF Li 2026"="Mid","GSE290578"="Late","GSE333257"="Late","GSE298602"="Late","GSE298119"="Late","GSE173193"="Late")
seu$tri <- setNames(tri_map[as.character(seu$dataset)], colnames(seu))
seu$disease_clean <- as.character(seu$disease_final)
seu$disease_clean[seu$disease_final=="Normal" & seu$tri=="Early"] <- "Normal_Early"
seu$disease_clean[seu$disease_final=="Normal" & seu$tri=="Late"] <- "Normal_Late"
seu$disease_clean[labels$disease_detail=="TL"] <- "Excluded"
keep <- seu$disease_clean %in% c("Normal_Early","Normal_Late","Miscarriage","PE","Infection","Preterm")
seu <- subset(seu, cells=colnames(seu)[keep])
seu$disease_clean <- as.character(seu$disease_clean)

run_gsea <- function(disease, control, label) {
  cells <- colnames(seu)[seu$disease_clean %in% c(disease, control)]
  sub <- subset(seu, cells=cells)
  sub$group <- ifelse(sub$disease_clean==disease, "Disease", "Control")
  Idents(sub) <- "group"; sub <- JoinLayers(sub)
  deg <- FindMarkers(sub, ident.1="Disease", ident.2="Control", logfc.threshold=0, min.pct=0.1, test.use="wilcox")
  deg$symbol <- rownames(deg)
  gene_conv <- bitr(deg$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  deg_merged <- inner_join(deg, gene_conv, by=c("symbol"="SYMBOL"))
  gene_list <- setNames(deg_merged$avg_log2FC, deg_merged$ENTREZID)
  gene_list <- gene_list[!is.na(names(gene_list))]
  gene_list <- gene_list[!duplicated(names(gene_list))]
  gene_list <- sort(gene_list, decreasing=TRUE)
  gse <- gseGO(geneList=gene_list, ont="BP", OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, verbose=FALSE)
  if(nrow(as.data.frame(gse)) > 0) {
    gse <- setReadable(gse, OrgDb="org.Hs.eg.db", keyType="ENTREZID")
    res <- as.data.frame(gse)
    res$comparison <- label
    res$Count <- sapply(strsplit(as.character(res$core_enrichment), "/"), length)
    return(res)
  }
  return(NULL)
}

cat("Running GO:BP GSEA...\n")
res_mis <- run_gsea("Miscarriage","Normal_Early","Miscarriage")
res_inf <- run_gsea("Infection","Normal_Early","Infection")
res_pe  <- run_gsea("PE","Normal_Late","PE")
res_pt  <- run_gsea("Preterm","Normal_Late","Preterm")

combined <- rbind(res_mis, res_inf, res_pe, res_pt)
combined$comparison <- factor(combined$comparison, levels=c("Miscarriage","Infection","PE","Preterm"))
cat(sprintf("Total: %d results\n", nrow(combined)))

# Top pathways
top_paths <- combined %>% filter(p.adjust < 0.1) %>%
  group_by(comparison) %>% slice_max(abs(NES), n=15) %>% pull(Description) %>% unique()
plot_data <- combined %>% filter(Description %in% top_paths) %>%
  mutate(Description = str_wrap(Description, width=50))

cat(sprintf("Plotting %d pathways\n", length(unique(plot_data$Description))))

p <- ggplot(plot_data, aes(x=comparison, y=Description)) +
  geom_hline(yintercept=seq_along(unique(plot_data$Description)), color="grey92", linewidth=0.3) +
  geom_point(aes(size=Count, color=NES), stroke=0) +
  scale_color_gradient2(low="#4575B4", mid="white", high="#D73027", midpoint=0, name="NES") +
  scale_size_continuous(range=c(2, 7), name="Core\nGenes") +
  labs(title="GO Biological Process — GSEA by Disease", x="", y="",
       subtitle="Disease vs trimester-matched Normal") +
  theme_bw(base_size=10) +
  theme(axis.text.y=element_text(size=8, color="black"),
        axis.text.x=element_text(size=11, color="black", face="bold"),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        panel.grid.major=element_line(color="grey92", linewidth=0.2),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black", linewidth=0.5),
        plot.title=element_text(hjust=0.5, size=13, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=10, face="plain"))

ggsave(file.path(FIGDIR,"Fig5c_GSEA_dotplot.png"), p, w=11, h=10, dpi=300, bg="white")
cat("Saved\n")
