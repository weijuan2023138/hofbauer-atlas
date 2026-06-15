#!/usr/bin/env Rscript
# Fig1E: GSEA dotplot for ALL THREE trimesters (Early, Mid, Late)
library(Seurat); library(clusterProfiler); library(org.Hs.eg.db)
library(ggplot2); library(dplyr); library(stringr)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
seu$disease_final <- factor(labels$disease_final[1:ncol(seu)],
  levels=c("Normal","PE","Miscarriage","Infection","Preterm"))
seu <- subset(seu, disease_final=="Normal")

tri_map <- c("E-MTAB-12421"="Early","E-MTAB-6701"="Early","E-MTAB-12795"="Early",
             "GSE214607"="Early","UCSF Li 2026"="Mid","GSE290578"="Late",
             "GSE333257"="Late","GSE298602"="Late","GSE298119"="Late","GSE173193"="Late")
seu$tri <- factor(setNames(tri_map[as.character(seu$dataset)], colnames(seu)),
                  levels=c("Early","Mid","Late"))
print(table(seu$tri))

Idents(seu) <- "tri"
seu <- JoinLayers(seu)
all_gsea <- list()

for(tri in c("Early","Mid","Late")) {
  cat(sprintf("\n=== GSEA: %s vs others ===\n", tri))
  deg <- FindMarkers(seu, ident.1=tri, ident.2=NULL, logfc.threshold=0, min.pct=0.1, test.use="t")
  deg$symbol <- rownames(deg)
  gene_conv <- bitr(deg$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  deg_merged <- inner_join(deg, gene_conv, by=c("symbol"="SYMBOL"))
  gene_list <- setNames(deg_merged$avg_log2FC, deg_merged$ENTREZID)
  gene_list <- gene_list[!is.na(names(gene_list))]
  gene_list <- gene_list[!duplicated(names(gene_list))]
  gene_list <- sort(gene_list, decreasing=TRUE)
  cat(sprintf("  Genes: %d\n", length(gene_list)))
  gse <- gseGO(geneList=gene_list, ont="BP", OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, verbose=FALSE)
  if(nrow(as.data.frame(gse)) > 0) {
    gse <- setReadable(gse, OrgDb="org.Hs.eg.db", keyType="ENTREZID")
    all_gsea[[tri]] <- gse
    cat(sprintf("  Gene sets: %d\n", nrow(as.data.frame(gse))))
  }
}

combined <- data.frame()
for(tri in c("Early","Mid","Late")) {
  if(!is.null(all_gsea[[tri]])) {
    res <- as.data.frame(all_gsea[[tri]])
    if(nrow(res) > 0) {
      res$Trimester <- tri
      res$Count <- sapply(strsplit(as.character(res$core_enrichment), "/"), length)
      combined <- rbind(combined, res)
    }
  }
}
cat(sprintf("\nTotal: %d\n", nrow(combined)))
write.csv(combined, file.path(OUTDIR, "dev_trimester_GSEA_3tri.csv"), row.names=FALSE)

# Top 15 per trimester
top_pathways <- combined %>% group_by(Trimester) %>%
  slice_max(order_by=abs(NES), n=15) %>% pull(Description) %>% unique()

plot_data <- combined %>% filter(Description %in% top_pathways) %>%
  mutate(Description=str_wrap(Description,width=55),
         Trimester=factor(Trimester,levels=c("Early","Mid","Late")))

p <- ggplot(plot_data, aes(x=Trimester, y=reorder(Description, NES))) +
  geom_hline(yintercept=seq_along(unique(plot_data$Description)), color="grey92") +
  geom_point(aes(size=Count, color=NES)) +
  scale_color_gradient2(low="#377EB8",mid="white",high="#E41A1C",midpoint=0,name="NES") +
  scale_size_continuous(range=c(3,9),name="Core Genes") +
  theme_bw() + labs(title="GO Biological Process — GSEA by Trimester",x="",y="") +
  theme(text=element_text(face="bold"),axis.text.x=element_text(size=12,color="black"),
    axis.text.y=element_text(size=8.5,color="black"),legend.title=element_text(size=10),
    legend.text=element_text(size=9),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),panel.border=element_rect(color="black",linewidth=0.6),
    plot.title=element_text(hjust=0.5,size=13))
ggsave(file.path(FIGDIR,"Fig1E_dev_GSEA_dotplot.png"), p, w=12, h=16, dpi=300, bg="white")
cat("\nSaved Fig1E_dev_GSEA_dotplot.png\n")
