#!/usr/bin/env Rscript
# GSEA-style GO dot plot across three trimesters
# Reference style: 母体和胎儿细胞差异分析.R lines 306-373
library(Seurat); library(clusterProfiler); library(org.Hs.eg.db)
library(ggplot2); library(dplyr); library(stringr)

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

# ---- Run GSEA per trimester ----
Idents(seu_normal) <- "tri"
all_gsea <- list()

for(tri in levels(seu_normal$tri)) {
  cat("\n=== GSEA for", tri, "===\n")
  
  # DEG: this trimester vs all others (full gene list, no filtering)
  deg <- FindMarkers(seu_normal, ident.1=tri, ident.2=NULL,
                     logfc.threshold=0, min.pct=0.1, test.use="t")
  deg$symbol <- rownames(deg)
  
  # Symbol → Entrez
  gene_conv <- bitr(deg$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  deg_merged <- inner_join(deg, gene_conv, by=c("symbol"="SYMBOL"))
  
  # Ranked gene list for GSEA
  gene_list <- setNames(deg_merged$avg_log2FC, deg_merged$ENTREZID)
  gene_list <- gene_list[!is.na(names(gene_list))]
  gene_list <- gene_list[!duplicated(names(gene_list))]
  gene_list <- sort(gene_list, decreasing=TRUE)
  
  cat(sprintf("  Genes in ranked list: %d\n", length(gene_list)))
  
  gse <- gseGO(geneList=gene_list, ont="BP", OrgDb=org.Hs.eg.db,
               pvalueCutoff=0.05, verbose=FALSE)
  
  if(nrow(as.data.frame(gse)) > 0) {
    gse <- setReadable(gse, OrgDb="org.Hs.eg.db", keyType="ENTREZID")
    all_gsea[[tri]] <- gse
    cat(sprintf("  Significant gene sets: %d\n", nrow(as.data.frame(gse))))
  } else {
    cat("  No significant enrichment\n")
  }
}

# ---- Combine into one dot plot (reference style) ----
combined <- data.frame()
trimesters <- c("Early","Mid","Late")

for(tri in trimesters) {
  if(!is.null(all_gsea[[tri]])) {
    res <- as.data.frame(all_gsea[[tri]])
    if(nrow(res) > 0) {
      res$Trimester <- tri
      res$Count <- sapply(strsplit(as.character(res$core_enrichment), "/"), length)
      combined <- rbind(combined, res)
    }
  }
}

cat(sprintf("\nTotal GSEA results across all trimesters: %d\n", nrow(combined)))

# Pick top pathways: for each trimester, take top NES terms
# Keep both activated (NES>0) and suppressed (NES<0)
top_pathways <- combined %>%
  group_by(Trimester) %>%
  slice_max(order_by=abs(NES), n=15) %>%
  pull(Description) %>%
  unique()

plot_data <- combined %>%
  filter(Description %in% top_pathways) %>%
  mutate(
    Trimester = factor(Trimester, levels=trimesters),
    Description = str_wrap(Description, width=55)
  )

# Remove pathways that appear in only 1 trimester (keep cross-trimester comparison)
pathway_freq <- table(plot_data$Description)
shared_pathways <- names(pathway_freq[pathway_freq >= 2])

if(length(shared_pathways) >= 8) {
  plot_data <- plot_data %>% filter(Description %in% shared_pathways)
  cat(sprintf("Filtered to %d shared pathways (appear in ≥2 trimesters)\n", length(shared_pathways)))
} else {
  cat("Keeping all pathways (not enough shared across trimesters)\n")
}

# ---- Plot ----
p <- ggplot(plot_data, aes(x=Trimester, y=reorder(Description, NES))) +
  geom_hline(yintercept=seq_along(unique(plot_data$Description)), color="grey92") +
  geom_point(aes(size=Count, color=NES)) +
  scale_color_gradient2(low="#377EB8", mid="white", high="#E41A1C",
                        midpoint=0, name="NES") +
  scale_size_continuous(range=c(3, 9), name="Core Genes") +
  theme_bw() +
  labs(title="GO Biological Process — GSEA by Trimester",
       subtitle="Early / Mid / Late vs. Other Trimesters",
       x="", y="") +
  theme(
    text=element_text(face="bold"),
    axis.text.x=element_text(size=12, color="black"),
    axis.text.y=element_text(size=8.5, color="black"),
    legend.title=element_text(size=10),
    legend.text=element_text(size=9),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.border=element_rect(color="black", linewidth=0.6),
    plot.title=element_text(hjust=0.5, size=13),
    plot.subtitle=element_text(hjust=0.5, face="plain", size=10)
  )

FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'
ggsave(file.path(FIGDIR,"Fig_dev_GSEA_dotplot.png"), p, w=10, h=14, dpi=300, bg="white")
cat("\nSaved Fig_dev_GSEA_dotplot.png\n")

# Save data
write.csv(combined, '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/dev_trimester_GSEA.csv', row.names=FALSE)
cat("Done\n")
