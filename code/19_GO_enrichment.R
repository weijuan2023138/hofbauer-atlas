#!/usr/bin/env Rscript
# GO enrichment by trimester — GSEA-style dot plot
# Pathways on left (y-axis), trimesters on bottom (x-axis)
library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2); library(dplyr)

deg_all <- read.csv('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/dev_trimester_DEGs.csv')
FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'

# Helper: gene symbol → ENTREZID
symbol2entrez <- function(genes) {
  map <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  return(map$ENTREZID)
}

tri_cols <- c("Early"="Early", "Mid"="Mid", "Late"="Late")
go_all <- data.frame()

for(tri in c("Early","Mid","Late")) {
  deg_tri <- deg_all[deg_all$trimester == tri,]
  # Top 500 DEGs by log2FC for enrichment
  deg_tri <- deg_tri %>% arrange(desc(avg_log2FC)) %>% head(500)
  entrez <- symbol2entrez(deg_tri$gene)
  cat(sprintf("%s: %d/%d genes mapped to ENTREZ\n", tri, length(entrez), nrow(deg_tri)))
  
  ego <- enrichGO(gene=entrez, OrgDb=org.Hs.eg.db, ont="BP",
    pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2, readable=TRUE)
  
  if(nrow(as.data.frame(ego)) > 0) {
    ego_df <- as.data.frame(ego)
    # Keep top 10 per trimester
    ego_df <- ego_df %>% arrange(p.adjust) %>% head(10)
    ego_df$Trimester <- tri
    go_all <- rbind(go_all, ego_df)
  }
}

# Calculate GeneRatio for dot size
go_all$GeneRatio <- sapply(go_all$GeneRatio, function(x) eval(parse(text=x)))
go_all$Trimester <- factor(go_all$Trimester, levels=c("Early","Mid","Late"))
go_all$negLogP <- -log10(go_all$p.adjust)

# Pick top pathways across all trimesters (top 20 overall)
top_pathways <- go_all %>% group_by(Description) %>% 
  summarise(min_p=min(p.adjust)) %>% arrange(min_p) %>% head(30) %>% pull(Description)
go_plot <- go_all[go_all$Description %in% top_pathways,]

# Dot plot: pathways on left, trimesters on bottom
p <- ggplot(go_plot, aes(x=Trimester, y=Description)) +
  geom_point(aes(size=GeneRatio, color=negLogP), stroke=0) +
  scale_color_gradientn(colors=c("#4DBBD5","#E18727","#C62828"), 
                         name="-log10(P.adj)") +
  scale_size_continuous(range=c(3,10), name="GeneRatio") +
  theme_bw() +
  theme(
    panel.grid.major=element_line(color="grey92", linewidth=0.4),
    panel.grid.minor=element_blank(),
    panel.border=element_rect(color="black", linewidth=0.6),
    axis.text.y=element_text(size=9, color="black"),
    axis.text.x=element_text(size=11, color="black", face="bold"),
    axis.title=element_blank(),
    legend.position="right",
    legend.box="vertical",
    legend.text=element_text(size=9),
    legend.title=element_text(size=10)
  )

ggsave(file.path(FIGDIR,"Fig_dev_GO_dotplot.png"), p, w=12, h=9, dpi=300, bg="white")
cat("Saved GO dot plot\n")

# Also save the GO table for reference
write.csv(go_all, '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/dev_trimester_GO.csv', row.names=FALSE)
cat("Done\n")
