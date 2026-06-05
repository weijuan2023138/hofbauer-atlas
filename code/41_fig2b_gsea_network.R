#!/usr/bin/env Rscript
# Fig2B: GSEA enrichment curves + TF-pathway network
library(Seurat); library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2); library(dplyr); library(enrichplot)

seu <- readRDS('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/Hofbauer_Atlas_Final.rds')
cond <- read.csv("/tmp/gse290578_conditions.csv", stringsAsFactors=FALSE)
rownames(cond) <- cond$barcode; seu$gse <- NA
for(bc in colnames(seu)) if(bc %in% cond$barcode) seu$gse[bc] <- cond[bc,"condition"]
seu$tri <- NA
seu$tri[seu$dataset=="Arutyunyan"] <- "Early"
seu$tri[seu$dataset=="GSE290578" & seu$gse=="Normal"] <- "Late"
seu_normal <- subset(seu, cells=colnames(seu)[!is.na(seu$tri)])

# ---- GSEA for Late vs Early ----
Idents(seu_normal) <- "tri"
deg <- FindMarkers(seu_normal, ident.1="Late", ident.2="Early", logfc.threshold=0, min.pct=0.1)
deg$symbol <- rownames(deg)
gene_conv <- bitr(deg$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
deg_m <- inner_join(deg, gene_conv, by=c("symbol"="SYMBOL"))
gene_list <- setNames(deg_m$avg_log2FC, deg_m$ENTREZID)
gene_list <- sort(gene_list[!duplicated(names(gene_list))], decreasing=TRUE)

gse <- gseGO(geneList=gene_list, ont="BP", OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, verbose=FALSE)

FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'

# Pick 4 key pathways
pathways <- c("GO:0006956","GO:0002376","GO:0051301","GO:0006260")  # complement, immune, cell division, DNA replication
names(pathways) <- c("Complement Activation","Immune System Process","Cell Division","DNA Replication")

for(i in seq_along(pathways)) {
  pid <- pathways[i]; pname <- names(pathways)[i]
  png(file.path(FIGDIR, paste0("Fig2B_GSEA_",i,".png")), w=5, h=4, units="in", res=300)
  print(gseaplot2(gse, geneSetID=pid, title=pname, color="#D73027"))
  dev.off()
  cat(pname, "\n")
}

# ---- TF-Pathway Network ----
# Get leading edge genes for each GO term, map to TFs
go_df <- as.data.frame(gse)
top_go <- go_df %>% filter(NES > 1.5) %>% head(10)

early_tfs <- c("CEBPA","ID2","MYC","HMGA1","SOX4")
late_tfs  <- c("NR4A3","NFKB1","RELB","EGR1","KLF4","BCL6")

network <- data.frame()
for(i in 1:nrow(top_go)) {
  genes <- strsplit(top_go$core_enrichment[i], "/")[[1]]
  gene_conv <- tryCatch(bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"), error=function(e) NULL)
  if(is.null(gene_conv)) next
  pathway_genes <- gene_conv$SYMBOL
  
  for(tf in c(early_tfs, late_tfs)) {
    if(tf %in% pathway_genes) {
      network <- rbind(network, data.frame(TF=tf, Pathway=top_go$Description[i],
        NES=top_go$NES[i], TF_type=ifelse(tf %in% early_tfs, "Early","Late")))
    }
  }
}
cat(sprintf("\nTF-Pathway edges: %d\n", nrow(network)))
write.csv(network, file.path(FIGDIR, "../results/TF_pathway_network.csv"), row.names=FALSE)

# Network visualization
if(nrow(network) > 3) {
  p_net <- ggplot(network, aes(x=TF, y=Pathway, color=TF_type, size=abs(NES))) +
    geom_point() +
    scale_color_manual(values=c("Early"="#4575B4","Late"="#D73027")) +
    scale_size_continuous(range=c(3,10)) +
    labs(title="TF-Pathway Regulatory Network", x="", y="") +
    theme_bw() + theme(panel.grid=element_blank(),
      axis.text.x=element_text(angle=45, hjust=1, size=9),
      plot.title=element_text(face="bold",hjust=0.5))
  ggsave(file.path(FIGDIR, "Fig2F_TF_network.png"), p_net, w=8, h=5, dpi=200)
  cat("Saved Fig2F_TF_network.png\n")
}
cat("Done\n")
