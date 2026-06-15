#!/usr/bin/env Rscript
# GSE298602 internal: PreE_SF (severe PE) vs Control — same-dataset analysis
library(Seurat); library(clusterProfiler); library(org.Hs.eg.db)
library(ggplot2); library(dplyr); library(stringr)
set.seed(42)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
detail <- labels$disease_detail[1:ncol(seu)]

# GSE298602 only: PreE_SF vs Control
mask <- seu$dataset=="GSE298602" & detail %in% c("PreE_SF","Control")
seu602 <- subset(seu, cells=colnames(seu)[mask])
seu602$grp <- ifelse(detail[mask]=="PreE_SF", "PreE_SF", "Control")
Idents(seu602) <- "grp"; seu602 <- JoinLayers(seu602)

cat(sprintf("GSE298602: PreE_SF=%d, Control=%d\n", sum(seu602$grp=="PreE_SF"), sum(seu602$grp=="Control")))

# ── 1. Subtype proportions ──
prop <- prop.table(table(seu602$subtype_pred, seu602$grp), margin=2)*100
cat("=== Subtype proportions ===\n")
print(round(prop, 1))

# ── 2. Key gene expression ──
key_genes <- c("FLT1","FN1","PAPPA","HLA-DRA","CEBPA","SPP1","COL1A2","IL1B","TNF")
for(g in intersect(key_genes, rownames(seu602))) {
  expr <- FetchData(seu602, g)
  pe <- mean(expr[seu602$grp=="PreE_SF",1]); ct <- mean(expr[seu602$grp=="Control",1])
  cat(sprintf("  %-10s  PreE_SF=%.3f  Control=%.3f  FC=%.2f\n", g, pe, ct, pe/ct))
}

# ── 3. DEG + GSEA ──
deg <- FindMarkers(seu602, ident.1="PreE_SF", ident.2="Control", logfc.threshold=0, min.pct=0.05, test.use="wilcox")
write.csv(deg, file.path(OUTDIR,"deg_GSE298602_PreESF_vs_Control.csv"))
deg$symbol <- rownames(deg)
conv <- bitr(deg$symbol, from="SYMBOL", to="ENTREZID", Org=org.Hs.eg.db)
dm <- inner_join(deg, conv, by=c("symbol"="SYMBOL"))
gl <- setNames(dm$avg_log2FC, dm$ENTREZID); gl <- gl[!is.na(names(gl))]; gl <- gl[!duplicated(names(gl))]; gl <- sort(gl, decreasing=TRUE)
gse <- gseGO(geneList=gl, ont="BP", OrgDb=org.Hs.eg.db, pvalueCutoff=0.5, verbose=FALSE)
gse <- setReadable(gse, OrgDb="org.Hs.eg.db", keyType="ENTREZID")
res <- as.data.frame(gse)

sig <- res %>% filter(p.adjust<0.1)
cat(sprintf("\nDEG: %d genes | GSEA sig: %d\n", nrow(deg), nrow(sig)))

cat("=== Top PreE_SF UP (NES>0) ===\n")
for(i in 1:min(6, sum(sig$NES>0))) cat(sprintf("  %.2f  %s\n", sig$NES[sig$NES>0][i], substr(sig$Description[sig$NES>0][i],1,60)))
cat("=== Top PreE_SF DOWN (NES<0) ===\n")
for(i in 1:min(6, sum(sig$NES<0))) cat(sprintf("  %.2f  %s\n", sig$NES[sig$NES<0][i], substr(sig$Description[sig$NES<0][i],1,60)))

write.csv(res, file.path(OUTDIR,"GSE298602_PreESF_vs_Control_GSEA.csv"), row.names=FALSE)

# ── 4. Bar chart ──
top20 <- sig %>% slice_max(abs(NES), n=20) %>% mutate(Description=str_wrap(Description,width=50))
top20 <- top20 %>% arrange(-NES)
top20$Description <- factor(top20$Description, levels=rev(top20$Description))
top20$direction <- ifelse(top20$NES>0, "PreE_SF UP", "PreE_SF DOWN")

p <- ggplot(top20, aes(x=Description, y=NES, fill=direction)) +
  geom_col(width=0.6, alpha=0.9) +
  scale_fill_manual(values=c("PreE_SF UP"="#C62828","PreE_SF DOWN"="#1565C0")) +
  geom_hline(yintercept=0, color="black", linewidth=0.6) + coord_flip() +
  labs(title="GO:BP — GSE298602 PreE_SF vs Control", x="", y="NES") +
  theme_minimal(12) +
  theme(axis.text.y=element_text(size=9.5,color="black",face="bold"),
    axis.text.x=element_text(size=9,color="black"),
    panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
    panel.grid.major.x=element_line(color="grey92",linewidth=0.3),
    plot.title=element_text(hjust=0.5,size=13,face="bold"),
    legend.position="top",legend.title=element_blank(),legend.text=element_text(size=11))
ggsave(file.path(FIGDIR,"GSE298602_PreESF_vs_Control_bar.png"), p, w=8, h=6, dpi=300, bg="white")
cat("Saved\n")
