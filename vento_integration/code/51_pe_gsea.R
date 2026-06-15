#!/usr/bin/env Rscript
# GSEA: Early PE vs Late PE — each vs Normal_Late
library(Seurat); library(clusterProfiler); library(org.Hs.eg.db)
library(ggplot2); library(dplyr); library(stringr)
set.seed(42)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
detail <- labels$disease_detail[1:ncol(seu)]

# Groups: Early PE, Late PE, Normal_Late
seu$group <- NA_character_
seu$group[seu$dataset=="GSE290578" & detail=="PE"] <- "Early_PE"
seu$group[detail %in% c("PreE_SF","gHTN","GSE173193","GSE298119")] <- "Late_PE"
seu$group[detail=="Normal"] <- "Normal_Late"

seu <- subset(seu, group %in% c("Early_PE","Late_PE","Normal_Late"))

run_gsea <- function(disease, control, label) {
  sub <- subset(seu, group %in% c(disease, control))
  sub$grp <- ifelse(sub$group==disease, "Disease", "Control")
  Idents(sub) <- "grp"; sub <- JoinLayers(sub)
  deg <- FindMarkers(sub, ident.1="Disease", ident.2="Control", logfc.threshold=0, min.pct=0.1, test.use="wilcox")
  write.csv(deg, file.path(OUTDIR, paste0("deg_",label,"_vs_NormalLate.csv")))
  deg$symbol <- rownames(deg)
  conv <- bitr(deg$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  dm <- inner_join(deg, conv, by=c("symbol"="SYMBOL"))
  gl <- setNames(dm$avg_log2FC, dm$ENTREZID); gl <- gl[!is.na(names(gl))]; gl <- gl[!duplicated(names(gl))]; gl <- sort(gl, decreasing=TRUE)
  gse <- gseGO(geneList=gl, ont="BP", OrgDb=org.Hs.eg.db, pvalueCutoff=0.5, verbose=FALSE)
  if(nrow(as.data.frame(gse))>0) {
    gse <- setReadable(gse, OrgDb="org.Hs.eg.db", keyType="ENTREZID")
    res <- as.data.frame(gse); res$comparison <- label
    res$Count <- sapply(strsplit(as.character(res$core_enrichment), "/"), length)
    cat(sprintf("%s: %d DEG, %d terms (padj<0.1=%d)\n", label, nrow(deg), nrow(res), sum(res$p.adjust<0.1)))
    return(res)
  }
  NULL
}

cat("Running GSEA...\n")
d <- rbind(
  run_gsea("Early_PE","Normal_Late","Early_PE"),
  run_gsea("Late_PE","Normal_Late","Late_PE"))
d$comparison <- factor(d$comparison, levels=c("Early_PE","Late_PE"))
write.csv(d, file.path(OUTDIR,"PE_subgroup_GSEA.csv"), row.names=FALSE)

# Top pathways
top <- d %>% filter(p.adjust<0.1) %>% group_by(comparison) %>% slice_max(abs(NES), n=20) %>% pull(Description) %>% unique()
pd <- d %>% filter(Description %in% top) %>% mutate(Description=ifelse(nchar(Description)>55, paste0(substr(Description,1,52),"..."), Description))
pd$comparison <- factor(pd$comparison, levels=c("Early_PE","Late_PE"))
cat(sprintf("Plotting %d pathways\n", length(unique(pd$Description))))

# Order by functional group
dl <- unique(pd$Description)
func <- setNames(rep(5,length(dl)),dl)
func[grep("translat|ribosom|cytoplas",dl)] <- 1
func[grep("mitochond|oxidative|atp|aerobic|respir|proton|electron|NADH",dl)] <- 2
func[grep("actin|cytoskel|fiber|filament",dl)] <- 3
func[grep("antigen|mhc|immune|t cell|leukocyt|lymphocyte|interferon|inflammat|cytokin|chemokin",dl)] <- 4
pg <- pd %>% group_by(Description) %>% summarise(g=first(func[Description]),.groups="drop") %>% arrange(g,Description)
pd$Description <- factor(pd$Description, levels=rev(pg$Description))

p <- ggplot(pd, aes(x=comparison, y=Description)) +
  geom_hline(yintercept=seq_along(levels(pd$Description)), color="grey92", linewidth=0.3) +
  geom_point(aes(size=Count, color=NES), stroke=0) +
  scale_color_gradient2(low="#313695",mid="#FFFFBF",high="#A50026",midpoint=0,name="NES",limits=c(-3.2,3.2)) +
  scale_size_continuous(range=c(2,8),name="Core\nGenes") + theme_bw() +
  labs(title="GO:BP — Early PE vs Late PE (each vs Normal_Late)",x="",y="") +
  theme(text=element_text(face="bold"),axis.text.x=element_text(size=13,color="black"),
    axis.text.y=element_text(size=9,color="black"),legend.title=element_text(size=10),
    legend.text=element_text(size=9),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    panel.border=element_rect(color="black",linewidth=0.6),
    plot.title=element_text(hjust=0.5,size=12))

h <- max(6, length(dl)*0.28)
ggsave(file.path(FIGDIR,"PE_early_vs_late_GSEA.png"), p, w=9, h=h, dpi=300, bg="white")
cat("Saved\n")
