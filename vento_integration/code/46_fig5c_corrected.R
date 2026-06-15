#!/usr/bin/env Rscript
# Fig5c corrected: E-MTAB-12795 normal → Normal_Early (not Infection)
# Save to new files: deg_*_corrected.csv, disease_GSEA_corrected.csv, Fig5c_corrected.png
library(Seurat); library(clusterProfiler); library(org.Hs.eg.db)
library(ggplot2); library(dplyr); library(stringr)
set.seed(42)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
detail <- labels$disease_detail[1:ncol(seu)]

# Fix: E-MTAB-12795 normal → Normal_Early
detail_fixed <- detail
hoo_normal <- seu$dataset == "E-MTAB-12795" & detail == "normal"
cat(sprintf("Fixing %d E-MTAB-12795 normal cells → Normal_Early\n", sum(hoo_normal)))
detail_fixed[hoo_normal] <- "Normal_Early"

tri_map <- c("E-MTAB-12421"="Early","E-MTAB-6701"="Early","E-MTAB-12795"="Early",
             "GSE214607"="Early","UCSF Li 2026"="Mid","GSE290578"="Late",
             "GSE333257"="Late","GSE298602"="Late","GSE298119"="Late","GSE173193"="Late")
tri <- tri_map[as.character(seu$dataset)]

# Build disease_clean
dc <- ifelse(detail_fixed %in% c("Normal","Control","E-MTAB-12421","E-MTAB-6701","UCSF Li 2026","Normal_Early"), "Normal",
      ifelse(detail_fixed %in% c("PE","PreE_SF","gHTN","GSE173193","GSE298119"), "PE",
      ifelse(detail_fixed %in% c("RM","GSE214607"), "Miscarriage",
      ifelse(detail_fixed %in% c("toxoplasmosis","listeriosis","Plasmodium malariae malaria"), "Infection",
      ifelse(detail_fixed %in% c("PTL","PTNL"), "Preterm",
      ifelse(detail_fixed=="TL", "Excluded", "Other"))))))

dc[dc=="Normal" & tri=="Early"] <- "Normal_Early"
dc[dc=="Normal" & tri=="Late"] <- "Normal_Late"

# Save corrected labels
labels_corrected <- data.frame(dataset=seu$dataset, disease_corrected=dc, original_detail=detail_fixed)
write.csv(labels_corrected, file.path(OUTDIR, "per_cell_disease_corrected.csv"), row.names=FALSE)

cat("Corrected groups:\n")
for(g in c("Normal_Early","Normal_Late","Miscarriage","Infection","PE","Preterm","Excluded"))
  cat(sprintf("  %s: %d\n", g, sum(dc==g)))

# Run DEG + GSEA
keep <- dc %in% c("Normal_Early","Normal_Late","Miscarriage","PE","Infection","Preterm")
seu <- subset(seu, cells=colnames(seu)[keep]); dc <- dc[keep]
names(dc) <- colnames(seu)

run_gsea <- function(disease, control, label) {
  cells <- names(dc)[dc %in% c(disease, control)]
  sub <- subset(seu, cells=cells)
  sub$group <- ifelse(dc[cells]==disease, "Disease", "Control")
  Idents(sub) <- "group"; sub <- JoinLayers(sub)
  deg <- FindMarkers(sub, ident.1="Disease", ident.2="Control", logfc.threshold=0, min.pct=0.1, test.use="wilcox")
  write.csv(deg, file.path(OUTDIR, paste0("deg_",label,"_corrected.csv")))
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

cat("\nRunning corrected DEG + GSEA...\n")
d <- rbind(
  run_gsea("Miscarriage","Normal_Early","Miscarriage"),
  run_gsea("Infection","Normal_Early","Infection"),
  run_gsea("PE","Normal_Late","PE"),
  run_gsea("Preterm","Normal_Late","Preterm"))
d$comparison <- factor(d$comparison, levels=c("Miscarriage","Infection","PE","Preterm"))
write.csv(d, file.path(OUTDIR,"disease_GSEA_corrected.csv"), row.names=FALSE)

# Top cross-disease pathways
pc <- d %>% filter(p.adjust<0.1) %>% group_by(Description) %>% summarise(nd=n()) %>% filter(nd>=2)
top <- d %>% filter(Description %in% pc$Description) %>% group_by(comparison) %>% slice_max(abs(NES), n=12) %>% pull(Description) %>% unique()
pd <- d %>% filter(Description %in% top) %>% mutate(Description=ifelse(nchar(Description)>55, paste0(substr(Description,1,52),"..."), Description))
pd$comparison <- factor(pd$comparison, levels=c("Miscarriage","Infection","PE","Preterm"))

cat(sprintf("\nCross-disease pathways: %d plotted\n", length(unique(pd$Description))))

# Order
dl <- unique(pd$Description)
func <- setNames(rep(5,length(dl)),dl)
func[grep("translat|ribosom|spliceosom|cytoplas",dl)] <- 1
func[grep("mitochond|oxidative|atp|aerobic|respir|proton|electron|NADH|nucleoside",dl)] <- 2
func[grep("actin|cytoskel",dl)] <- 3
func[grep("antigen|mhc|immune|t cell|leukocyt|cytotox|humoral|antimicrob|interferon|viral|symbion|lipopolysacch|bacterial|interleukin|lymphocyte",dl)] <- 4
pg <- pd %>% group_by(Description) %>% summarise(g=first(func[Description]),.groups="drop") %>% arrange(g,Description)
pd$Description <- factor(pd$Description, levels=rev(pg$Description))

p <- ggplot(pd, aes(x=comparison, y=Description)) +
  geom_hline(yintercept=seq_along(levels(pd$Description)), color="grey92", linewidth=0.3) +
  geom_point(aes(size=Count, color=NES), stroke=0) +
  scale_color_gradient2(low="#313695",mid="#FFFFBF",high="#A50026",midpoint=0,name="NES",limits=c(-3.2,3.2)) +
  scale_size_continuous(range=c(2,7),name="Core\nGenes") + theme_bw() +
  labs(title="GO Biological Process — GSEA by Disease (corrected)",x="",y="",subtitle="Disease vs trimester-matched Normal") +
  theme(text=element_text(face="bold"),axis.text.x=element_text(size=12,color="black"),
    axis.text.y=element_text(size=9,color="black"),legend.title=element_text(size=10),
    legend.text=element_text(size=9),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    panel.border=element_rect(color="black",linewidth=0.6),
    plot.title=element_text(hjust=0.5,size=13),plot.subtitle=element_text(hjust=0.5,face="plain",size=10))
ggsave(file.path(FIGDIR,"Fig5c_GSEA_dotplot_corrected.png"), p, w=9, h=max(8, length(dl)*0.25), dpi=300, bg="white")
cat("Saved Fig5c_GSEA_dotplot_corrected.png\n")
