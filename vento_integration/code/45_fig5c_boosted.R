#!/usr/bin/env Rscript
# Fig5c boosted: PTNL→Preterm, t-test, relaxed GSEA cutoff
library(Seurat); library(clusterProfiler); library(org.Hs.eg.db)
library(ggplot2); library(dplyr); library(stringr)
set.seed(42)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
seu$disease_final <- factor(labels$disease_final[1:ncol(seu)], levels=c("Normal","PE","Miscarriage","Infection","Preterm"))

# Fix: PTNL → Preterm (not Normal)
tri_map <- c("E-MTAB-12421"="Early","E-MTAB-6701"="Early","E-MTAB-12795"="Early","GSE214607"="Early","UCSF Li 2026"="Mid","GSE290578"="Late","GSE333257"="Late","GSE298602"="Late","GSE298119"="Late","GSE173193"="Late")
seu$tri <- setNames(tri_map[as.character(seu$dataset)], colnames(seu))

# Build disease_clean with PTNL in Preterm
detail <- labels$disease_detail[1:ncol(seu)]
dc <- as.character(labels$disease_final)[1:ncol(seu)]
dc[dc=="Normal" & seu$tri=="Early"] <- "Normal_Early"
dc[dc=="Normal" & seu$tri=="Late"] <- "Normal_Late"
dc[detail=="TL"] <- "Excluded"           # TL excluded
dc[detail=="PTNL"] <- "Preterm"          # PTNL → Preterm (old behavior)
dc[dc=="Preterm"] <- "Preterm"
keep <- dc %in% c("Normal_Early","Normal_Late","Miscarriage","PE","Infection","Preterm")
seu <- subset(seu, cells=colnames(seu)[keep])
dc <- dc[keep]
cat("Groups:"); for(g in c("Miscarriage","Infection","PE","Preterm","Normal_Early","Normal_Late")) cat(sprintf(" %s=%d", g, sum(dc==g))); cat("\n")

# Run GSEA with t-test and relaxed cutoff
run_gsea <- function(disease, control, label) {
  cells <- colnames(seu)[dc %in% c(disease, control)]
  sub <- subset(seu, cells=cells); sub$group <- ifelse(dc[cells]==disease, "Disease", "Control")
  Idents(sub) <- "group"; sub <- JoinLayers(sub)
  deg <- FindMarkers(sub, ident.1="Disease", ident.2="Control", logfc.threshold=0, min.pct=0.1, test.use="t")
  deg$symbol <- rownames(deg)
  conv <- bitr(deg$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  deg_m <- inner_join(deg, conv, by=c("symbol"="SYMBOL"))
  gl <- setNames(deg_m$avg_log2FC, deg_m$ENTREZID); gl <- gl[!is.na(names(gl))]; gl <- gl[!duplicated(names(gl))]; gl <- sort(gl, decreasing=TRUE)
  # Relaxed p-value cutoff
  gse <- gseGO(geneList=gl, ont="BP", OrgDb=org.Hs.eg.db, pvalueCutoff=0.5, verbose=FALSE)
  if(nrow(as.data.frame(gse)) > 0) {
    gse <- setReadable(gse, OrgDb="org.Hs.eg.db", keyType="ENTREZID")
    res <- as.data.frame(gse); res$comparison <- label
    res$Count <- sapply(strsplit(as.character(res$core_enrichment), "/"), length)
    cat(sprintf("%s: %d terms (padj<0.1=%d)\n", label, nrow(res), sum(res$p.adjust<0.1)))
    return(res)
  }
  NULL
}

cat("Running boosted GSEA...\n")
combined <- rbind(
  run_gsea("Miscarriage","Normal_Early","Miscarriage"),
  run_gsea("Infection","Normal_Early","Infection"),
  run_gsea("PE","Normal_Late","PE"),
  run_gsea("Preterm","Normal_Late","Preterm"))
combined$comparison <- factor(combined$comparison, levels=c("Miscarriage","Infection","PE","Preterm"))
write.csv(combined, file.path(OUTDIR,"disease_GSEA_boosted.csv"), row.names=FALSE)

# Top 15 per disease by |NES|, prefer significant
top <- combined %>% filter(p.adjust < 0.1) %>%
  group_by(comparison) %>% slice_max(abs(NES), n=15) %>%
  pull(Description) %>% unique()
if(length(top) < 40) {
  extra <- combined %>% group_by(comparison) %>% slice_max(abs(NES), n=15) %>%
    pull(Description) %>% unique()
  top <- unique(c(top, extra))
}

plot_data <- combined %>% filter(Description %in% top) %>%
  mutate(Description=str_wrap(Description, width=50),
         comparison=factor(comparison, levels=c("Miscarriage","Infection","PE","Preterm")))
cat(sprintf("Plotting %d pathways\n", length(unique(plot_data$Description))))

# Order
desc_list <- unique(as.character(plot_data$Description))
func <- setNames(rep(5,length(desc_list)), desc_list)
func[grep("translat|ribosom|rna|spliceosom|cytoplas",desc_list)] <- 1
func[grep("mitochond|oxidative|atp|aerobic|respir|proton motive|electron",desc_list)] <- 2
func[grep("nuclear division|spin|chromosom|metaphase|mitotic|cell cycle|segregat",desc_list)] <- 3
func[grep("antigen|mhc|immune|t cell|leukocyt|cytotox|humoral|antimicrob|interferon|viral|symbion|lipopolysacch|bacterial|interleukin",desc_list)] <- 4
pg <- plot_data %>% group_by(Description) %>% summarise(g=first(func[Description]),.groups="drop") %>% arrange(g,Description)
plot_data$Description <- factor(plot_data$Description, levels=rev(pg$Description))

p <- ggplot(plot_data, aes(x=comparison, y=Description)) +
  geom_hline(yintercept=seq_along(levels(plot_data$Description)), color="grey92", linewidth=0.3) +
  geom_point(aes(size=Count, color=NES), stroke=0) +
  scale_color_gradient2(low="#313695",mid="#FFFFBF",high="#A50026",midpoint=0,name="NES",limits=c(-3.2,3.2)) +
  scale_size_continuous(range=c(2,7),name="Core\nGenes") + theme_bw() +
  labs(title="GO Biological Process — GSEA by Disease",x="",y="",
       subtitle="Disease vs trimester-matched Normal (t-test, PTNL→Preterm)") +
  theme(text=element_text(face="bold"),axis.text.x=element_text(size=12,color="black"),
    axis.text.y=element_text(size=9,color="black"),legend.title=element_text(size=10),
    legend.text=element_text(size=9),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    panel.border=element_rect(color="black",linewidth=0.6),
    plot.title=element_text(hjust=0.5,size=13),plot.subtitle=element_text(hjust=0.5,face="plain",size=10))

h <- max(8, length(unique(plot_data$Description))*0.25)
ggsave(file.path(FIGDIR,"Fig5c_GSEA_dotplot_boosted.png"), p, w=11, h=h, dpi=300, bg="white")
cat("Saved Fig5c_GSEA_dotplot_boosted.png\n")
