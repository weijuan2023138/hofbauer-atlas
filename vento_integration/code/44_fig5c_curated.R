#!/usr/bin/env Rscript
# Fig5c: Curated GSEA dotplot — old 31 + new unique pathways, Fig1E format
library(Seurat); library(clusterProfiler); library(org.Hs.eg.db)
library(ggplot2); library(dplyr); library(stringr)
set.seed(42)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"

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

# ── Run GSEA ──
run_gsea <- function(disease, control, label) {
  cells <- colnames(seu)[seu$disease_clean %in% c(disease, control)]
  sub <- subset(seu, cells=cells); sub$group <- ifelse(sub$disease_clean==disease, "Disease", "Control")
  Idents(sub) <- "group"; sub <- JoinLayers(sub)
  deg <- FindMarkers(sub, ident.1="Disease", ident.2="Control", logfc.threshold=0, min.pct=0.1, test.use="wilcox")
  deg$symbol <- rownames(deg); conv <- bitr(deg$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  deg_m <- inner_join(deg, conv, by=c("symbol"="SYMBOL"))
  gl <- setNames(deg_m$avg_log2FC, deg_m$ENTREZID); gl <- gl[!is.na(names(gl))]; gl <- gl[!duplicated(names(gl))]; gl <- sort(gl, decreasing=TRUE)
  gse <- gseGO(geneList=gl, ont="BP", OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, verbose=FALSE)
  if(nrow(as.data.frame(gse)) > 0) {
    gse <- setReadable(gse, OrgDb="org.Hs.eg.db", keyType="ENTREZID")
    res <- as.data.frame(gse); res$comparison <- label
    res$Count <- sapply(strsplit(as.character(res$core_enrichment), "/"), length)
    return(res)
  }
  NULL
}

cat("Running 4 GSEAs...\n")
combined <- rbind(
  run_gsea("Miscarriage","Normal_Early","Miscarriage"),
  run_gsea("Infection","Normal_Early","Infection"),
  run_gsea("PE","Normal_Late","PE"),
  run_gsea("Preterm","Normal_Late","Preterm")
)
combined$comparison <- factor(combined$comparison, levels=c("Miscarriage","Infection","PE","Preterm"))
write.csv(combined, file.path(OUTDIR, "disease_GSEA_GOBP.csv"), row.names=FALSE)
cat(sprintf("Total results: %d\n", nrow(combined)))

# ── Old 31 pathways ──
old31 <- c("translation","t cell mediated immunity","t cell activation involved in immune response",
  "spliceosomal snrnp assembly","regulation of viral life cycle","regulation of t cell mediated cytotoxicity",
  "regulation of spindle checkpoint","regulation of mitotic nuclear division","proton motive force driven atp synthesis",
  "positive regulation of t cell mediated cytotoxicity","positive regulation of leukocyte cell cell adhesion",
  "positive regulation of cell cell adhesion","oxidative phosphorylation","negative regulation of viral process",
  "negative regulation of viral genome replication","multi multicellular organism process",
  "mitotic metaphase chromosome alignment","mitochondrial translation","mitochondrial respiratory chain complex assembly",
  "mitochondrial gene expression","cytoplasmic translation","atp synthesis coupled electron transport",
  "antimicrobial humoral response","antigen processing and presentation of peptide or polysaccharide antigen via mhc class ii",
  "antigen processing and presentation of peptide antigen","antigen processing and presentation of exogenous peptide antigen",
  "antigen processing and presentation of exogenous antigen","antigen processing and presentation of endogenous peptide antigen",
  "antigen processing and presentation of endogenous antigen","antigen processing and presentation","aerobic respiration")

# New unique pathways (padj<0.1, not in old31)
new_sig <- combined %>% filter(p.adjust < 0.1)
new_unique_desc <- setdiff(unique(new_sig$Description), old31)
# Keep top 10 by |NES|
new_top <- new_sig %>% filter(Description %in% new_unique_desc) %>%
  group_by(Description) %>% summarise(maxNES=max(abs(NES))) %>%
  arrange(-maxNES) %>% head(10) %>% pull(Description)

# Combine
all_curated <- unique(c(old31, new_top))
cat(sprintf("Curated pathways: %d (old=%d + new=%d)\n", length(all_curated), length(old31), length(new_top)))

plot_data <- combined %>% filter(Description %in% all_curated) %>%
  mutate(Description=str_wrap(Description, width=50),
         comparison=factor(comparison, levels=c("Miscarriage","Infection","PE","Preterm")))

# Order: functional groups
func_group <- setNames(rep(1, length(old31)), old31)
func_group[grep("translation|ribosom|rna|spliceosom|cytoplasmic", old31)] <- 1  # Translation
func_group[grep("mitochond|oxidative|atp|aerobic|respir|proton", old31)] <- 2  # Mitochondria
func_group[grep("cell cycle|mitotic|spindle|chromosom|nuclear division", old31)] <- 3  # Cell cycle
func_group[grep("antigen|mhc|immune|t cell|leukocyte|cytotox|humoral|antimicrobial", old31)] <- 4  # Immune
func_group[grep("viral|multi multicell", old31)] <- 5  # Other
# New pathways → group 4 (immune-related mostly)
for(p in new_top) func_group[p] <- 6

pg <- plot_data %>% group_by(Description) %>% summarise(g=first(func_group[Description]), .groups="drop") %>% arrange(g, Description)
plot_data$Description <- factor(plot_data$Description, levels=rev(pg$Description))

# Plot
p <- ggplot(plot_data, aes(x=comparison, y=Description)) +
  geom_hline(yintercept=seq_along(levels(plot_data$Description)), color="grey92", linewidth=0.3) +
  geom_point(aes(size=Count, color=NES), stroke=0) +
  scale_color_gradient2(low="#313695", mid="#FFFFBF", high="#A50026", midpoint=0, name="NES", limits=c(-3.2,3.2)) +
  scale_size_continuous(range=c(2, 7), name="Core\nGenes") + theme_bw() +
  labs(title="GO Biological Process — GSEA by Disease", x="", y="",
       subtitle="Disease vs trimester-matched Normal") +
  theme(text=element_text(face="bold"), axis.text.x=element_text(size=12,color="black"),
    axis.text.y=element_text(size=9,color="black"), legend.title=element_text(size=10),
    legend.text=element_text(size=9), panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(), panel.border=element_rect(color="black",linewidth=0.6),
    plot.title=element_text(hjust=0.5,size=13), plot.subtitle=element_text(hjust=0.5,face="plain",size=10))

ggsave(file.path(FIGDIR,"Fig5c_GSEA_dotplot.png"), p, w=11, h=12, dpi=300, bg="white")
cat("Saved Fig5c_GSEA_dotplot.png\n")
