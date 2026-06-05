#!/usr/bin/env Rscript
# TF activity: differential TF expression Early vs Late
# Using curated human TF list from Lambert et al. 2018
library(Seurat); library(dplyr); library(ggplot2)

seu <- readRDS('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/Hofbauer_Atlas_Final.rds')
cond <- read.csv("/tmp/gse290578_conditions.csv", stringsAsFactors=FALSE)
rownames(cond) <- cond$barcode; seu$gse <- NA
for(bc in colnames(seu)) if(bc %in% cond$barcode) seu$gse[bc] <- cond[bc,"condition"]
seu$tri <- NA
seu$tri[seu$dataset=="Arutyunyan"] <- "Early"
seu$tri[seu$dataset=="GSE290578" & seu$gse=="Normal"] <- "Late"
seu$tri <- factor(seu$tri, levels=c("Early","Late"))
seu_normal <- subset(seu, cells=colnames(seu)[!is.na(seu$tri)])

# Load full DEG results
deg <- read.csv('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/dev_trimester_DEGs.csv')
early_deg <- deg[deg$trimester=="Early",]; late_deg <- deg[deg$trimester=="Late",]

# ---- Human TFs (Lambert et al. 2018 + GO:0003700) ----
# Use a manually curated list of well-known TFs
human_tfs <- c(
  # bZIP
  "FOS","FOSB","FOSL1","FOSL2","JUN","JUNB","JUND","ATF1","ATF2","ATF3","ATF4","ATF6",
  "BATF","BATF3","CEBPA","CEBPB","CEBPD","CEBPE","CEBPG","CREB1","CREB3","MAF","MAFB",
  "MAFF","MAFG","MAFK","NFE2","NFE2L1","NFE2L2","XBP1",
  # bHLH
  "MYC","MYCN","MYCL","MAX","MXD1","MXD3","MXD4","MXI1","TFEB","TFE3","TFEC","MITF",
  "USF1","USF2","ID1","ID2","ID3","ID4","HES1","HEY1","ARNT","ARNT2","AHR",
  # Homeodomain
  "HOXA5","HOXA9","HOXA10","HOXB2","HOXB7","HOXC6","HOXD4","PBX1","PBX3","MEIS1","MEIS2",
  "CDX2","LHX2","ISL1","EN1","OTX2","PITX1","NKX2-1","NKX2-3","NKX3-1","SIX1","SIX2",
  # Forkhead
  "FOXA1","FOXA2","FOXA3","FOXC1","FOXC2","FOXD1","FOXF1","FOXK2","FOXL2","FOXM1",
  "FOXN3","FOXO1","FOXO3","FOXO4","FOXP1","FOXP2","FOXP3",
  # ETS
  "ETS1","ETS2","ELK1","ELK3","ELK4","ELF1","ELF2","ELF3","ELF4","ELF5",
  "ERG","FLI1","ETV1","ETV4","ETV5","ETV6","SPI1","SPIB","GABPA",
  # IRF
  "IRF1","IRF2","IRF3","IRF4","IRF5","IRF6","IRF7","IRF8","IRF9",
  # NF-kB
  "NFKB1","NFKB2","RELA","RELB","REL","NFKBIA","NFKBIZ",
  # STAT
  "STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6",
  # Nuclear receptors
  "NR3C1","NR4A1","NR4A2","NR4A3","PPARG","PPARA","PPARD","RXRA","RXRB","VDR",
  "ESR1","ESR2","AR","RARA","RARB","RARG","RORA","RORC","NR1H3","NR1H4",
  # SMAD
  "SMAD1","SMAD2","SMAD3","SMAD4","SMAD5","SMAD7",
  # Others
  "TP53","TP63","TP73","E2F1","E2F2","E2F3","E2F4","E2F6","E2F7",
  "KLF4","KLF5","KLF6","KLF7","KLF10","KLF13","SP1","SP2","SP3","SP4",
  "GATA1","GATA2","GATA3","GATA4","GATA6",
  "RUNX1","RUNX2","RUNX3","TBX21","EOMES","TBX2","TBX3",
  "SOX2","SOX4","SOX5","SOX6","SOX9","SOX10","SOX11","SOX17",
  "TCF3","TCF4","TCF7","TCF7L1","TCF7L2","LEF1",
  "HIF1A","HIF2A","SREBF1","SREBF2","MEF2A","MEF2C","MEF2D",
  "NOTCH1","NOTCH2","NOTCH3","RBPJ","CTCF","YY1","YY2","EGR1","EGR2","EGR3",
  "ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","TWIST2","PRDM1","BCL6","BCL6B",
  "HMGA1","HMGA2","PAX5","PAX6","WT1","TEAD1","TEAD4","TAZ","YAP1"
)

# Find TFs that are DEG
early_tfs <- early_deg[early_deg$gene %in% human_tfs,]
late_tfs  <- late_deg[late_deg$gene %in% human_tfs,]

# Combine and pick top by |avg_log2FC|
all_tfs <- rbind(
  data.frame(TF=early_tfs$gene, logFC=early_tfs$avg_log2FC, p_adj=early_tfs$p_val_adj, Direction="Early-up"),
  data.frame(TF=late_tfs$gene, logFC=late_tfs$avg_log2FC, p_adj=late_tfs$p_val_adj, Direction="Late-up")
)
all_tfs <- all_tfs[!duplicated(all_tfs$TF),]
all_tfs <- all_tfs[order(-abs(all_tfs$logFC)),]

# Keep top TFs (p_adj < 0.05 and |logFC| > 0.3)
top_tfs <- all_tfs %>% filter(p_adj < 0.05, abs(logFC) > 0.3)
cat(sprintf("Significant TFs: %d (Early-up=%d, Late-up=%d)\n",
  nrow(top_tfs), sum(top_tfs$Direction=="Early-up"), sum(top_tfs$Direction=="Late-up")))

# Show top
cat("\nTop Early-up TFs:\n")
print(head(top_tfs[top_tfs$Direction=="Early-up", c("TF","logFC","p_adj")], 20))
cat("\nTop Late-up TFs:\n")
print(head(top_tfs[top_tfs$Direction=="Late-up", c("TF","logFC","p_adj")], 20))

# ---- Plot: TF dot plot ----
plot_data <- head(top_tfs[order(-abs(top_tfs$logFC)),], 30)
plot_data$TF <- factor(plot_data$TF, levels=plot_data$TF[order(plot_data$logFC)])

p <- ggplot(plot_data, aes(x=logFC, y=TF, color=Direction)) +
  geom_vline(xintercept=0, color="grey70", linewidth=0.4) +
  geom_point(aes(size=-log10(p_adj)), stroke=0) +
  scale_color_manual(values=c("Early-up"="#4575B4","Late-up"="#D73027")) +
  scale_size_continuous(range=c(2,7), name="-log10(P.adj)") +
  labs(x="log2 Fold Change (Trimester-specific)", y="",
       title="Differential Transcription Factors",
       subtitle="Early (GW4.5\u201310) vs Late (GW32\u201338)") +
  theme_bw() +
  theme(
    axis.text=element_text(color="black", size=9),
    axis.text.y=element_text(face="bold.italic", size=9),
    legend.position="right",
    legend.title=element_text(size=9), legend.text=element_text(size=8),
    panel.grid.major=element_line(color="grey92", linewidth=0.3),
    panel.grid.minor=element_blank(),
    panel.border=element_rect(color="black", linewidth=0.5),
    plot.title=element_text(hjust=0.5, face="bold", size=12),
    plot.subtitle=element_text(hjust=0.5, size=9)
  )

FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'
ggsave(file.path(FIGDIR,"Fig1_TF_dotplot.png"), p, w=6.5, h=8, dpi=300, bg="white")
write.csv(all_tfs, file.path(FIGDIR, "../results/dev_TF_DEGs.csv"), row.names=FALSE)
cat("\nSaved Fig1_TF_dotplot.png\n")
