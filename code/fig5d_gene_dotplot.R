#!/usr/bin/env Rscript
# Fig5d: Gene Z-score dotplot — 40 genes across 6 disease groups
library(Seurat); library(ggplot2); library(dplyr); library(tidyr)

seu <- readRDS("results/Hofbauer_Atlas_Final.rds")
mask_gse <- seu$dataset == "GSE290578"; bcs <- colnames(seu)
is_norm_gse <- grepl("_Norm_", bcs) & mask_gse; is_pt_gse <- grepl("_Pt_", bcs) & mask_gse
seu$disease_clean <- NA
seu$disease_clean[seu$disease == "Normal_1st"] <- "Normal_Early"
seu$disease_clean[seu$disease == "Normal"] <- "Normal_Late"
seu$disease_clean[seu$disease == "RM/NC"] <- "Miscarriage"
seu$disease_clean[seu$disease == "Normal/Listeria/Toxoplasma/Malaria"] <- "Infection"
seu$disease_clean[is_norm_gse] <- "Normal_Late"; seu$disease_clean[is_pt_gse] <- "PE"
seu$disease_clean[seu$disease %in% c("PTL","PTNL")] <- "Preterm"
seu$disease_clean[seu$disease == "TL"] <- "Normal_Late"
keep <- seu$disease_clean %in% c("Normal_Early","Normal_Late","Miscarriage","PE","Infection","Preterm")
seu_sub <- subset(seu, cells=colnames(seu)[keep])

genes <- c(
  "SQSTM1","BNIP3","PRKN","TREM2","RELB","NFKB1","STAT3",
  "STAT1","ID2","MAFB",
  "CEBPA","IFI27",
  "FCGR2B","HLA-DQA1","HLA-DQB1",
  "NNMT","HLA-G",
  "CFI","AEBP1","HIGD2A","COX14","ATP5MF",
  "PRG2","PSG9","TNF","IL1B","CXCL8","HLA-DRA","FCGR3A","C1QC","C1QB","C1QA","ITGB1",
  "ITGAV","CD44",
  "PTPRM","TGFB1","COL1A2","FN1",
  "SPP1"
)
genes <- intersect(genes, rownames(seu_sub))

groups <- c("Normal_Early","Miscarriage","Infection","Normal_Late","PE","Preterm")
group_labels <- c("Normal\n(Early)","Miscarriage","Infection","Normal\n(Late)","PE","Preterm")

expr <- FetchData(seu_sub, vars=c(genes, "disease_clean"))
plot_data <- data.frame()
for(g in genes) {
  grp_means <- sapply(groups, function(gr) mean(expr[expr$disease_clean==gr, g]))
  z <- (grp_means - mean(grp_means)) / sd(grp_means)
  for(i in seq_along(groups)) {
    plot_data <- rbind(plot_data, data.frame(
      Gene=g, Group=groups[i], GroupLabel=group_labels[i],
      Zscore=z[i], AbsZ=abs(z[i])))
  }
}
plot_data$Gene <- factor(plot_data$Gene, levels=genes)
plot_data$GroupLabel <- factor(plot_data$GroupLabel, levels=group_labels)

# Horizontal version
ph <- ggplot(plot_data, aes(x=Gene, y=GroupLabel)) +
  geom_point(aes(fill=Zscore, size=AbsZ), shape=21, stroke=0.3, color="grey80") +
  scale_fill_gradient2(low="#313695", mid="#FFFFBF", high="#A50026", midpoint=0, name="Z-score") +
  scale_size_continuous(range=c(2, 8), name="|Z-score|", breaks=c(0.5,1,1.5,2)) +
  labs(title="Key gene expression", x="", y="") +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9, color="black"),
        axis.text.y=element_text(size=10, color="black"),
        axis.line=element_line(color="black", linewidth=0.4),
        panel.grid.major=element_line(color="grey92", linewidth=0.3),
        panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold", size=14, hjust=0.5),
        legend.position="right",
        legend.title=element_text(size=10, face="bold"),
        legend.text=element_text(size=9))
ggsave("figures/Fig5/Fig5d_gene_dotplot.png", ph, w=12, h=4.5, dpi=300, bg="white")

# Vertical version
pv <- ggplot(plot_data, aes(x=GroupLabel, y=Gene)) +
  geom_point(aes(fill=Zscore, size=AbsZ), shape=21, stroke=0.3, color="grey80") +
  scale_fill_gradient2(low="#313695", mid="#FFFFBF", high="#A50026", midpoint=0, name="Z-score") +
  scale_size_continuous(range=c(2, 8), name="|Z-score|", breaks=c(0.5,1,1.5,2)) +
  labs(title="Key gene expression", x="", y="") +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(size=10, color="black"),
        axis.text.y=element_text(size=9, color="black", face="italic"),
        axis.line=element_line(color="black", linewidth=0.4),
        panel.grid.major=element_line(color="grey92", linewidth=0.3),
        panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold", size=14, hjust=0.5),
        legend.position="right",
        legend.title=element_text(size=10, face="bold"),
        legend.text=element_text(size=9))
ggsave("figures/Fig5/Fig5d_gene_dotplot_vertical.png", pv, w=4.5, h=10, dpi=300, bg="white")
cat("Fig5d done\n")
