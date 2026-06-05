#!/usr/bin/env Rscript
# Fig5: Disease analysis — use pre-computed DEGs
library(Seurat); library(ggplot2); library(dplyr); library(ggrepel)

FIGDIR <- "figures/Fig5"; dir.create(FIGDIR, showWarnings=FALSE, recursive=TRUE)

seu <- readRDS("results/Hofbauer_Atlas_Final.rds")
cond <- read.csv("/tmp/gse290578_conditions.csv", stringsAsFactors=FALSE)
rownames(cond) <- cond$barcode
seu$gse_cond <- NA
for(bc in colnames(seu)) if(bc %in% cond$barcode) seu$gse_cond[bc] <- cond[bc,"condition"]
seu$disease_clean <- as.character(seu$disease)
for(i in 1:ncol(seu)) {
  if(seu$disease[i] == "Normal/PE" && !is.na(seu$gse_cond[i]))
    seu$disease_clean[i] <- ifelse(seu$gse_cond[i] == "PE", "PE", "Normal")
}
seu$disease_clean[seu$disease_clean %in% c("Normal","Normal_1st")] <- "Normal"
seu$disease_clean[seu$disease_clean == "RM/NC"] <- "Miscarriage"
seu$disease_clean[seu$disease_clean == "Normal/Listeria/Toxoplasma/Malaria"] <- "Infection"
target <- c("Normal","PE","Miscarriage","Infection")
seu_disease <- subset(seu, disease_clean %in% target)
seu_disease$disease_clean <- factor(seu_disease$disease_clean, levels=target)

st_cols <- c("Pro-inflammatory"="#C62828","MHCII+ Antigen-presenting"="#BF4E1A",
  "Homeostatic"="#1B6B93","PRKN+ Autophagy"="#7B4FA0",
  "Vascular remodeling"="#2D8B57","MKI67+ Proliferating"="#37474F")

# ── Panel A: Subtype proportions ──
prop <- prop.table(table(seu_disease$subtype, seu_disease$disease_clean), margin=2)
prop_df <- as.data.frame(prop); colnames(prop_df) <- c("Subtype","Disease","Proportion")
pA <- ggplot(prop_df, aes(x=Disease, y=Proportion, fill=Subtype)) +
  geom_col(width=0.7) + scale_fill_manual(values=st_cols) +
  labs(title="Subtype Proportions by Disease") +
  theme_bw() + theme(panel.grid=element_blank(), legend.position="right")
ggsave(file.path(FIGDIR,"Fig5a_subtype_proportions.png"), pA, w=8, h=5, dpi=300)

# ── Volcano function ──
make_volcano <- function(deg_file, title, up_label, key_genes) {
  m <- read.csv(deg_file, row.names=1)
  m$gene <- rownames(m)
  m$sig <- ifelse(m$p_val_adj<0.05 & abs(m$avg_log2FC)>0.5,
    ifelse(m$avg_log2FC>0, up_label, "Normal-up"), "NS")
  m$label <- ifelse(m$gene %in% key_genes & m$sig!="NS", m$gene, "")
  n_up <- sum(m$sig==up_label); n_dn <- sum(m$sig=="Normal-up")
  ggplot(m, aes(x=avg_log2FC, y=-log10(p_val_adj), color=sig)) +
    geom_point(size=0.4, alpha=0.5) +
    geom_text_repel(aes(label=label), size=3, max.overlaps=20, box.padding=0.3) +
    scale_color_manual(values=setNames(c("#D73027","#4575B4","grey80"), c(up_label,"Normal-up","NS"))) +
    geom_vline(xintercept=c(-0.5,0.5), lty="dashed", color="grey50", lwd=0.3) +
    labs(title=sprintf("%s (%d up, %d down)", title, n_up, n_dn), x="log2FC", y="-log10(P.adj)") +
    theme_bw() + theme(panel.grid=element_blank(), legend.position="none")
}

# Panel B-D
key_pe <- c("NFKB1","RELB","TNF","IL1B","CCL8","CXCL8","NLRP3","FCGR3A","CD36","TREM2","SOD2","C1QA")
pB <- make_volcano("results/deg_PE_vs_Normal.csv", "PE vs Normal", "PE-up", key_pe)
ggsave(file.path(FIGDIR,"Fig5b_PE_volcano.png"), pB, w=7, h=6, dpi=300)

key_mis <- c("NFKB1","TNF","IL1B","CXCL8","HLA-DRA","FCGR3A","C1QA","CD36")
pC <- make_volcano("results/deg_Miscarriage_vs_Normal.csv", "Miscarriage vs Normal", "Misc-up", key_mis)
ggsave(file.path(FIGDIR,"Fig5c_Miscarriage_volcano.png"), pC, w=7, h=6, dpi=300)

key_inf <- c("NFKB1","RELB","TNF","IL1B","CXCL8","CCL8","NLRP3","IFI30","HLA-DRA","CD74")
pD <- make_volcano("results/deg_Infection_vs_Normal.csv", "Infection vs Normal", "Inf-up", key_inf)
ggsave(file.path(FIGDIR,"Fig5d_Infection_volcano.png"), pD, w=7, h=6, dpi=300)

cat("Fig5 a-d saved (300 dpi)\n")
