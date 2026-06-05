#!/usr/bin/env Rscript
# Fig4 Disease Analysis: PE, Miscarriage, Infection vs Normal
library(Seurat); library(ggplot2); library(dplyr); library(patchwork)

seu <- readRDS('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/Hofbauer_Atlas_Final.rds')
cond <- read.csv("/tmp/gse290578_conditions.csv", stringsAsFactors=FALSE)
rownames(cond) <- cond$barcode
seu$gse_cond <- NA
for(bc in colnames(seu)) if(bc %in% cond$barcode) seu$gse_cond[bc] <- cond[bc,"condition"]

# ---- Build clean disease groups ----
seu$disease_clean <- as.character(seu$disease)
# Split Normal/PE into Normal and PE
for(i in 1:ncol(seu)) {
  if(seu$disease[i] == "Normal/PE" && !is.na(seu$gse_cond[i])) {
    seu$disease_clean[i] <- ifelse(seu$gse_cond[i] == "PE", "PE", "Normal")
  }
}
# Merge Normal + Normal_1st
seu$disease_clean[seu$disease_clean %in% c("Normal","Normal_1st")] <- "Normal"
# Rename RM/NC to Miscarriage
seu$disease_clean[seu$disease_clean == "RM/NC"] <- "Miscarriage"
# Rename infection
seu$disease_clean[seu$disease_clean == "Normal/Listeria/Toxoplasma/Malaria"] <- "Infection"

cat("Disease groups:\n"); print(table(seu$disease_clean))

# Focus on PE, Miscarriage, Infection vs Normal
target <- c("Normal","PE","Miscarriage","Infection")
seu_disease <- subset(seu, disease_clean %in% target)
seu_disease$disease_clean <- factor(seu_disease$disease_clean, levels=target)

# ---- Subtype proportions by disease ----
prop <- prop.table(table(seu_disease$subtype, seu_disease$disease_clean), margin=2)
prop_df <- as.data.frame(prop)
colnames(prop_df) <- c("Subtype","Disease","Proportion")

st_cols <- c("Pro-inflammatory"="#C62828","MHCII+ Antigen-presenting"="#BF4E1A",
  "Homeostatic"="#1B6B93","PRKN+ Autophagy"="#7B4FA0",
  "Vascular remodeling"="#2D8B57","MKI67+ Proliferating"="#37474F")

p_prop <- ggplot(prop_df, aes(x=Disease, y=Proportion, fill=Subtype)) +
  geom_col(width=0.7) + scale_fill_manual(values=st_cols) +
  labs(title="Subtype Proportions by Disease", y="Proportion") +
  theme_bw() + theme(panel.grid=element_blank())

FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'
ggsave(file.path(FIGDIR,"Fig4a_subtype_proportions.png"), p_prop, w=8, h=5, dpi=200)

# ---- DEG: PE vs Normal ----
Idents(seu_disease) <- "disease_clean"
pe_markers <- FindMarkers(seu_disease, ident.1="PE", ident.2="Normal", logfc.threshold=0.25)
pe_markers$gene <- rownames(pe_markers)
write.csv(pe_markers, file.path(FIGDIR, "../results/deg_PE_vs_Normal.csv"))

# Volcano PE
pe_markers$sig <- ifelse(pe_markers$p_val_adj<0.05 & abs(pe_markers$avg_log2FC)>0.5, 
  ifelse(pe_markers$avg_log2FC>0, "PE-up","Normal-up"), "NS")
key_pe <- c("NFKB1","RELB","TNF","IL1B","CCL8","CXCL8","NLRP3","FCGR3A")
pe_markers$label <- ifelse(pe_markers$gene %in% key_pe & pe_markers$sig!="NS", pe_markers$gene, "")

p_pe <- ggplot(pe_markers, aes(x=avg_log2FC, y=-log10(p_val_adj), color=sig)) +
  geom_point(size=0.5, alpha=0.5) + ggrepel::geom_text_repel(aes(label=label), size=3, max.overlaps=15) +
  scale_color_manual(values=c("PE-up"="#D73027","Normal-up"="#4575B4","NS"="grey80")) +
  geom_vline(xintercept=c(-0.5,0.5), lty="dashed", color="grey50", lwd=0.3) +
  labs(title="PE vs Normal", x="log2FC", y="-log10(P.adj)") +
  theme_bw() + theme(panel.grid=element_blank(), legend.position="right", legend.title=element_blank())
ggsave(file.path(FIGDIR,"Fig4b_PE_volcano.png"), p_pe, w=7, h=6, dpi=200)

# ---- DEG: Miscarriage vs Normal ----
mis_markers <- FindMarkers(seu_disease, ident.1="Miscarriage", ident.2="Normal", logfc.threshold=0.25)
mis_markers$gene <- rownames(mis_markers)
write.csv(mis_markers, file.path(FIGDIR, "../results/deg_Miscarriage_vs_Normal.csv"))

mis_markers$sig <- ifelse(mis_markers$p_val_adj<0.05 & abs(mis_markers$avg_log2FC)>0.5, 
  ifelse(mis_markers$avg_log2FC>0, "Misc-up","Normal-up"), "NS")
p_mis <- ggplot(mis_markers, aes(x=avg_log2FC, y=-log10(p_val_adj), color=sig)) +
  geom_point(size=0.5, alpha=0.5) +
  scale_color_manual(values=c("Misc-up"="#D73027","Normal-up"="#4575B4","NS"="grey80")) +
  geom_vline(xintercept=c(-0.5,0.5), lty="dashed", color="grey50", lwd=0.3) +
  labs(title="Miscarriage vs Normal", x="log2FC", y="-log10(P.adj)") +
  theme_bw() + theme(panel.grid=element_blank(), legend.position="right", legend.title=element_blank())
ggsave(file.path(FIGDIR,"Fig4c_Miscarriage_volcano.png"), p_mis, w=7, h=6, dpi=200)

# ---- DEG: Infection vs Normal ----
inf_markers <- FindMarkers(seu_disease, ident.1="Infection", ident.2="Normal", logfc.threshold=0.25)
inf_markers$gene <- rownames(inf_markers)
write.csv(inf_markers, file.path(FIGDIR, "../results/deg_Infection_vs_Normal.csv"))

cat(sprintf("DEG: PE=%d, Misc=%d, Infect=%d\n", 
  sum(pe_markers$sig!="NS"), sum(mis_markers$sig!="NS"), sum(inf_markers$p_val_adj<0.05 & abs(inf_markers$avg_log2FC)>0.5)))
cat("Saved Fig4a-c\n")
