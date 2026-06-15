#!/usr/bin/env Rscript
# PE subgroup analysis: early vs late onset, severe vs control
library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
set.seed(42)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"
OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
detail <- labels$disease_detail[1:ncol(seu)]

# ── PE subtype assignment ──
# Early-onset PE: GSE290578 PE (GW 29-34)
# Late-onset PE: GSE298602 PreE_SF+gHTN + GSE173193 + GSE298119 (term)
# Severe: GSE298602 PreE_SF
# Control (within GSE298602): GSE298602 Control

seu$pe_subtype <- "Other"
seu$pe_subtype[seu$dataset=="GSE290578" & detail=="PE"] <- "Early_PE"
seu$pe_subtype[seu$dataset=="GSE298602" & detail=="PreE_SF"] <- "Late_Severe_PE"
seu$pe_subtype[seu$dataset=="GSE298602" & detail=="gHTN"] <- "Late_Mild_PE"
seu$pe_subtype[seu$dataset=="GSE298602" & detail=="Control"] <- "GSE298602_Control"
seu$pe_subtype[detail %in% c("GSE173193","GSE298119")] <- "Late_PE"
seu$pe_subtype[detail=="Normal"] <- "Normal_Late"

# ── 1. Subtype proportions ──
target <- c("Early_PE","Late_Severe_PE","Late_Mild_PE","Late_PE","Normal_Late","GSE298602_Control")
seu_sub <- subset(seu, pe_subtype %in% target)
seu_sub$pe_subtype <- factor(seu_sub$pe_subtype,
  levels=c("Early_PE","Late_Severe_PE","Late_Mild_PE","Late_PE","Normal_Late","GSE298602_Control"))

prop <- prop.table(table(seu_sub$subtype_pred, seu_sub$pe_subtype), margin=2)*100

st_cols <- c("Pro-inflammatory"="#C62828","MHCII+ Antigen-presenting"="#E65100",
  "Homeostatic"="#1565C0","PRKN+ Autophagy"="#6A1B9A",
  "Vascular remodeling"="#2E7D32","MKI67+ Proliferating"="#455A64")

prop_df <- as.data.frame(prop); colnames(prop_df) <- c("Subtype","Group","Proportion")

p1 <- ggplot(prop_df, aes(x=Group, y=Proportion, fill=Subtype)) +
  geom_col(width=0.65) + scale_fill_manual(values=st_cols) +
  labs(title="Subtype proportions: PE subgroups", y="Proportion (%)", x="") +
  theme_bw(11) + theme(panel.grid=element_blank(), plot.title=element_text(face="bold",size=13,hjust=0.5),
    axis.text.x=element_text(angle=30, hjust=1, size=9, color="black"),
    axis.text.y=element_text(size=9, color="black"), legend.text=element_text(size=8))
ggsave(file.path(FIGDIR,"PE_subtype_proportions.png"), p1, w=9, h=5, dpi=300, bg="white")

# Print proportions
cat("=== PE Subtype Proportions ===\n")
print(round(prop, 1))

# ── 2. Key gene expression: CEBPA, FLT1, FN1, PAPPA, HLA-DRA ──
key_genes <- c("CEBPA","FLT1","FN1","PAPPA","HLA-DRA","SPP1","IL1B")
seu_sub <- JoinLayers(seu_sub)
for(g in key_genes) {
  if(!g %in% rownames(seu_sub)) next
  expr <- FetchData(seu_sub, c(g, "pe_subtype"))
  means <- tapply(expr[,1], expr$pe_subtype, mean)
  cat(sprintf("\n%s means:\n", g))
  for(grp in levels(seu_sub$pe_subtype)) cat(sprintf("  %-25s %.3f\n", grp, means[grp]))
}

cat("\n=== Key: Early PE vs Late Severe PE ===\n")
early_genes <- seu_sub$pe_subtype=="Early_PE"
late_sev_genes <- seu_sub$pe_subtype=="Late_Severe_PE"
for(g in key_genes) {
  if(!g %in% rownames(seu_sub)) next
  expr <- FetchData(seu_sub, g)
  e <- mean(expr[early_genes,1]); l <- mean(expr[late_sev_genes,1])
  cat(sprintf("  %-10s  Early_PE=%.3f  Late_Severe_PE=%.3f  FC=%.2f\n", g, e, l, l/e))
}

# Save subset object
saveRDS(seu_sub, file.path(OUTDIR, "pe_subgroup_analysis.rds"))
cat("\nSaved: pe_subgroup_analysis.rds\n")
