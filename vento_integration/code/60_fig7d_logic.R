#!/usr/bin/env Rscript
# Fig7d: TF-regulon coupling slope via AddModuleScore — 6-group + PE 3-group
library(Seurat); library(ggplot2); library(dplyr)
set.seed(42)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
detail <- labels$disease_detail[1:ncol(seu)]

# ── Disease labels ──
dc <- rep(NA_character_, ncol(seu))
dc[seu$dataset %in% c("E-MTAB-12421","E-MTAB-6701") | (seu$dataset=="E-MTAB-12795" & detail=="normal")] <- "Normal Early"
dc[seu$dataset=="GSE214607"] <- "Miscarriage"
dc[detail %in% c("toxoplasmosis","listeriosis","Plasmodium malariae malaria")] <- "Infection"
dc[(seu$dataset=="GSE290578" & detail=="Normal") | (seu$dataset=="GSE298602" & detail=="Control")] <- "Normal Late"
dc[detail %in% c("PE","PreE_SF","gHTN","GSE173193","GSE298119")] <- "PE"
dc[detail %in% c("PTL","PTNL")] <- "Preterm"

# ── PE 3-group labels ──
dc3 <- rep(NA_character_, ncol(seu))
dc3[(seu$dataset=="GSE290578" & detail=="Normal") | (seu$dataset=="GSE298602" & detail=="Control")] <- "Normal Late"
dc3[seu$dataset=="GSE290578" & detail=="PE"] <- "Early PE"
dc3[detail %in% c("PreE_SF","gHTN","GSE173193","GSE298119")] <- "Late PE"

# ── Regulons ──
regulons <- list(
  CEBPA = c("SPP1","FN1","COL1A2","PAPPA","FLT1","CD44","AEBP1","NOTUM","COL1A1","IGF1"),
  IRF1  = c("HLA-DRA","CD74","STAT1","IFITM1","CXCL10","GBP2","FCGR3A","IFI27","HLA-DQB1","HLA-DPA1","IRF8"),
  KLF4  = c("TGFB1","VEGFA","THBS1","ITGAV","PTPRM","MMP9","BMP2","COL4A1","CD47","ITGB1","ITGB5")
)
for(nm in names(regulons)) regulons[[nm]] <- intersect(regulons[[nm]], rownames(seu))

# Add module scores
seu <- AddModuleScore(seu, features=regulons, name="regulon_", ctrl=min(100, nrow(seu)/2))
# Rename module columns to TF_AUC
for(i in seq_along(regulons)) {
  colnames(seu@meta.data)[colnames(seu@meta.data) == paste0("regulon_", i)] <- paste0(names(regulons)[i], "_AUC")
}
# Add TF expression to metadata
for(tf in names(regulons)) {
  seu@meta.data[[tf]] <- FetchData(seu, vars=tf)[,1]
}

compute_slopes <- function(meta, groups_vec, target_groups) {
  results <- data.frame()
  for(tf in names(regulons)) {
    auc_col <- paste0(tf, "_AUC")
    for(g in target_groups) {
      mask <- groups_vec == g & !is.na(groups_vec)
      if(sum(mask) < 50) next
      r <- cor(meta[[tf]][mask], meta[[auc_col]][mask], method="spearman")
      results <- rbind(results, data.frame(Disease=g, TF=tf, Slope=r))
    }
  }
  results$TF <- factor(results$TF, levels=names(regulons))
  results$Disease <- factor(results$Disease, levels=target_groups)
  results
}

# ── Figure 1: 6 groups ──
groups6 <- c("Normal Early","Miscarriage","Infection","Normal Late","PE","Preterm")
slope6 <- compute_slopes(seu@meta.data, dc, groups6)
cat("\n=== 6-group slopes ===\n")
for(i in 1:nrow(slope6)) cat(sprintf("  %-10s %-6s r=%.3f\n", slope6$Disease[i], slope6$TF[i], slope6$Slope[i]))

p1 <- ggplot(slope6, aes(x=Disease, y=Slope, fill=TF)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.7, color="black", linewidth=0.3) +
  scale_fill_manual(values=c("CEBPA"="#4A7BB0","IRF1"="#D93829","KLF4"="#7B3294")) +
  geom_hline(yintercept=0, color="grey50", linewidth=0.3) +
  labs(title="Regulatory coupling strength", y="TF-regulon coupling (Spearman r)", x="") +
  theme_bw(12) + theme(axis.text.x=element_text(angle=30,hjust=1,size=11,color="black",face="bold"),
    axis.text.y=element_text(size=10,color="black"), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
    plot.title=element_text(hjust=0.5,size=14,face="bold"),
    legend.position="top",legend.title=element_blank(),legend.text=element_text(size=12))
ggsave(file.path(FIGDIR,"Fig7d_slope_6groups.png"), p1, w=10, h=5, dpi=300, bg="white")

# ── Figure 2: PE 3-group ──
groups3 <- c("Normal Late","Early PE","Late PE")
slope3 <- compute_slopes(seu@meta.data, dc3, groups3)
cat("\n=== PE 3-group slopes ===\n")
for(i in 1:nrow(slope3)) cat(sprintf("  %-10s %-6s r=%.3f\n", slope3$Disease[i], slope3$TF[i], slope3$Slope[i]))

p2 <- ggplot(slope3, aes(x=Disease, y=Slope, fill=TF)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.7, color="black", linewidth=0.3) +
  scale_fill_manual(values=c("CEBPA"="#4A7BB0","IRF1"="#D93829","KLF4"="#7B3294")) +
  geom_hline(yintercept=0, color="grey50", linewidth=0.3) +
  labs(title="Regulatory coupling — PE subtypes", y="TF-regulon coupling (Spearman r)", x="") +
  theme_bw(12) + theme(axis.text.x=element_text(angle=0,hjust=0.5,size=12,color="black",face="bold"),
    axis.text.y=element_text(size=10,color="black"), panel.grid.major.x=element_blank(), panel.grid.minor=element_blank(),
    plot.title=element_text(hjust=0.5,size=14,face="bold"),
    legend.position="top",legend.title=element_blank(),legend.text=element_text(size=12))
ggsave(file.path(FIGDIR,"Fig7d_slope_PE3.png"), p2, w=6, h=5, dpi=300, bg="white")

cat("\nBoth figures saved\n")
