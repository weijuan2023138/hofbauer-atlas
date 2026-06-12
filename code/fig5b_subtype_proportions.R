#!/usr/bin/env Rscript
# Fig5b: Subtype proportions stacked bar chart (6 disease groups)
library(Seurat); library(ggplot2); library(dplyr)

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

st_cols <- c("Pro-inflammatory"="#C62828","MHCII+ Antigen-presenting"="#E65100",
             "Homeostatic"="#1565C0","PRKN+ Autophagy"="#6A1B9A",
             "Vascular remodeling"="#2E7D32","MKI67+ Proliferating"="#455A64")
seu_sub$disease_clean <- factor(seu_sub$disease_clean,
  levels=c("Normal_Early","Miscarriage","Infection","Normal_Late","PE","Preterm"))
prop <- prop.table(table(seu_sub$subtype, seu_sub$disease_clean), margin=2)
prop_df <- as.data.frame(prop); colnames(prop_df) <- c("Subtype","Disease","Proportion")

# with TL (TL merged into Normal_Late)
pb <- ggplot(prop_df, aes(x=Disease, y=Proportion*100, fill=Subtype)) +
  geom_col(width=0.65) + scale_fill_manual(values=st_cols) +
  labs(title="Subtype proportions", y="Proportion (%)", x="") +
  theme_bw(base_size=11) +
  theme(panel.grid=element_blank(), plot.title=element_text(face="bold",size=13,hjust=0.5,color="black"),
        legend.position="right",
        axis.text.x=element_text(face="bold",size=9,color="black"),
        axis.text.y=element_text(size=9,color="black"),
        axis.title.y=element_text(size=10,color="black"),
        legend.text=element_text(size=9,color="black"),
        legend.title=element_text(size=10,face="bold",color="black"))
ggsave("figures/Fig5/Fig5b_subtype_proportions.png", pb, w=7, h=5, dpi=300, bg="white")

# noTL (exclude TL from Normal_Late)
keep2 <- keep & seu$disease != "TL"
seu_sub2 <- subset(seu, cells=colnames(seu)[keep2])
seu_sub2$disease_clean <- factor(seu_sub2$disease_clean,
  levels=c("Normal_Early","Miscarriage","Infection","Normal_Late","PE","Preterm"))
prop2 <- prop.table(table(seu_sub2$subtype, seu_sub2$disease_clean), margin=2)
prop_df2 <- as.data.frame(prop2); colnames(prop_df2) <- c("Subtype","Disease","Proportion")
pb2 <- ggplot(prop_df2, aes(x=Disease, y=Proportion*100, fill=Subtype)) +
  geom_col(width=0.65) + scale_fill_manual(values=st_cols) +
  labs(title="Subtype proportions", y="Proportion (%)", x="") +
  theme_bw(base_size=11) +
  theme(panel.grid=element_blank(), plot.title=element_text(face="bold",size=13,hjust=0.5,color="black"),
        legend.position="right",
        axis.text.x=element_text(face="bold",size=9,color="black"),
        axis.text.y=element_text(size=9,color="black"),
        axis.title.y=element_text(size=10,color="black"),
        legend.text=element_text(size=9,color="black"),
        legend.title=element_text(size=10,face="bold",color="black"))
ggsave("figures/Fig5/Fig5b_subtype_proportions_noTL.png", pb2, w=7, h=5, dpi=300, bg="white")
cat("Fig5b done\n")
