#!/usr/bin/env Rscript
# Fig5S1: UpSet plot — DEG overlap across 4 diseases
library(UpSetR)

mis <- read.csv("results/deg_Miscarriage_vs_NormalEarly.csv", row.names=1)
inf <- read.csv("results/deg_Infection_vs_NormalEarly.csv", row.names=1)
pe  <- read.csv("results/deg_PE_vs_NormalLate.csv", row.names=1)
pt  <- read.csv("results/deg_Preterm_vs_NormalLate.csv", row.names=1)

get_sig <- function(deg) {
  rownames(deg)[!is.na(deg$p_val_adj) & deg$p_val_adj < 0.05 & abs(deg$avg_log2FC) > 0.5]
}

all_genes <- unique(c(get_sig(mis), get_sig(inf), get_sig(pe), get_sig(pt)))
mat <- data.frame(
  Miscarriage = as.integer(all_genes %in% get_sig(mis)),
  Infection   = as.integer(all_genes %in% get_sig(inf)),
  PE          = as.integer(all_genes %in% get_sig(pe)),
  Preterm     = as.integer(all_genes %in% get_sig(pt)),
  row.names = all_genes
)

png("figures/Fig5/补充图Fig5A_UpSet_DEG_overlap.png", w=10, h=6, units="in", res=300, bg="white")
upset(mat, sets=c("Miscarriage","Infection","PE","Preterm"),
      keep.order=TRUE, order.by="freq",
      main.bar.color="#37474F",
      sets.bar.color=c("#D73027","#4575B4","#FC8D59","#66C2A5"),
      text.scale=c(1.5,1.3,1.3,1.1,1.5,1.3))
dev.off()
cat("Fig5S1 UpSet done\n")
