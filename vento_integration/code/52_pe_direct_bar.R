#!/usr/bin/env Rscript
# PE early vs late direct GSEA — bidirectional bar chart
library(ggplot2); library(dplyr); library(stringr)
set.seed(42)

d <- read.csv("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/PE_EarlyDirect_vs_Late_GSEA.csv")
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

# Top 20 pathways by abs(NES), padj<0.1
sig <- d %>% filter(p.adjust<0.1) %>% slice_max(abs(NES), n=20) %>%
  mutate(Description=str_wrap(Description, width=50))

# Order by NES
sig <- sig %>% arrange(-NES)
sig$Description <- factor(sig$Description, levels=rev(sig$Description))

# Color: Early red, Late blue
sig$direction <- ifelse(sig$NES>0, "Early PE", "Late PE")

p <- ggplot(sig, aes(x=Description, y=NES, fill=direction)) +
  geom_col(width=0.6, alpha=0.9) +
  scale_fill_manual(values=c("Early PE"="#C62828","Late PE"="#1565C0")) +
  geom_hline(yintercept=0, color="black", linewidth=0.6) +
  coord_flip() +
  labs(title="GO:BP — Early PE vs Late PE", x="", y="NES (positive = Early PE biased)") +
  theme_minimal(12) +
  theme(axis.text.y=element_text(size=9.5,color="black",face="bold"),
    axis.text.x=element_text(size=9,color="black"),
    panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
    panel.grid.major.x=element_line(color="grey92",linewidth=0.3),
    plot.title=element_text(hjust=0.5,size=13,face="bold"),
    legend.position="top",legend.title=element_blank(),legend.text=element_text(size=11))
ggsave(file.path(FIGDIR,"PE_early_vs_late_direct_bar.png"), p, w=8, h=6.5, dpi=300, bg="white")
cat("Saved\n")
