#!/usr/bin/env Rscript
# Fig6c: TF Motif Enrichment — grouped bar with comparison brackets
setwd("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration")
library(ggplot2); library(dplyr); library(tidyr)

motif <- read.csv("results/atac_motif_enrichment.csv")
tf_order <- c("STAT3","STAT1","RELB","CEBPA","NFKB1")
motif_5 <- motif[motif$name %in% tf_order,]
motif_5 <- motif_5[match(tf_order, motif_5$name),]

motif_long <- motif_5 %>% select(name, term_pct, mid_pct) %>%
  rename(Term=term_pct, Mid=mid_pct) %>%
  pivot_longer(-name, names_to="Condition", values_to="pct")
motif_long$name <- factor(motif_long$name, levels=tf_order)
motif_long$Condition <- factor(motif_long$Condition, levels=c("Mid","Term"))

# Brackets: between Mid and Term bars per TF
dodge_w <- 0.7
bracket_h <- 5
motif_5$xpos <- 1:nrow(motif_5)
motif_5$ymax <- pmax(motif_5$term_pct, motif_5$mid_pct)
motif_5$p_label <- ifelse(motif_5$pval<0.001, "***",
                    ifelse(motif_5$pval<0.01, "**",
                    ifelse(motif_5$pval<0.05, "*", "ns")))

p <- ggplot(motif_long, aes(x=name, y=pct, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(width=dodge_w), width=0.6,
           color="black", linewidth=0.3) +
  # Percentage labels above bars
  geom_text(data=motif_long, aes(label=sprintf("%.0f%%", pct), y=pct+1.5),
            position=position_dodge(width=dodge_w), size=2.8, fontface="bold", color="black") +
  # Bracket horizontal line between bars
  geom_segment(data=motif_5, aes(x=xpos-0.15, xend=xpos+0.15,
              y=ymax+bracket_h, yend=ymax+bracket_h),
              inherit.aes=FALSE, linewidth=0.5, color="black") +
  # Significance label above bracket — ns smaller
  geom_text(data=subset(motif_5, p_label!="ns"),
            aes(x=xpos, y=ymax+bracket_h+3, label=p_label),
            inherit.aes=FALSE, size=5, fontface="bold") +
  geom_text(data=subset(motif_5, p_label=="ns"),
            aes(x=xpos, y=ymax+bracket_h+3, label=p_label),
            inherit.aes=FALSE, size=3.5, fontface="plain") +
  scale_fill_manual(values=c("Mid"="#4575B4","Term"="#D73027")) +
  labs(x="", y="Peaks with Motif (%)",
       title="TF Binding Motifs in Differential ATAC Peaks",
       subtitle="Mid vs Term Hofbauer cells") +
  ylim(0, max(motif_5$ymax)+bracket_h+8) +
  theme_classic(base_size=11) +
  theme(legend.position="top", legend.title=element_text(face="bold",size=10),
        legend.text=element_text(size=9),
        plot.title=element_text(face="bold",size=14,hjust=0.5),
        plot.subtitle=element_text(size=10,hjust=0.5,color="grey40"),
        axis.text.x=element_text(size=12,face="bold",color="black"),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=11))

ggsave("figures/Fig6/Fig6c_motif_enrichment.png", p, w=7, h=5, dpi=300, bg="white")
message("Done: figures/Fig6/Fig6c_motif_enrichment.png")
