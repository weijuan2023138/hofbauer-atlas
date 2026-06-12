#!/usr/bin/env Rscript
# Test multiple high-impact journal color schemes for GSEA dot plot
library(ggplot2); library(dplyr); library(stringr); library(patchwork)

d <- read.csv('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/dev_trimester_GSEA.csv')
FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'

curated <- c(
  "oxidative phosphorylation","mitochondrial translation",
  "ATP synthesis coupled electron transport","ribosome biogenesis","rRNA processing",
  "chromosome segregation","mitotic spindle assembly checkpoint signaling",
  "cytoplasmic translation","inflammatory response",
  "positive regulation of cytokine production","positive regulation of T cell activation",
  "leukocyte migration","leukocyte mediated immunity",
  "immune response-activating signaling pathway","humoral immune response",
  "antigen processing and presentation of peptide antigen via MHC class II",
  "regulation of small GTPase mediated signal transduction","cilium organization",
  "cell-substrate adhesion","enzyme-linked receptor protein signaling pathway",
  "cell morphogenesis involved in neuron differentiation"
)

plot_data <- d %>% 
  filter(Description %in% curated) %>%
  mutate(
    Trimester = factor(Trimester, levels=c("Early","Mid","Late")),
    Description = str_wrap(Description, width=45)
  )

desc_order <- plot_data %>% group_by(Description) %>% 
  summarise(m=sum(NES, na.rm=TRUE)) %>% arrange(desc(m)) %>% pull(Description)
plot_data$Description <- factor(plot_data$Description, levels=rev(desc_order))

make_plot <- function(low, mid, high, name, subtitle) {
  ggplot(plot_data, aes(x=Trimester, y=Description)) +
    geom_hline(yintercept=seq_along(levels(plot_data$Description)), color="grey92", linewidth=0.3) +
    geom_point(aes(size=Count, color=NES), stroke=0) +
    scale_color_gradient2(low=low, mid=mid, high=high, midpoint=0, name="NES", limits=c(-3.2,3.2)) +
    scale_size_continuous(range=c(3,9), name="Core\nGenes") +
    theme_bw() +
    labs(title=name, subtitle=subtitle, x="", y="") +
    theme(
      text=element_text(face="bold"),
      axis.text.x=element_text(size=11, color="black"),
      axis.text.y=element_text(size=8, color="black"),
      legend.title=element_text(size=9), legend.text=element_text(size=8),
      panel.grid=element_blank(),
      panel.border=element_rect(color="black", linewidth=0.5),
      plot.title=element_text(hjust=0.5, size=12),
      plot.subtitle=element_text(hjust=0.5, face="plain", size=9)
    )
}

# ---- 4 schemes ----

# 1. Nature Comms — RdYlBu (blue-yellow-red, Nature's go-to diverging)
p1 <- make_plot("#313695","#FFFFBF","#A50026","Nature Communications","RdYlBu — blue / yellow / red")
# 2. Cell — deep blue-white-deep red (classic RdBu)
p2 <- make_plot("#2166AC","#F7F7F7","#B2182B","Cell / Science","RdBu — deep blue / white / deep red")
# 3. PNAS — cool-warm (teal-white-coral, modern)
p3 <- make_plot("#018571","#F5F5F5","#A6611A","PNAS / Modern","BrBG — teal / white / brown-orange")
# 4. Cell Systems — purple-white-orange (popular in 2024-2025)
p4 <- make_plot("#5E3C99","#F7F7F7","#E66101","Cell Systems / Trends","PuOr — purple / white / orange")

combined <- (p1 + p2) / (p3 + p4) + 
  plot_annotation(title="GSEA Color Scheme Comparison", theme=theme(plot.title=element_text(hjust=0.5,face="bold",size=14)))

ggsave(file.path(FIGDIR,"Fig_dev_GSEA_color_comparison.png"), combined, w=16, h=14, dpi=300, bg="white")
cat("Saved color comparison\n")
