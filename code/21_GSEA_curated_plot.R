#!/usr/bin/env Rscript
# Refined GSEA dot plot — curated pathways, no redundancy
library(ggplot2); library(dplyr); library(stringr)

d <- read.csv('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/dev_trimester_GSEA.csv')
FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'

# ---- Curate: keep one representative pathway per functional group ----
# Remove redundant terms (keep the most specific / highest |NES|)
curated <- c(
  # Early-enriched: mitochondrial & energy
  "oxidative phosphorylation",
  "mitochondrial translation",
  "ATP synthesis coupled electron transport",
  "ribosome biogenesis",
  "rRNA processing",
  # Early-enriched: cell cycle
  "chromosome segregation",
  "mitotic spindle assembly checkpoint signaling",
  # Late-enriched: immune activation
  "cytoplasmic translation",
  "inflammatory response",
  "positive regulation of cytokine production",
  "positive regulation of T cell activation",
  "leukocyte migration",
  "leukocyte mediated immunity",
  "immune response-activating signaling pathway",
  "humoral immune response",
  "antigen processing and presentation of peptide antigen via MHC class II",
  # Mid-enriched or signaling
  "regulation of small GTPase mediated signal transduction",
  "cilium organization",
  "cell-substrate adhesion",
  "enzyme-linked receptor protein signaling pathway",
  "cell morphogenesis involved in neuron differentiation"
)

plot_data <- d %>% 
  filter(Description %in% curated) %>%
  mutate(
    Trimester = factor(Trimester, levels=c("Early","Mid","Late")),
    Description = str_wrap(Description, width=45)
  )

# Sort by functional group: Early-enriched → Mid → Late-enriched → others
# Assign functional group order
func_order <- c(
  "oxidative phosphorylation"=1,"ATP synthesis coupled electron transport"=1,
  "mitochondrial translation"=1,"ribosome biogenesis"=1,"rRNA processing"=1,
  "chromosome segregation"=1,"mitotic spindle assembly checkpoint signaling"=1,
  "regulation of small GTPase mediated signal transduction"=2,"cell-substrate adhesion"=2,
  "enzyme-linked receptor protein signaling pathway"=2,
  "cilium organization"=2,"cell morphogenesis involved in neuron differentiation"=2,
  "cytoplasmic translation"=3,"positive regulation of T cell activation"=3,
  "leukocyte mediated immunity"=3,"inflammatory response"=3,
  "positive regulation of cytokine production"=3,"leukocyte migration"=3,
  "humoral immune response"=3,"immune response-activating signaling pathway"=3,
  "antigen processing and presentation of peptide antigen via MHC class II"=3
)
plot_data$group <- func_order[as.character(plot_data$Description)]
desc_order <- plot_data %>% group_by(Description) %>% 
  summarise(g=first(group), .groups="drop") %>% arrange(g, Description) %>% pull(Description)
plot_data$Description <- factor(plot_data$Description, levels=rev(desc_order))

cat(sprintf("%d curated pathways, %d data points\n", 
  length(unique(plot_data$Description)), nrow(plot_data)))

p <- ggplot(plot_data, aes(x=Trimester, y=Description)) +
  geom_hline(yintercept=seq_along(levels(plot_data$Description)), color="grey92", linewidth=0.3) +
  geom_point(aes(size=Count, color=NES), stroke=0) +
  scale_color_gradient2(low="#313695", mid="#FFFFBF", high="#A50026",
                        midpoint=0, name="NES",
                        limits=c(-3.2, 3.2)) +
  scale_size_continuous(range=c(3, 9), name="Core\nGenes") +
  theme_bw() +
  labs(title="GO Biological Process — GSEA by Trimester",
       subtitle="Early (GW4.5–10) / Mid (GW11–24) / Late (GW32–38)",
       x="", y="") +
  theme(
    text=element_text(face="bold"),
    axis.text.x=element_text(size=12, color="black"),
    axis.text.y=element_text(size=9, color="black"),
    legend.title=element_text(size=10),
    legend.text=element_text(size=9),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.border=element_rect(color="black", linewidth=0.6),
    plot.title=element_text(hjust=0.5, size=13),
    plot.subtitle=element_text(hjust=0.5, face="plain", size=10)
  )

ggsave(file.path(FIGDIR,"Fig_dev_GSEA_dotplot_curated.png"), p, w=9, h=7.5, dpi=300, bg="white")
cat("\nSaved Fig_dev_GSEA_dotplot_curated.png\n")
