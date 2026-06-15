#!/usr/bin/env Rscript
# Figure 1A: 10 datasets, 5 groups, colors from Shiny app.R disease_cols
library(ggplot2); library(dplyr)

OUTDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

datasets <- data.frame(
  name     = c("E-MTAB-6701","E-MTAB-12421","E-MTAB-12795","GSE214607",
               "UCSF Li 2026","GSE290578","GSE298602","GSE333257",
               "GSE298119","GSE329173","GSE173193"),
  gw_start = c(6, 4.5, 4, 6, 11, 29, 37, 32, 37, 37, 37),
  gw_end   = c(12, 12.5, 8.5, 8, 39, 40, 40, 39, 40, 40, 40),
  type     = c("Normal","Normal","Infection","Miscarriage","Normal",
               "PE","PE","Preterm","PE","PE","PE"),
  stringsAsFactors = FALSE
)

datasets$type <- factor(datasets$type,
  levels=c("Normal","Miscarriage","Infection","PE","Preterm"))

# Exact colors from Shiny app.R disease_cols (no Normal 3rd)
type_cols <- c(
  "Normal"      = "#4575B4",
  "PE"          = "#FC8D59",
  "Miscarriage" = "#D73027",
  "Infection"   = "#FDB462",
  "Preterm"     = "#E41A1C"
)

datasets <- datasets[order(datasets$gw_start),]
datasets$y_pos <- nrow(datasets):1

p <- ggplot(datasets) +
  geom_segment(aes(x=gw_start, xend=gw_end, y=y_pos, yend=y_pos, color=type),
               linewidth=4.5, alpha=0.88) +
  geom_text(aes(x=3.5, y=y_pos, label=name), hjust=1, size=3, lineheight=0.9, fontface="bold") +
  scale_color_manual(values=type_cols, name="Condition") +
  scale_x_continuous(breaks=seq(5,40,5), limits=c(0,42),
                     labels=paste0("GW", seq(5,40,5))) +
  scale_y_continuous(limits=c(0.5, nrow(datasets)+0.5)) +
  labs(x="Gestational Age", y="") +
  theme_minimal() + theme(
    text=element_text(size=11), axis.text.y=element_blank(),
    axis.ticks.y=element_blank(), axis.text.x=element_text(color="black", size=9),
    axis.title.x=element_text(size=10), legend.position="top",
    legend.title=element_text(size=9), legend.text=element_text(size=8),
    legend.key.size=unit(0.3,"cm"), panel.grid.major.y=element_blank(),
    panel.grid.minor=element_blank(),
    panel.grid.major.x=element_line(color="grey90", linewidth=0.3))

ggsave(file.path(OUTDIR, "Fig1A_data_overview_10datasets.png"), p, w=8.5, h=4.5, dpi=300, bg="white")
cat("Saved: Fig1A_data_overview_10datasets.png\n")
