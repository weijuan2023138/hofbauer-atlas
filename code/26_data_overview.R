#!/usr/bin/env Rscript
# Figure 1A: Data overview schematic — gestational timeline + dataset bar chart
library(ggplot2); library(dplyr)

# ---- Dataset info ----
datasets <- data.frame(
  name = c("Arutyunyan\n2023","Vento-Tormo\n2018","UCSF Li\n2026","GSE290578\n(Normal)","GSE290578\n(PE)","GSE214607\n(Miscarriage)","Hoo\n2024","GSE183338\n(PE)","my_preterm\ncohort"),
  gw_start = c(4.5, 6, 11, 32, 32, 6, 6, 27, 28),
  gw_end   = c(10, 12, 24, 38, 38, 9, 12, 37, 37),
  n_cells  = c(3172, 1200, 461, 1148, 4998, 1800, 2100, 800, 1100),
  type     = c("Normal","Normal","Normal","Normal","PE","Miscarriage","Infection","PE","Preterm"),
  tech     = c("snRNA-seq","scRNA-seq","snRNA-seq","scRNA-seq","scRNA-seq","scRNA-seq","scRNA-seq","snRNA-seq","scRNA-seq"),
  stringsAsFactors = FALSE
)

datasets$type <- factor(datasets$type, levels=c("Normal","PE","Miscarriage","Infection","Preterm"))
type_cols <- c("Normal"="#4575B4","PE"="#D73027","Miscarriage"="#FDAE61","Infection"="#7B4FA0","Preterm"="#2D8B57")

# Sort by gestational age
datasets <- datasets[order(datasets$gw_start),]
datasets$y_pos <- nrow(datasets):1

p <- ggplot(datasets) +
  # Gestational week axis segments
  geom_segment(aes(x=gw_start, xend=gw_end, y=y_pos, yend=y_pos, color=type),
               linewidth=4, alpha=0.85) +
  # Cell count labels
  geom_text(aes(x=(gw_start+gw_end)/2, y=y_pos, label=paste0("n=", n_cells)),
            color="white", size=3, fontface="bold") +
  # Dataset name on the left
  geom_text(aes(x=3.5, y=y_pos, label=name), hjust=1, size=3.2, lineheight=0.9) +
  scale_color_manual(values=type_cols, name="Condition") +
  scale_x_continuous(breaks=seq(5,40,5), limits=c(0,42),
                     labels=paste0("GW", seq(5,40,5))) +
  scale_y_continuous(limits=c(0.5, nrow(datasets)+0.5)) +
  labs(x="Gestational Age", y="") +
  theme_minimal() +
  theme(
    text=element_text(size=11),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_text(color="black", size=9),
    axis.title.x=element_text(size=10),
    legend.position="top",
    legend.title=element_text(size=9),
    legend.text=element_text(size=8),
    legend.key.size=unit(0.3,"cm"),
    panel.grid.major.y=element_blank(),
    panel.grid.minor=element_blank(),
    panel.grid.major.x=element_line(color="grey90", linewidth=0.3)
  )

FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'
ggsave(file.path(FIGDIR,"Fig1A_data_overview.png"), p, w=8, h=4.5, dpi=300, bg="white")
cat("Saved Fig1A_data_overview.png\n")
