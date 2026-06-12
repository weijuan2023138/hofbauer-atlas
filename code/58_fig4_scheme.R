#!/usr/bin/env Rscript
# 58_fig4_scheme.R
# CellChat analysis workflow schematic for Figure 4

library(ggplot2)
library(grid)

# Create a schematic diagram
scheme_data <- data.frame(
  x = c(1, 3, 5, 7, 9),
  y = c(1, 1, 1, 1, 1),
  label = c(
    "Single-cell\nRNA-seq data",
    "CellChat\nAnalysis",
    "Ligand-Receptor\nPairs",
    "Communication\nProbability",
    "Significant\nInteractions"
  ),
  stringsAsFactors = FALSE
)

# Arrow data
arrows <- data.frame(
  x = c(1.8, 3.8, 5.8, 7.8),
  xend = c(2.2, 4.2, 6.2, 8.2),
  y = c(1, 1, 1, 1),
  yend = c(1, 1, 1, 1)
)

# Color scheme
box_colors <- c("#4575B4", "#FDAE61", "#D73027", "#7570B3", "#1B9E77")

p <- ggplot() +
  # Draw boxes
  geom_rect(data=scheme_data,
    aes(xmin=x-0.8, xmax=x+0.8, ymin=y-0.4, ymax=y+0.4),
    fill=box_colors, color="black", linewidth=0.8, alpha=0.9) +
  # Add text labels
  geom_text(data=scheme_data,
    aes(x=x, y=y, label=label),
    size=3.5, fontface="bold", color="white") +
  # Draw arrows
  geom_segment(data=arrows,
    aes(x=x, y=y, xend=xend, yend=yend),
    arrow=arrow(length=unit(0.3, "cm"), type="closed"),
    linewidth=1.2, color="black") +
  # Add step numbers
  annotate("text", x=scheme_data$x, y=scheme_data$y+0.6,
    label=paste0("(", 1:5, ")"),
    size=3.5, fontface="bold", color="black") +
  # Add explanation text below
  annotate("text", x=5, y=-0.3,
    label=paste0(
      "(1) Single-cell RNA-seq data from placental tissues  ",
      "(2) CellChat identifies cell-cell communication based on L-R pairs  ",
      "(3) Ligand-receptor pairs: specific molecular interactions  ",
      "(4) Communication probability: strength of interaction  ",
      "(5) Statistical testing (p < 0.01) identifies significant interactions"
    ),
    size=3, fontface="plain", hjust=0.5, lineheight=1.3) +
  # Theme
  theme_void() +
  theme(
    plot.margin=margin(t=20, r=20, b=40, l=20),
    plot.title=element_text(face="bold", size=14, hjust=0.5, margin=margin(b=15))
  ) +
  labs(title="CellChat Analysis Workflow") +
  coord_cartesian(xlim=c(0, 10), ylim=c(-0.8, 2.5))

ggsave("figures/Fig4/Fig4_cellchat_scheme.png",
  p, width=10, height=4, dpi=300, bg="white")

cat("Done: Fig4_cellchat_scheme.png\n")
