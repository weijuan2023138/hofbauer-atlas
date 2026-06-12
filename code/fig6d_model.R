#!/usr/bin/env Rscript
# Fig6d: Dual-track TF → communication gene regulation model — polished
setwd("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration")
library(ggplot2); library(grid)

# === Track boxes (background) ===
boxes <- data.frame(
  xmin = c(0.3, 0.3), xmax = c(9.7, 9.7),
  ymin = c(3.6, 0.8), ymax = c(5.6, 2.8),
  track = c("De-repression Track", "Activation Track"),
  col  = c("#D73027", "#4575B4"),
  fill = c("#FFF0F0", "#F0F4FF")
)

# === Main nodes ===
nodes <- data.frame(
  x    = c(1.5, 4.0, 6.5, 8.5,  1.5, 4.0, 6.5, 8.5),
  y    = c(4.6, 4.6, 4.6, 4.6,  2.2, 2.2, 2.2, 2.2),
  label= c("CEBPA ↓", "Chromatin\nclosing", "SPP1 / FN1", "ECM secretion\n& remodeling",
           "NFKB1/STAT3 ↑", "Chromatin\nopening", "PTPRM/CD44\nCD47", "Adhesion\n& homing"),
  fill = c("#4575B4","#FC8D59","#D73027","#D73027",
           "#D73027","#91BFDB","#4575B4","#4575B4"),
  track= rep(c("De-repression","Activation"), each=4)
)

# === Arrows ===
arrows <- data.frame(
  x    = c(2.2, 5.0, 7.5,  2.2, 5.0, 7.5),
  xend = c(3.4, 5.8, 8.1,  3.4, 5.8, 8.1),
  y    = c(4.6, 4.6, 4.6,  2.2, 2.2, 2.2),
  yend = c(4.6, 4.6, 4.6,  2.2, 2.2, 2.2)
)

# === P-value callouts ===
callouts <- data.frame(
  x = c(5.0, 5.0), y = c(5.3, 2.9),
  label = c("Motif: CEBPA ns   |   SPP1 fold=1.16  FN1 fold=1.49",
            "Motif: STAT3 ***   |   CD44 fold=1.98  CD47 fold=1.21"),
  color = c("#D73027", "#4575B4")
)

# === Mechanism labels ===
mech <- data.frame(
  x = c(2.8, 2.8), y = c(5.0, 2.6),
  label = c("Indirect de-repression\n(no motif enrichment)", "Direct motif-driven\nactivation"),
  color = c("#FC8D59", "#4575B4")
)

p <- ggplot() +
  # Track backgrounds
  geom_rect(data=boxes, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill=boxes$fill, color=boxes$col, linewidth=0.8, linetype="dashed", alpha=0.5) +
  # Track titles
  geom_text(data=boxes, aes(x=5, y=ymax-0.15, label=track), size=4.5, fontface="bold", color=boxes$col) +
  # Arrows
  geom_segment(data=arrows, aes(x=x, xend=xend, y=y, yend=yend),
               arrow=arrow(length=unit(0.15,"inch"), type="closed"),
               color="grey40", linewidth=1.2) +
  # Mechanism text
  geom_text(data=mech, aes(x=x, y=y, label=label), size=3.2, color=mech$color,
            fontface="italic", hjust=0.5, lineheight=0.9) +
  # Main nodes
  geom_label(data=nodes, aes(x=x, y=y, label=label, fill=fill),
             color="white", size=4, fontface="bold", linewidth=0.3,
             label.padding=unit(0.5,"lines"), label.r=unit(0.3,"lines")) +
  scale_fill_identity() +
  # Stats callouts
  geom_text(data=callouts, aes(x=x, y=y, label=label), size=3.5, color=callouts$color,
            fontface="italic", hjust=0.5) +
  # Central title
  annotate("text", x=5, y=6.4, label="Dual-track TF Regulation of Hofbauer Cell Communication Identity",
           size=6, fontface="bold") +
  # "Hofbauer Cell" badge
  annotate("rect", xmin=3.8, xmax=6.2, ymin=6.8, ymax=7.3, fill="#37474F", color=NA) +
  annotate("text", x=5, y=7.05, label="Hofbauer Cell", size=5, fontface="bold", color="white") +
  # Divider lines connecting tracks
  annotate("segment", x=5, xend=5, y=5.6, yend=6.4, color="grey70", linewidth=0.4, linetype="dotted") +
  xlim(0, 10) + ylim(0.2, 7.6) +
  theme_void()

ggsave("figures/Fig6/Fig6d_model.png", p, w=14, h=6.5, dpi=300, bg="white")
message("Fig6d done")
