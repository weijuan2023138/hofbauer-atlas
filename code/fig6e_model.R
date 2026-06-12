#!/usr/bin/env Rscript
# Fig6e: Dual-track model — biological-style schematic
library(ggplot2); library(grid)

dir.create("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures/Fig6", showWarnings=FALSE)

p <- ggplot() + theme_void()

# === Left: Mid-gestation cell ===
# Cell body (ellipse)
p <- p + annotate("point", x=2.5, y=5.5, size=60, shape=21, fill="#F0F4FF", color="#4575B4", stroke=1.5)
# Nucleus
p <- p + annotate("point", x=2.5, y=5.5, size=25, shape=21, fill="#E0E8F8", color="#4575B4", stroke=0.5)
# CEBPA label inside nucleus
p <- p + annotate("text", x=2.5, y=5.5, label="CEBPA\nHIGH", size=4, fontface="bold", color="#1565C0")
# Repressed target genes below cell
p <- p + annotate("text", x=2.5, y=3.8, label="SPP1  FN1", size=4, fontface="italic", color="#1565C0")
p <- p + annotate("segment", x=2.5, xend=2.5, y=4.2, yend=3.95, linewidth=1, color="#1565C0",
  arrow=arrow(length=unit(0.08,"inch"),type="closed"))

# === Center: Developmental TF switch label ===
p <- p + annotate("segment", x=3.5, xend=6, y=5.5, yend=5.5, linewidth=2, color="grey50",
  arrow=arrow(length=unit(0.15,"inch"),type="closed"))
p <- p + annotate("text", x=4.75, y=6.5, label="Developmental\nTF switch", size=4, fontface="bold", 
  color="grey50", hjust=0.5)

# === Right: Term cell ===
# Cell body
p <- p + annotate("point", x=7, y=5.5, size=60, shape=21, fill="#FFF0F0", color="#C62828", stroke=1.5)
# Nucleus (smaller = CEBPA lower)
p <- p + annotate("point", x=6.8, y=6.0, size=18, shape=21, fill="#F8E0E0", color="#C62828", stroke=0.5)
p <- p + annotate("text", x=6.8, y=6.0, label="CEBPA\nlow", size=3.5, fontface="bold", color="#C62828")
# STAT/NFKB in cytoplasm
p <- p + annotate("point", x=7.2, y=5.0, size=22, shape=21, fill="#FFE0E0", color="#C62828", stroke=0.5)
p <- p + annotate("text", x=7.2, y=5.0, label="STAT3\nNFKB1", size=3.5, fontface="bold", color="#C62828")

# Arrows: CEBPA stops inhibiting, STAT/NFKB activate
p <- p + annotate("segment", x=6.8, xend=6.8, y=4.5, yend=3.6, linewidth=1.5, color="#1565C0",
  arrow=arrow(length=unit(0.1,"inch"),type="closed"))
p <- p + annotate("text", x=5.5, y=3.3, label="De-repression", size=3.5, fontface="bold.italic", color="#1565C0")

p <- p + annotate("segment", x=7.2, xend=7.2, y=4.5, yend=2.6, linewidth=1.5, color="#C62828",
  arrow=arrow(length=unit(0.1,"inch"),type="closed"))
p <- p + annotate("text", x=8.3, y=2.3, label="Direct activation\n(motif-driven)", size=3.5, fontface="bold.italic", color="#C62828")

# Output labels
p <- p + annotate("text", x=5.8, y=2.8, label="SPP1↑  FN1↑\n(ECM remodeling)", size=4, fontface="bold", color="#1565C0", hjust=1)
p <- p + annotate("text", x=8.5, y=1.8, label="CD44↑  PTPRM↑\nCD47↑  IL1B↑\n(Immune adhesion)", size=4, fontface="bold", color="#C62828", hjust=0)

# Chromatin below (ATAC tracks - simple representation)
p <- p + annotate("rect", xmin=4.5, xmax=9.5, ymin=0.2, ymax=1.2, fill="#FFF8F0", color="grey70", size=0.5)
p <- p + annotate("text", x=7, y=0.9, label="Open chromatin at STAT/NF-κB motifs", size=3.5, fontface="italic")
p <- p + annotate("text", x=7, y=0.5, label="Closed chromatin at CEBPA targets (SPP1)", size=3.5, fontface="italic")

# Title
p <- p + annotate("text", x=4.75, y=7.8, 
  label="Dual-track TF model programming Hofbauer cell communication identity",
  size=6, fontface="bold")

# Sub-labels
p <- p + annotate("text", x=2.5, y=7.2, label="Mid-gestation", size=4.5, fontface="bold", color="#4575B4")
p <- p + annotate("text", x=7, y=7.2, label="Term", size=4.5, fontface="bold", color="#C62828")

p <- p + xlim(1,10) + ylim(0,8.5)

ggsave("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures/Fig6/Fig6e_model.png", 
       p, w=12, h=9, dpi=300, bg="white")
message("Fig6e done")
