#!/usr/bin/env Rscript
# Fig1A: Study design + model schematic
library(ggplot2); library(dplyr); library(grid); library(ggrepel)

# ── Colour palette ──
st_cols <- c("Pro-inflammatory"="#C62828","MHCII+ Antigen-presenting"="#BF4E1A",
  "Homeostatic"="#1B6B93","PRKN+ Autophagy"="#7B4FA0",
  "Vascular remodeling"="#2D8B57","MKI67+ Proliferating"="#37474F")
tri_cols <- c("Early"="#4575B4","Mid"="#FDAE61","Late"="#D73027")

# ── Layout helpers ──
blank <- theme_void() + theme(plot.margin=margin(0,0,0,0))

# ── Panel 1: Datasets ──
datasets <- data.frame(
  name = c("Arutyunyan\n2024","UCSF\nLi 2026","GSE290578\n(PTB)","GSE290578\n(Normal)",
           "Hoo 2024","GSE214607","My Cohort\n(PTB)","My Cohort\n(Normal)"),
  x = rep(1:4, 2), y = rep(c(2,1), each=4),
  type = c(rep("Public",6), rep("In-house",2)),
  stage = c("Early","Mid","Late","Late","Mid","Early","Mid/Late","Mid/Late")
)
p1 <- ggplot(datasets, aes(x,y)) +
  geom_tile(aes(fill=stage), width=0.85, height=0.7, colour="grey30", size=0.4, radius=unit(4,"pt")) +
  geom_text(aes(label=name), size=2.3, lineheight=0.85, colour="white", fontface="bold") +
  scale_fill_manual(values=c("Early"="#4575B4","Mid"="#FDAE61","Late"="#D73027","Mid/Late"="grey60")) +
  annotate("text", x=2.5, y=2.8, label="8 Datasets", size=4, fontface="bold") +
  annotate("text", x=2.5, y=-0.2, label="GW4.5 – GW38", size=3, colour="grey40") +
  blank + theme(legend.position="none") + xlim(0.3,4.7) + ylim(-0.5,3.2)

# ── Panel 2: Integration + Harmony ──
p2 <- ggplot() + 
  annotate("text", x=1, y=1, label="Integration\nHarmony", size=3.5, fontface="bold", colour="#1B6B93") +
  annotate("segment", x=0.5, xend=1.5, y=1.15, yend=1.15, arrow=arrow(length=unit(3,"pt")), colour="grey50") +
  blank + xlim(0,2) + ylim(0.7,1.5)

# ── Panel 3: 6 Subtype UMAP (approximate circles) ──
circles <- data.frame(
  subtype = names(st_cols),
  x = c(1,3,5,2,4,2.5), y = c(2.5,2.5,2.5,1,1,3.8),
  size = c(1.2,0.9,1.1,0.7,0.8,0.5),
  label = c("Pro-inflammatory","MHCII+ AP","Homeostatic","PRKN+\nAutophagy","Vascular\nRemodeling","MKI67+\nProlif.")
)
p3 <- ggplot(circles, aes(x,y)) +
  geom_point(aes(size=size, fill=subtype), shape=21, colour="grey20", stroke=0.3) +
  geom_text_repel(aes(label=label, colour=subtype), size=2.5, lineheight=0.85, fontface="bold",
    box.padding=0.5, point.padding=0.3, segment.size=0.3, max.overlaps=20) +
  scale_fill_manual(values=st_cols) + scale_colour_manual(values=st_cols) +
  scale_size_continuous(range=c(8,20)) +
  annotate("text", x=3, y=4.5, label="6 Hofbauer Subtypes\n17,896 cells", size=3.5, fontface="bold") +
  blank + theme(legend.position="none") + xlim(-1,7) + ylim(-0.5,5.2)

# ── Panel 4: Developmental axis ──
times <- data.frame(
  stage = c("Early","Mid","Late"),
  x = c(1,3.5,6),
  label = c("GW4.5–10\nFetal\nProgenitor-like","GW11–24\nMid-gestation\nMetabolic active","GW32–38\nTerm\nImmune effector"),
  func = c("Proliferation\nTissue remodeling","Glycolysis\nLipid metabolism","Complement\nCytokine/chemokine\nAntigen presentation")
)
p4 <- ggplot(times, aes(x, y=0)) +
  geom_segment(aes(x=0.3, xend=7, y=0, yend=0), colour="grey70", size=1.5) +
  geom_point(aes(fill=stage), shape=21, size=8, colour="grey30", stroke=0.6) +
  geom_text(aes(y=-1.2, label=label), size=2.5, lineheight=0.9, fontface="bold") +
  geom_text(aes(y=1.2, label=func), size=2.3, lineheight=0.9, colour="grey30") +
  scale_fill_manual(values=tri_cols) +
  annotate("text", x=3.5, y=2.2, label="Developmental Trajectory", size=3.5, fontface="bold") +
  blank + theme(legend.position="none") + ylim(-2.5,2.8) + xlim(0,7.5)

# ── Panel 5: Key TFs ──
tfs <- data.frame(
  name = c("CEBPA","MAFB","ID2","SOX4","NFKB1","RELB","NR4A3","IRF5"),
  stage = c("Early","Early","Early","Early","Late","Late","Late","Late"),
  x = c(1,2,3,4,1,2,3,4), y = c(rep(2,4), rep(1,4))
)
p5 <- ggplot(tfs, aes(x,y)) +
  geom_tile(aes(fill=stage), width=0.8, height=0.6, colour="white", size=1, radius=unit(3,"pt")) +
  geom_text(aes(label=name), size=2.8, colour="white", fontface="bold") +
  scale_fill_manual(values=c("Early"="#4575B4","Late"="#D73027")) +
  annotate("text", x=2.5, y=2.8, label="Key Transcription Factors", size=3.5, fontface="bold") +
  blank + theme(legend.position="none") + xlim(0.3,4.7) + ylim(0.3,3.2)

# ── Assemble ──
library(patchwork)

pdf("figures/Fig1/Fig1A_model.pdf", w=14, h=6)
print((p1 | p3) / (p4 | p5) + 
  plot_annotation(title="Fetal Hofbauer Cell Atlas — Study Design & Model",
    theme=theme(plot.title=element_text(hjust=0.5, face="bold", size=14))))
dev.off()

png("figures/Fig1/Fig1A_model.png", w=14, h=6, units="in", res=300)
print((p1 | p3) / (p4 | p5) + 
  plot_annotation(title="Fetal Hofbauer Cell Atlas — Study Design & Model",
    theme=theme(plot.title=element_text(hjust=0.5, face="bold", size=14))))
dev.off()
cat("Fig1A saved\n")
