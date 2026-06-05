#!/usr/bin/env Rscript
# Fig4: CellChat — Hofbauer-stromal ligand-receptor interactions
library(CellChat); library(ggplot2); library(patchwork); library(grid)

cellchat <- readRDS("results/cellchat_ucsf_mid.rds")

# ── Panel A: Bubble — incoming signals to HB ──
pdf(NULL)  # suppress device
p_bubble <- netVisual_bubble(cellchat, 
  sources.use=c("FB","fEC","vEC","VCT","SCT","dNK","CD14_M"),
  targets.use="HB", remove.isolate=FALSE, return.data=TRUE)

# ── Panel B: Chord — FB → HB ──
# Generate as grob
png("/tmp/fig4_chord.png", w=6, h=6, units="in", res=300)
netVisual_chord_gene(cellchat, sources.use="FB", targets.use="HB",
  legend.pos.x=8, title.name=NULL)
dev.off()

# ── Panel C: Dotplot — key L-R pairs ──
hb_fb <- subsetCommunication(cellchat, sources.use="FB", targets.use="HB")
hb_fec <- subsetCommunication(cellchat, sources.use="fEC", targets.use="HB")
hb_vec <- subsetCommunication(cellchat, sources.use="vEC", targets.use="HB")

all_pairs <- bind_rows(
  hb_fb %>% arrange(desc(prob)) %>% head(10) %>% mutate(source="FB"),
  hb_fec %>% arrange(desc(prob)) %>% head(8) %>% mutate(source="fEC"),
  hb_vec %>% arrange(desc(prob)) %>% head(6) %>% mutate(source="vEC")
) %>% mutate(pair = paste0(ligand," \u2192 ",receptor))

# Highlight key pairs
all_pairs$highlight <- ifelse(
  all_pairs$pair %in% c("TGFB1 → TGFBR2","THBS1 → CD36",
    "COL1A1 → ITGA9+ITGB1","COL1A2 → ITGA9+ITGB1",
    "FN1 → ITGAV+ITGB1","SPP1 → ITGAV+ITGB1"),
  "bold","plain")

p_dot <- ggplot(all_pairs, aes(x=source, y=reorder(pair, prob), size=prob, fill=pathway_name)) +
  geom_point(shape=21, colour="grey30", stroke=0.3) +
  scale_size_continuous(range=c(3,11), name="Probability") +
  scale_fill_manual(values=c(
    "COLLAGEN"="#4575B4","FN1"="#D73027","TGFb"="#FDAE61",
    "SPP1"="#7B4FA0","LAMININ"="#2D8B57","PTPRM"="grey60",
    "CD45"="grey50","GDF"="grey70","NOTCH"="grey70",
    "BMP"="grey70","PECAM1"="#1B6B93","THBS"="#C62828","JAM"="grey70")) +
  theme_minimal(base_size=11) +
  theme(
    panel.grid.major=element_line(colour="grey92"),
    panel.grid.minor=element_blank(),
    axis.text.y=element_text(size=9.5),
    legend.position="right",
    plot.title=element_text(face="bold", size=13)
  ) +
  labs(x="Source cell type", y="",
    title="Key stromal → Hofbauer ligand-receptor pairs",
    fill="Pathway")

# ── Assemble Fig4 ──
library(png); library(grid)

# Read chord image
chord_img <- readPNG("/tmp/fig4_chord.png")
chord_grob <- rasterGrob(chord_img, interpolate=TRUE)

# Create bubble plot using CellChat's built-in
png("/tmp/fig4_bubble.png", w=10, h=8, units="in", res=300)
netVisual_bubble(cellchat,
  sources.use=c("FB","fEC","vEC","VCT","SCT"),
  targets.use="HB", remove.isolate=FALSE,
  title.name=NULL)
dev.off()
bubble_img <- readPNG("/tmp/fig4_bubble.png")
bubble_grob <- rasterGrob(bubble_img, interpolate=TRUE)

# Layout: top=bubble wide, bottom=chord + dotplot side by side
top <- wrap_elements(bubble_grob) + 
  labs(tag="a") + theme(plot.tag=element_text(face="bold", size=16))

bottom_left <- wrap_elements(chord_grob) +
  labs(tag="b") + theme(plot.tag=element_text(face="bold", size=16))

bottom_right <- wrap_elements(p_dot) +
  labs(tag="c") + theme(plot.tag=element_text(face="bold", size=16))

fig4 <- top / (bottom_left | bottom_right) +
  plot_annotation(
    title="Hofbauer-stromal cell communication in mid-gestation placenta",
    theme=theme(plot.title=element_text(hjust=0.5, face="bold", size=15))
  )

png("figures/Fig4_cellchat.png", w=18, h=14, units="in", res=300)
print(fig4)
dev.off()

cat("Fig4 saved: figures/Fig4_cellchat.png\n")
