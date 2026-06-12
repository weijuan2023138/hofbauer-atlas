#!/usr/bin/env Rscript
# 56_fig4_bubble_functional_category.R
# Add biologically meaningful functional category strips to CellChat bubble plots
# Categories: Immune / Inflammation / ECM-Adhesion / Angiogenesis-Growth / Development / Phagocytosis

library(CellChat)
library(ggplot2)
library(dplyr)

cellchat <- readRDS("results/cellchat_ucsf_mid.rds")
ct_levels <- c("CD14_M","dNK","FB","fEC","SCT","VCT","vEC")

# ── LR-pair -> Functional Category mapping ──
# Each entry: "LIGAND - RECEPTOR" -> category
lr_category <- list(
  # === Immune Regulation (免疫调控) ===
  "CD86 - CD28"             = "Immune Regulation",
  "HLA-DRA - CD4"           = "Immune Regulation",
  "HLA-DRB1 - CD4"          = "Immune Regulation",
  "PTPRC - MRC1"            = "Immune Regulation",
  "CD46 - JAG1"             = "Immune Regulation",
  "MIF - (CD74+CD44)"       = "Immune Regulation",
  "HLA-E - CD94:NKG2A"      = "Immune Regulation",
  "HLA-E - CD94:NKG2C"      = "Immune Regulation",
  "HLA-E - CD94:NKG2E"      = "Immune Regulation",
  "HLA-E - KLRC1"           = "Immune Regulation",
  "HLA-E - KLRC2"           = "Immune Regulation",
  "HLA-E - KLRK1"           = "Immune Regulation",
  "VCAM1 - (ITGA9+ITGB1)"   = "Immune Regulation",
  "APP - CD74"              = "Immune Regulation",
  "ADGRE5 - CD55"           = "Immune Regulation",

  # === Inflammation / Chemotaxis (炎症/趋化) ===
  "RARRES2 - CMKLR1"        = "Inflammation",
  "GDF15 - TGFBR2"          = "Inflammation",
  "NAMPT - INSR"            = "Inflammation",

  # === ECM / Adhesion (ECM/黏附) ===
  "COL1A1 - (ITGA9+ITGB1)"  = "ECM & Adhesion",
  "COL1A2 - (ITGA1+ITGB1)"  = "ECM & Adhesion",
  "COL1A2 - (ITGA9+ITGB1)"  = "ECM & Adhesion",
  "COL1A2 - (ITGAV+ITGB8)"  = "ECM & Adhesion",
  "COL1A2 - CD44"           = "ECM & Adhesion",
  "COL1A2 - SDC1"           = "ECM & Adhesion",
  "COL4A1 - (ITGA1+ITGB1)"  = "ECM & Adhesion",
  "COL4A1 - (ITGA9+ITGB1)"  = "ECM & Adhesion",
  "COL4A1 - (ITGAV+ITGB8)"  = "ECM & Adhesion",
  "COL4A1 - CD44"           = "ECM & Adhesion",
  "COL4A1 - SDC1"           = "ECM & Adhesion",
  "COL4A2 - (ITGA9+ITGB1)"  = "ECM & Adhesion",
  "COL4A5 - (ITGA9+ITGB1)"  = "ECM & Adhesion",
  "COL6A1 - (ITGA9+ITGB1)"  = "ECM & Adhesion",
  "COL6A2 - (ITGA9+ITGB1)"  = "ECM & Adhesion",
  "COL6A3 - (ITGA9+ITGB1)"  = "ECM & Adhesion",
  "FN1 - (ITGA4+ITGB1)"     = "ECM & Adhesion",
  "FN1 - (ITGA5+ITGB1)"     = "ECM & Adhesion",
  "FN1 - (ITGAV+ITGB1)"     = "ECM & Adhesion",
  "FN1 - (ITGAV+ITGB8)"     = "ECM & Adhesion",
  "FN1 - CD44"              = "ECM & Adhesion",
  "FN1 - SDC1"              = "ECM & Adhesion",
  "LAMA2 - (ITGA9+ITGB1)"   = "ECM & Adhesion",
  "LAMA3 - (ITGA9+ITGB1)"   = "ECM & Adhesion",
  "LAMA4 - (ITGA9+ITGB1)"   = "ECM & Adhesion",
  "LAMB1 - (ITGA9+ITGB1)"   = "ECM & Adhesion",
  "LAMC1 - (ITGA9+ITGB1)"   = "ECM & Adhesion",
  "TNXB - (ITGA9+ITGB1)"    = "ECM & Adhesion",
  "POSTN - (ITGAV+ITGB5)"   = "ECM & Adhesion",
  "SPP1 - (ITGA4+ITGB1)"    = "ECM & Adhesion",
  "SPP1 - (ITGA5+ITGB1)"    = "ECM & Adhesion",
  "SPP1 - (ITGA9+ITGB1)"    = "ECM & Adhesion",
  "SPP1 - (ITGAV+ITGB1)"    = "ECM & Adhesion",
  "SPP1 - (ITGAV+ITGB5)"    = "ECM & Adhesion",
  "SPP1 - CD44"             = "ECM & Adhesion",
  "THBS1 - CD36"            = "ECM & Adhesion",
  "THBS1 - CD47"            = "ECM & Adhesion",
  "CD99 - CD99"             = "ECM & Adhesion",
  "CD99 - CD99L2"           = "ECM & Adhesion",
  "PECAM1 - PECAM1"         = "ECM & Adhesion",
  "NECTIN3 - NECTIN2"       = "ECM & Adhesion",
  "JAM2 - (ITGAV+ITGB1)"    = "ECM & Adhesion",
  "NCAM1 - FGFR1"           = "ECM & Adhesion",

  # === Angiogenesis & Growth (血管生成/生长因子) ===
  "PGF - VEGFR1"                    = "Angiogenesis & Growth",
  "ANGPT2 - (ITGA5+ITGB1)"          = "Angiogenesis & Growth",
  "ANGPT2 - TEK"                    = "Angiogenesis & Growth",
  "PDGFC - PDGFRA"                  = "Angiogenesis & Growth",
  "FGF10 - FGFR1"                   = "Angiogenesis & Growth",
  "FGF7 - FGFR1"                    = "Angiogenesis & Growth",
  "IGF1 - (ITGA6+ITGB4)"            = "Angiogenesis & Growth",
  "IGF1 - IGF1R"                    = "Angiogenesis & Growth",

  # === Development & Patterning (发育/模式发生) ===
  "BMP5 - (BMPR1A+ACVR2A)"  = "Development",
  "BMP5 - (BMPR1A+BMPR2)"   = "Development",
  "BMP6 - (BMPR1A+ACVR2A)"  = "Development",
  "BMP6 - (BMPR1A+BMPR2)"   = "Development",
  "BMP7 - (BMPR1A+ACVR2A)"  = "Development",
  "BMP7 - (BMPR1A+BMPR2)"   = "Development",
  "DLK1 - NOTCH2"           = "Development",
  "JAG1 - NOTCH2"           = "Development",
  "SEMA3D - (NRP1+PLXNA2)"  = "Development",
  "SEMA3D - (NRP1+PLXNA4)"  = "Development",
  "SEMA3D - (NRP2+PLXNA2)"  = "Development",
  "SEMA3D - (NRP2+PLXNA4)"  = "Development",
  "SEMA6A - PLXNA2"         = "Development",
  "SEMA6A - PLXNA4"         = "Development",
  "PTPRM - PTPRM"           = "Development",

  # === Phagocytosis & Clearance (吞噬/胞葬) ===
  "GAS6 - MERTK"            = "Phagocytosis",
  "GAS6 - AXL"              = "Phagocytosis",
  "PROS1 - AXL"             = "Phagocytosis"
)

cat_colors <- c(
  "Immune Regulation"       = "#1B9E77",
  "Inflammation"            = "#D95F02",
  "ECM & Adhesion"          = "#7570B3",
  "Angiogenesis & Growth"   = "#E7298A",
  "Development"             = "#66A61E",
  "Phagocytosis"            = "#E6AB02"
)

cat_labels <- c(
  "Immune Regulation"       = "Immune Regulation\n(免疫调控)",
  "Inflammation"            = "Inflammation\n(炎症/趋化)",
  "ECM & Adhesion"          = "ECM & Adhesion\n(基质/黏附)",
  "Angiogenesis & Growth"   = "Angiogenesis\n& Growth\n(血管生成/生长)",
  "Development"             = "Development\n(发育/模式)",
  "Phagocytosis"            = "Phagocytosis\n(吞噬/胞葬)"
)

# ── Helper: add category strip to a bubble plot ──
add_category_strip <- function(gg, lr_map, colors, labels, strip_width = 0.04) {
  y_levels <- levels(gg$data$interaction_name_2)

  y_df <- data.frame(
    interaction_name_2 = y_levels,
    y_pos = seq_along(y_levels),
    stringsAsFactors = FALSE
  )
  y_df$category <- sapply(y_df$interaction_name_2, function(x) {
    if (x %in% names(lr_map)) lr_map[[x]] else "Other"
  })
  y_df$category <- factor(y_df$category, levels = names(colors))

  # Place strip on the far left of the plot
  # The bubble plot's x-axis is categorical; we add a dummy column before the first source
  n_sources <- length(levels(gg$data$source))
  strip_data <- y_df
  strip_data$xmin <- 0.4           # left of first source column
  strip_data$xmax <- 0.4 + strip_width * n_sources  # narrow strip

  gg_out <- gg +
    geom_rect(
      data = strip_data,
      aes(xmin = xmin, xmax = xmax, ymin = y_pos - 0.5, ymax = y_pos + 0.5,
          fill = category),
      inherit.aes = FALSE, alpha = 0.85
    ) +
    scale_fill_manual(
      values = colors,
      labels = labels[names(colors)],
      name = "Functional Category",
      drop = FALSE
    ) +
    coord_cartesian(clip = "off") +
    theme(
      plot.margin = margin(t = 5, r = 5, b = 5, l = 0),
      legend.title = element_text(face = "bold", size = 8),
      legend.text  = element_text(size = 7),
      legend.key.size = unit(0.35, "cm")
    )

  return(gg_out)
}

# ═══════════════════════════════════════════════════════════
# INCOMING: other cells → HB
# ═══════════════════════════════════════════════════════════
cat("\n=== Generating INCOMING bubble (other cells -> HB) ===\n")
pdf(NULL)
p_data_in <- netVisual_bubble(
  cellchat,
  sources.use = ct_levels,
  targets.use = "HB",
  remove.isolate = FALSE,
  return.data = TRUE
)
dev.off()

gg_in <- p_data_in$gg.obj

# Adjust source labels to include "-> HB"
gg_in$data$source <- factor(
  paste0(gg_in$data$source, " → HB"),
  levels = paste0(ct_levels, " → HB")
)

gg_in_final <- add_category_strip(gg_in, lr_category, cat_colors, cat_labels, strip_width = 0.055)

n_lr_in <- length(levels(gg_in$data$interaction_name_2))
h_in <- max(6, n_lr_in * 0.15)
w_in <- 9

ggsave("figures/Fig4/Fig4_cellchat_bubble_incoming_category.png",
       gg_in_final, width = w_in, height = h_in, dpi = 300, bg = "white")
cat(sprintf("  Saved: incoming (%d LR pairs, %.1f x %.1f in)\n", n_lr_in, w_in, h_in))

# ═══════════════════════════════════════════════════════════
# OUTGOING: HB → other cells
# ═══════════════════════════════════════════════════════════
cat("\n=== Generating OUTGOING bubble (HB -> other cells) ===\n")
pdf(NULL)
p_data_out <- netVisual_bubble(
  cellchat,
  sources.use = "HB",
  targets.use = ct_levels,
  remove.isolate = FALSE,
  return.data = TRUE
)
dev.off()

gg_out <- p_data_out$gg.obj

# Adjust target labels to include "HB ->"
gg_out$data$target <- factor(
  paste0("HB → ", gg_out$data$target),
  levels = paste0("HB → ", ct_levels)
)

gg_out_final <- add_category_strip(gg_out, lr_category, cat_colors, cat_labels, strip_width = 0.055)

n_lr_out <- length(levels(gg_out$data$interaction_name_2))
h_out <- max(6, n_lr_out * 0.14)
w_out <- 9

ggsave("figures/Fig4/Fig4_cellchat_bubble_outgoing_category.png",
       gg_out_final, width = w_out, height = h_out, dpi = 300, bg = "white")
cat(sprintf("  Saved: outgoing (%d LR pairs, %.1f x %.1f in)\n", n_lr_out, w_out, h_out))

# ═══════════════════════════════════════════════════════════
# Summary: count LR pairs per category
# ═══════════════════════════════════════════════════════════
cat("\n=== Category Summary ===\n")
for (direction in c("Incoming", "Outgoing")) {
  p_data <- if (direction == "Incoming") p_data_in else p_data_out
  y_lvls <- levels(p_data$gg.obj$data$interaction_name_2)
  cats <- sapply(y_lvls, function(x) if (x %in% names(lr_category)) lr_category[[x]] else "Other")
  cat(sprintf("\n%s:\n", direction))
  for (cn in names(cat_colors)) {
    n <- sum(cats == cn)
    cat(sprintf("  %-25s: %d\n", cn, n))
  }
}
cat("\nDone.\n")
