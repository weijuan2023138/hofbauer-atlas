#!/usr/bin/env Rscript
# Ridge plot v2: biology-driven gene selection + Nature color palette
library(Seurat); library(ggplot2); library(dplyr); library(ggridges)

seu <- readRDS('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/Hofbauer_Atlas_Final.rds')
cond <- read.csv("/tmp/gse290578_conditions.csv", stringsAsFactors=FALSE)
rownames(cond) <- cond$barcode
seu$gse <- NA
for(bc in colnames(seu)) if(bc %in% cond$barcode) seu$gse[bc] <- cond[bc,"condition"]

seu$tri <- NA
seu$tri[seu$dataset=="Arutyunyan"] <- "Early"
seu$tri[seu$dataset=="UCSF_Li_2026"] <- "Mid"
seu$tri[seu$dataset=="GSE290578" & seu$gse=="Normal"] <- "Late"
seu$tri <- factor(seu$tri, levels=c("Early","Mid","Late"))
seu_normal <- subset(seu, cells=colnames(seu)[!is.na(seu$tri)])

# ---- Gene selection by biological module ----
early_genes <- c(
  "TREM2",      # Hofbauer/microglia identity, phagocytosis
  "AXL",        # efferocytosis, tissue-resident macrophage
  "CEBPA",      # myeloid lineage TF
  "ID2",        # developmental TF, progenitor maintenance
  "CD36",       # scavenger receptor, lipid uptake
  "TIMP1",      # tissue remodeling, MMP inhibitor
  "FOLR2"       # Hofbauer-specific, folate receptor
)

mid_genes <- c(
  "RHOA",       # small GTPase, cytoskeleton
  "ITGAV",      # integrin, cell adhesion
  "TGFB1",      # tissue homeostasis, anti-inflammatory
  "VIM",        # intermediate filament, mesenchymal
  "ANGPT2"      # angiogenesis, vascular remodeling
)

late_genes <- c(
  "CD74",       # MHC-II chaperone, antigen processing
  "HLA-DRB1",   # MHC-II antigen presentation
  "S100A8",     # alarmin, inflammation
  "S100A9",     # alarmin, forms calprotectin
  "CCL2",       # monocyte/macrophage chemotaxis
  "NFKB1",      # NF-kB subunit, inflammatory signaling
  "CXCL8"       # neutrophil chemokine (IL-8)
)

genes <- c(early_genes, mid_genes, late_genes)
genes <- intersect(genes, rownames(seu_normal))
cat(sprintf("%d genes kept\n", length(genes)))

# ---- Extract expression ----
expr <- FetchData(seu_normal, vars=c(genes, "tri"))
expr_long <- tidyr::pivot_longer(expr, cols=all_of(genes),
                                  names_to="Gene", values_to="Expression")

module_map <- setNames(
  c(rep("Differentiation &\nDevelopment", length(early_genes)),
    rep("Tissue Remodeling\n& Signaling", length(mid_genes)),
    rep("Immune\nMaturation", length(late_genes))),
  c(early_genes, mid_genes, late_genes)
)
expr_long$Module <- module_map[expr_long$Gene]
expr_long$Module <- factor(expr_long$Module, levels=unique(module_map))

# Keep gene order within module
expr_long$Gene <- factor(expr_long$Gene, levels=intersect(genes, rownames(seu_normal)))

# ---- Nature-inspired developmental palette ----
tri_cols <- c("Early"="#4575B4", "Mid"="#FDAE61", "Late"="#D73027")

p <- ggplot(expr_long, aes(x=Expression, y=Gene, fill=tri)) +
  geom_density_ridges(alpha=0.75, scale=1.1, rel_min_height=0.01,
                       color=NA, panel_scaling=FALSE) +
  scale_fill_manual(values=tri_cols, name="Trimester") +
  facet_grid(Module ~ ., scales="free_y", space="free_y", switch="y") +
  labs(x="Expression Level", y="") +
  theme_bw() +
  theme(
    text=element_text(family="sans"),
    axis.text=element_text(color="black", size=10),
    axis.text.y=element_text(size=10, face="italic"),
    axis.title.x=element_text(size=11),
    legend.position="top",
    legend.title=element_text(size=10),
    legend.text=element_text(size=9),
    legend.key.size=unit(0.4, "cm"),
    panel.grid.major=element_line(color="grey92", linewidth=0.3),
    panel.grid.minor=element_blank(),
    panel.border=element_rect(color="black", linewidth=0.5),
    strip.background=element_rect(fill="white", color="black", linewidth=0.5),
    strip.text=element_text(size=10.5, face="bold"),
    strip.placement="outside"
  )

FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'
ggsave(file.path(FIGDIR,"Fig_dev_ridge.png"), p, w=7.5, h=8.5, dpi=300, bg="white")
cat("\nSaved Fig_dev_ridge.png\n")
