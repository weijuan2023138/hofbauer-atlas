#!/usr/bin/env Rscript
# Developmental analysis: pseudotime + GO enrichment by trimester
library(Seurat); library(ggplot2); library(dplyr); library(monocle3)

INPUT <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/Hofbauer_Atlas_Final.rds'
OUTDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results'
FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'

seu <- readRDS(INPUT)
cond <- read.csv("/tmp/gse290578_conditions.csv", stringsAsFactors=FALSE)
rownames(cond) <- cond$barcode
seu$gse <- NA
for(bc in colnames(seu)) if(bc %in% cond$barcode) seu$gse[bc] <- cond[bc,"condition"]

# ---- 1. Subset normal cells for development analysis ----
seu$tri <- NA
seu$tri[seu$dataset=="Arutyunyan"] <- "Early"
seu$tri[seu$dataset=="UCSF_Li_2026"] <- "Mid"
seu$tri[seu$dataset=="GSE290578" & seu$gse=="Normal"] <- "Late"
seu$tri <- factor(seu$tri, levels=c("Early","Mid","Late"))

normal_mask <- !is.na(seu$tri)
seu_normal <- subset(seu, cells=colnames(seu)[normal_mask])
cat("Normal cells for development:", ncol(seu_normal), "\n")
cat("Trimester distribution:\n")
print(table(seu_normal$tri))

# ---- 2. DEG by trimester (for GO) ----
Idents(seu_normal) <- "tri"

# Early vs Others, Mid vs Others, Late vs Others
deg_list <- list()
deg_all <- data.frame()
for(tri in levels(seu_normal$tri)) {
  deg <- FindMarkers(seu_normal, ident.1=tri, ident.2=NULL, 
                     only.pos=TRUE, min.pct=0.1, logfc.threshold=0.25, test.use="t")
  deg$gene <- rownames(deg)
  deg$trimester <- tri
  deg_list[[tri]] <- deg
  deg_all <- rbind(deg_all, deg)
  cat(sprintf("%s: %d DEGs\n", tri, nrow(deg)))
}

write.csv(deg_all, file.path(OUTDIR, "dev_trimester_DEGs.csv"), row.names=FALSE)
cat("Saved DEGs\n")

# ---- 3. Monocle3 pseudotime ----
# Convert to CellDataSet
expr <- GetAssayData(seu_normal, assay="RNA", layer="data")
cds <- new_cell_data_set(expr, cell_metadata=seu_normal@meta.data, gene_metadata=data.frame(
  gene_short_name=rownames(seu_normal), row.names=rownames(seu_normal)))

# Use DEGs + highly variable genes as ordering genes
ordering_genes <- unique(c(deg_all$gene, VariableFeatures(seu_normal)))
ordering_genes <- intersect(ordering_genes, rownames(cds))
cat("Ordering genes:", length(ordering_genes), "\n")

cds <- preprocess_cds(cds, num_dim=30)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# Root at Early GW4.5 cells
early_cells <- colnames(cds)[colData(cds)$tri=="Early"]
root_cells <- early_cells[1:min(500, length(early_cells))]
cds <- order_cells(cds, root_cells=root_cells)

# Save
saveRDS(cds, file.path(OUTDIR, "monocle_cds.rds"))
cat("Saved monocle_cds\n")

# ---- 4. Pseudotime heatmap (top variable genes along pseudotime) ----
# Get top genes changing along pseudotime
top_genes <- deg_all %>% group_by(trimester) %>% top_n(n=15, wt=avg_log2FC) %>% pull(gene) %>% unique()
top_genes <- intersect(top_genes, rownames(seu_normal))
cat("Top genes for heatmap:", length(top_genes), "\n")

# Get pseudotime values
pt <- pseudotime(cds)
pt <- pt[!is.na(pt)]
cells_ordered <- names(sort(pt))

# Heatmap data
hm_data <- GetAssayData(seu_normal, assay="RNA", layer="data")[top_genes, cells_ordered]
# Scale by row
hm_data <- t(scale(t(hm_data)))
hm_data <- pmin(pmax(hm_data, -2), 2)

# Annotations
ann <- data.frame(
  Trimester = seu_normal$tri[cells_ordered],
  Subtype = seu_normal$subtype[cells_ordered],
  Pseudotime = pt[cells_ordered]
)
rownames(ann) <- cells_ordered

# Save for Python plotting or plot with ComplexHeatmap
library(ComplexHeatmap); library(circlize)
st_cols <- c("Pro-inflammatory"="#C62828","MHCII+ Antigen-presenting"="#BF4E1A",
  "Homeostatic"="#1B6B93","PRKN+ Autophagy"="#7B4FA0",
  "Vascular remodeling"="#2D8B57","MKI67+ Proliferating"="#37474F")
tri_cols <- c("Early"="#4DBBD5","Mid"="#7E6148","Late"="#E18727")

ha <- HeatmapAnnotation(
  Trimester = ann$Trimester,
  Subtype = ann$Subtype,
  col = list(Trimester=tri_cols, Subtype=st_cols),
  show_legend=TRUE
)

png(file.path(FIGDIR,"Fig_dev_heatmap.png"), w=14, h=8, units="in", res=300, bg="white")
draw(Heatmap(hm_data, name="Expression", show_row_names=TRUE, show_column_names=FALSE,
  row_names_gp=gpar(fontsize=7), cluster_columns=FALSE, cluster_rows=TRUE,
  top_annotation=ha, use_raster=FALSE,
  col=colorRamp2(c(-2,0,2), c("#1B6B93","white","#C62828"))))
dev.off()
cat("Saved heatmap\n")

cat("\nDone! Go enrichment needs clusterProfiler\n")
