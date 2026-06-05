#!/usr/bin/env Rscript
# Figure 2a: Pseudotime trajectory on Harmony UMAP using Slingshot
library(Seurat); library(slingshot); library(ggplot2); library(dplyr)

seu <- readRDS('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/Hofbauer_Atlas_Final.rds')

# ---- Subset normal Early+Late ----
cond <- read.csv("/tmp/gse290578_conditions.csv", stringsAsFactors=FALSE)
rownames(cond) <- cond$barcode; seu$gse <- NA
for(bc in colnames(seu)) if(bc %in% cond$barcode) seu$gse[bc] <- cond[bc,"condition"]
seu$tri <- NA
seu$tri[seu$dataset=="Arutyunyan"] <- "Early"
seu$tri[seu$dataset=="GSE290578" & seu$gse=="Normal"] <- "Late"
seu$tri <- factor(seu$tri, levels=c("Early","Late"))
seu_normal <- subset(seu, cells=colnames(seu)[!is.na(seu$tri)])

# ---- Get Harmony UMAP coordinates ----
umap_coords <- Embeddings(seu_normal, "umap")
cat(sprintf("Cells: %d, UMAP dims: %d\n", nrow(umap_coords), ncol(umap_coords)))

# ---- Run Slingshot ----
# Use subtype as cluster labels for trajectory
cl <- seu_normal$subtype
sds <- slingshot(umap_coords, cl, start.clus="Pro-inflammatory",
                 stretch=2, approx_points=200)

cat(sprintf("Lineages found: %d\n", length(slingLineages(sds))))

# Extract pseudotime (use the first/primary lineage)
pt <- slingPseudotime(sds)
# Use the first non-NA lineage
pt_col <- which(colSums(!is.na(pt)) > 100)[1]
pseudotime <- pt[, pt_col]
cat(sprintf("Using lineage %d, %d cells with pseudotime\n", pt_col, sum(!is.na(pseudotime))))

# Add to Seurat
seu_normal$pseudotime <- pseudotime

# ---- Plot: UMAP colored by pseudotime ----
tri_cols <- c("Early"="#4575B4","Late"="#D73027")

# Pseudotime UMAP
umap_df <- as.data.frame(umap_coords)
# Keep original Seurat UMAP column names
colnames(umap_df) <- c("umap_1","umap_2")
umap_df$pseudotime <- seu_normal$pseudotime
umap_df$tri <- seu_normal$tri

p1 <- ggplot(umap_df, aes(x=umap_1, y=umap_2)) +
  geom_point(aes(color=pseudotime), size=0.3, alpha=0.8) +
  scale_color_viridis_c(option="D", name="Pseudotime", na.value="grey80") +
  theme_bw() + theme(panel.grid=element_blank(), panel.border=element_rect(color="black", linewidth=0.4),
    legend.position="right") +
  labs(title="Pseudotime Trajectory", x="UMAP1", y="UMAP2")

# Add curve overlay
curves <- slingCurves(sds)
p1 <- p1 + 
  geom_path(data=as.data.frame(curves[[1]]$s[curves[[1]]$ord, ]),
            aes(x=umap_1, y=umap_2), color="black", linewidth=0.5)

# Trimester UMAP
p2 <- ggplot(umap_df, aes(x=umap_1, y=umap_2)) +
  geom_point(aes(color=tri), size=0.3, alpha=0.8) +
  scale_color_manual(values=tri_cols) +
  theme_bw() + theme(panel.grid=element_blank(), panel.border=element_rect(color="black", linewidth=0.4),
    legend.position="right", legend.title=element_blank()) +
  labs(title="Trimester", x="UMAP1", y="UMAP2")

FIGDIR <- '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'
ggsave(file.path(FIGDIR,"Fig2a_pseudotime.png"), p1, w=5.5, h=4.5, dpi=300, bg="white")
ggsave(file.path(FIGDIR,"Fig2a_trimester.png"), p2, w=5.5, h=4.5, dpi=300, bg="white")
cat("Saved Fig2a_pseudotime.png + Fig2a_trimester.png\n")

# ---- Module scores for Fig2b ----
seu_normal <- AddModuleScore(seu_normal,
  features=list(
    c("TREM2","AXL","CEBPA","ID2","CD5L","NOTCH2"),
    c("TIMP1","VIM","MMP14","ENO1","PGK1","COL1A2","FN1","PDGFB","SOD2"),
    c("CCL8","IL18","IFI30","HLA-DRA","CTSS","FCGR3A")
  ),
  name=c("Dev","Remodeling","Immunity"),
  ctrl=50
)
colnames(seu_normal@meta.data)[grep("Dev\\d",colnames(seu_normal@meta.data))] <- "Score_Dev"
colnames(seu_normal@meta.data)[grep("Remodeling\\d",colnames(seu_normal@meta.data))] <- "Score_Remodeling"
colnames(seu_normal@meta.data)[grep("Immunity\\d",colnames(seu_normal@meta.data))] <- "Score_Immunity"

# Module score vs pseudotime scatter with smoothed lines
meta <- seu_normal@meta.data[!is.na(seu_normal$pseudotime),]
meta_long <- tidyr::pivot_longer(meta, cols=c("Score_Dev","Score_Remodeling","Score_Immunity"),
                                  names_to="Module", values_to="Score")
meta_long$Module <- gsub("Score_","", meta_long$Module)

mod_cols <- c("Dev"="#4575B4","Remodeling"="#2D8B57","Immunity"="#D73027")

p3 <- ggplot(meta_long, aes(x=pseudotime, y=Score, color=Module, fill=Module)) +
  geom_smooth(method="loess", span=0.3, alpha=0.15, linewidth=0.8) +
  scale_color_manual(values=mod_cols) +
  scale_fill_manual(values=mod_cols) +
  labs(x="Pseudotime", y="Module Score", title="Functional Modules along Pseudotime") +
  theme_bw() + theme(panel.grid=element_blank(), panel.border=element_rect(color="black", linewidth=0.4),
    legend.position="top", legend.title=element_blank())

ggsave(file.path(FIGDIR,"Fig2b_module_scores.png"), p3, w=5.5, h=4, dpi=300, bg="white")
cat("Saved Fig2b_module_scores.png\n")

saveRDS(seu_normal, file.path(FIGDIR, "../results/Hofbauer_normal_pseudotime.rds"))
cat("Done\n")
