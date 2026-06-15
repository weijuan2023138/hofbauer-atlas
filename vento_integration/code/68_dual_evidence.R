#!/usr/bin/env Rscript
# Multi-modality evidence: Vascular↔ECM, MHCII+↔Immune
library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
set.seed(42)

FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

# ====== 1. scRNA-seq: per-cell correlation ======
seu <- readRDS("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds")
labels <- read.csv("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/per_cell_disease_labels.csv")
detail <- labels$disease_detail[1:ncol(seu)]
dc <- rep(NA_character_, ncol(seu))
dc[(seu$dataset=="GSE290578" & detail=="Normal") | (seu$dataset=="GSE298602" & detail=="Control")] <- "Normal_Late"
dc[seu$dataset=="GSE290578" & detail=="PE"] <- "Early_PE"
dc[detail %in% c("PreE_SF","gHTN","GSE173193","GSE298119")] <- "Late_PE"
seu <- subset(seu, cells=colnames(seu)[dc %in% c("Normal_Late","Early_PE","Late_PE")])
seu$group <- factor(dc[dc %in% c("Normal_Late","Early_PE","Late_PE")], levels=c("Normal_Late","Early_PE","Late_PE"))
seu <- JoinLayers(seu)

# Add subtype module scores
ref <- readRDS("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds")
Idents(ref) <- "subtype_pred"
markers_list <- list(
  "Vascular remodeling" = c("FN1","COL1A2","SPP1","FLT1","PAPPA","AEBP1","NOTUM","SERPINE2","TIMP3","PTN"),
  "MHCII+ Antigen-presenting" = c("HLA-DRA","HLA-DRB1","CD74","FCGR3A","HLA-DPA1","HLA-DQB1")
)
for(st in names(markers_list)) {
  genes <- intersect(markers_list[[st]], rownames(seu))
  seu <- AddModuleScore(seu, features=list(genes), name=make.names(st), ctrl=min(100,nrow(seu)))
  old <- paste0(make.names(st), "1")
  if(old %in% colnames(seu@meta.data)) colnames(seu@meta.data)[colnames(seu@meta.data)==old] <- paste0("Module_",st)
}

# Per-cell scatter: Module vs gene expression
df_sc <- FetchData(seu, vars=c("Module_Vascular remodeling","Module_MHCII+ Antigen-presenting","FN1","FLT1","HLA-DRA","CD74","group"))
df_sc <- df_sc[sample(nrow(df_sc), min(8000,nrow(df_sc))),]

p1a <- ggplot(df_sc, aes(x=`Module_Vascular remodeling`, y=FN1, color=group)) +
  geom_point(size=0.4, alpha=0.4) + geom_smooth(method="lm", se=FALSE, linewidth=1.2) +
  scale_color_manual(values=c(Normal_Late="#66C2A5",Early_PE="#FC8D59",Late_PE="#8DA0CB")) +
  labs(title="scRNA-seq: Vascular module vs FN1", x="Vascular remodeling module", y="FN1 expression") +
  theme_bw(11) + theme(plot.title=element_text(face="bold",size=12,hjust=0.5), legend.position="none")

p1b <- ggplot(df_sc, aes(x=`Module_Vascular remodeling`, y=FLT1, color=group)) +
  geom_point(size=0.4, alpha=0.4) + geom_smooth(method="lm", se=FALSE, linewidth=1.2) +
  scale_color_manual(values=c(Normal_Late="#66C2A5",Early_PE="#FC8D59",Late_PE="#8DA0CB")) +
  labs(title="scRNA-seq: Vascular module vs FLT1", x="Vascular remodeling module", y="FLT1 expression") +
  theme_bw(11) + theme(plot.title=element_text(face="bold",size=12,hjust=0.5), legend.position="none")

p1c <- ggplot(df_sc, aes(x=`Module_MHCII+ Antigen-presenting`, y=HLA.DRA, color=group)) +
  geom_point(size=0.4, alpha=0.4) + geom_smooth(method="lm", se=FALSE, linewidth=1.2) +
  scale_color_manual(values=c(Normal_Late="#66C2A5",Early_PE="#FC8D59",Late_PE="#8DA0CB")) +
  labs(title="scRNA-seq: MHCII+ module vs HLA-DRA", x="MHCII+ module", y="HLA-DRA expression") +
  theme_bw(11) + theme(plot.title=element_text(face="bold",size=12,hjust=0.5), legend.position="none")

p1d <- ggplot(df_sc, aes(x=`Module_MHCII+ Antigen-presenting`, y=CD74, color=group)) +
  geom_point(size=0.4, alpha=0.4) + geom_smooth(method="lm", se=FALSE, linewidth=1.2) +
  scale_color_manual(values=c(Normal_Late="#66C2A5",Early_PE="#FC8D59",Late_PE="#8DA0CB")) +
  labs(title="scRNA-seq: MHCII+ module vs CD74", x="MHCII+ module", y="CD74 expression") +
  theme_bw(11) + theme(plot.title=element_text(face="bold",size=12,hjust=0.5), legend.position="none")

sc_evidence <- (p1a|p1b)/(p1c|p1d)

# ====== 2. Visium: per-spot correlation ======
SPDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/spatial"
visium_results <- data.frame()
for(sid in c("Pla_Camb9518737","Pla_HDBR9518710","WS_PLA_S9101764","WS_PLA_S9101765")) {
  rf <- file.path(SPDIR, paste0("spatial_modules_", sid, ".rds"))
  if(!file.exists(rf)) next
  sp <- readRDS(rf)
  vasc <- sp@meta.data[["Spatial_Vascular remodeling"]]
  mhcii <- sp@meta.data[["Spatial_MHCII+ Antigen-presenting"]]
  for(g in c("FN1","FLT1","HLA-DRA","CD74")) {
    if(!g %in% rownames(sp)) next
    expr <- FetchData(sp, g)[,1]
    r_vasc <- cor(vasc, expr, use="pairwise"); r_mhc <- cor(mhcii, expr, use="pairwise")
    visium_results <- rbind(visium_results, data.frame(Sample=sid, Gene=g, Vasc_r=r_vasc, MHCII_r=r_mhc))
  }
}
cat("\n=== Visium spatial correlation ===\n")
print(visium_results, row.names=FALSE)

# ====== 3. STOMICS: per-spot correlation ======
# (already computed in Python, re-run key comparisons)
cat("\n=== STOMICS spatial correlation ===\n")
stomics_genes <- c("FN1","FLT1","HLA-DRA","CD74")
for(g in stomics_genes) {
  r <- tryCatch({
    system(sprintf("cd /home/weijuan/文档/胎盘单细胞数据 && python3 -c \"
import scanpy as sc, numpy as np
ad=sc.read_h5ad('raw_data/UCSF_Li_2026/STOMICS.h5ad',backed='r')
mask=(ad.obs.celltype=='HB')&(ad.obs.sample_id.isin(['001','002','004']))
hb=ad[mask].to_memory()
sc.pp.normalize_total(hb,target_sum=1e4);sc.pp.log1p(hb)
for m in ['Vascular_score','Immune_score']:
  if m not in hb.obs: continue
  r=np.corrcoef(hb.obs[m],hb[:,'%s'].X.toarray().flatten())[0,1]
  print('%.3f'%%r)
\" 2>&1", g), intern=TRUE)
    paste(g, ":", r[1])
  }, error=function(e) paste(g, ": ERROR"))
  cat(sprintf("  %s\n", r))
}

# ====== 4. Summary correlation table ======
# Compute scRNA correlation
sc_r <- data.frame(
  Modality = "scRNA-seq",
  TF_Gene = c("Vascular↔FN1","Vascular↔FLT1","MHCII+↔HLA-DRA","MHCII+↔CD74"),
  Correlation = c(
    cor(df_sc$`Module_Vascular remodeling`, df_sc$FN1),
    cor(df_sc$`Module_Vascular remodeling`, df_sc$FLT1),
    cor(df_sc$`Module_MHCII+ Antigen-presenting`, df_sc$HLA.DRA),
    cor(df_sc$`Module_MHCII+ Antigen-presenting`, df_sc$CD74)
  )
)
visium_summary <- visium_results %>% group_by(Gene) %>% summarise(r=mean(case_when(Gene%in%c("FN1","FLT1")~Vasc_r, TRUE~MHCII_r),na.rm=TRUE), .groups="drop") %>% mutate(Modality="Visium", TF_Gene=paste0(ifelse(Gene%in%c("FN1","FLT1"),"Vascular↔","MHCII+↔"),Gene)) %>% select(Modality, TF_Gene, r)
colnames(visium_summary)[3] <- "Correlation"

cat("\n=== Multi-modality evidence summary ===\n")
print(rbind(sc_r, visium_summary), row.names=FALSE)

ggsave(file.path(FIGDIR,"Fig7S_evidence_scRNA.png"), sc_evidence, w=12, h=10, dpi=300, bg="white")
cat("\nEvidence figures saved\n")
