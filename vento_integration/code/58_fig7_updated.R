#!/usr/bin/env Rscript
# Fig7: Disease disruption of dual-track model — updated with corrected groups & IRF1
library(Seurat); library(ggplot2); library(dplyr); library(patchwork); library(tidyr)
set.seed(42)

INPUT <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/seurat_labeled_10datasets.rds"
FIGDIR <- "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures"

seu <- readRDS(INPUT)
labels <- read.csv(file.path(dirname(INPUT), "per_cell_disease_labels.csv"))
detail <- labels$disease_detail[1:ncol(seu)]

# ── Disease groups: PE-focused (Fig7 main), others → supplementary ──
disease_all <- c("Normal_Early","Miscarriage","Infection","Normal_Late","Early_PE","Late_PE","Preterm")
pe_groups <- c("Normal_Late","Early_PE","Late_PE")

dc <- rep(NA_character_, ncol(seu))
dc[seu$dataset %in% c("E-MTAB-12421","E-MTAB-6701") | (seu$dataset=="E-MTAB-12795" & detail=="normal")] <- "Normal_Early"
dc[seu$dataset=="GSE214607"] <- "Miscarriage"
dc[detail %in% c("toxoplasmosis","listeriosis","Plasmodium malariae malaria")] <- "Infection"
dc[(seu$dataset=="GSE290578" & detail=="Normal") | (seu$dataset=="GSE298602" & detail=="Control")] <- "Normal_Late"
dc[seu$dataset=="GSE290578" & detail=="PE"] <- "Early_PE"
dc[detail %in% c("PreE_SF","gHTN","GSE173193","GSE298119")] <- "Late_PE"
dc[detail %in% c("PTL","PTNL")] <- "Preterm"

keep <- dc %in% pe_groups
seu <- subset(seu, cells=colnames(seu)[keep])
seu <- JoinLayers(seu)
seu$disease <- factor(dc[keep], levels=pe_groups)

disease_colors <- c("Normal_Late"="#66C2A5","Early_PE"="#FC8D59","Late_PE"="#8DA0CB")

for(g in levels(seu$disease)) cat(sprintf("%s: %d  ", g, sum(seu$disease==g))); cat("\n")

# ── 7a: ECM vs Immune UMAP ──
ecm_genes <- c("FN1","SPP1","COL1A1","MMP9")
imm_genes <- c("IL1B","TNF","CXCL8","CD44","CD47")
ecm_score <- rowMeans(FetchData(seu, vars=intersect(ecm_genes, rownames(seu))))
imm_score <- rowMeans(FetchData(seu, vars=intersect(imm_genes, rownames(seu))))
seu$UMAP_1 <- Embeddings(seu,"umap")[,1]; seu$UMAP_2 <- Embeddings(seu,"umap")[,2]
umap_df <- data.frame(UMAP_1=seu$UMAP_1, UMAP_2=seu$UMAP_2, ECM=ecm_score, Immune=imm_score)

p1 <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color=ECM)) + geom_point(size=0.12) +
  scale_color_gradientn(colors=c("grey90","#BDD7E7","#6BAED6","#2171B5","#08306B")) +
  labs(title="ECM module") + theme_void() + theme(plot.title=element_text(face="bold",size=18,hjust=0.5), legend.position="right")
p2 <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color=Immune)) + geom_point(size=0.12) +
  scale_color_gradientn(colors=c("grey90","#FCBBA1","#FC9272","#CB181D","#67000D")) +
  labs(title="Immune module") + theme_void() + theme(plot.title=element_text(face="bold",size=18,hjust=0.5), legend.position="right")
ggsave(file.path(FIGDIR,"Fig7a_two_arms.png"), p1|p2, w=12, h=5.5, dpi=300, bg="white")
cat("7a done\n")

# ── 7b: TF expression — PE groups ──
main_tfs <- c("CEBPA","IRF1")
supp_tfs <- c("STAT1","STAT3","KLF4")

for(target_tf in c(main_tfs, supp_tfs)) {
  tf_data <- FetchData(seu, vars=c(target_tf, "disease")) %>% na.omit(); colnames(tf_data)[1] <- "expr"
  normal_vals <- tf_data$expr[tf_data$disease=="Normal_Late"]
  tf_pval <- data.frame()
  # EP vs NL, EP vs LP, LP vs NL
  comp_pairs <- list(c("Early_PE","Normal_Late"), c("Early_PE","Late_PE"), c("Late_PE","Normal_Late"))
  for(cp in comp_pairs) {
    v1 <- tf_data$expr[tf_data$disease==cp[1]]; v2 <- tf_data$expr[tf_data$disease==cp[2]]
    if(length(v1)>10 && length(v2)>10) {
      pv <- tryCatch(wilcox.test(v1, v2, exact=FALSE)$p.value, error=function(e) NA)
      if(!is.na(pv)) tf_pval <- rbind(tf_pval, data.frame(d1=cp[1], d2=cp[2], pval=pv))
    }
  }
  if(nrow(tf_pval)>0) {
    tf_pval$label <- ifelse(tf_pval$pval<0.001,"***",ifelse(tf_pval$pval<0.01,"**",ifelse(tf_pval$pval<0.05,"*","ns")))
    tf_pval$ypos <- max(tf_data$expr) * (1.05 + 0.10*(1:nrow(tf_pval)))
    tf_pval$x_start <- match(tf_pval$d1, pe_groups); tf_pval$x_end <- match(tf_pval$d2, pe_groups)
  }
  p <- ggplot(tf_data, aes(x=disease, y=expr, fill=disease)) +
    geom_violin(scale="width", linewidth=0.3, alpha=0.7) +
    geom_boxplot(width=0.04, outlier.size=0.02, alpha=0.5, fill="white", linewidth=0.2) +
    scale_fill_manual(values=disease_colors) +
    labs(y=sprintf("%s expression", target_tf), x="", title=target_tf) +
    theme_classic(11) + theme(plot.title=element_text(face="bold",size=14,hjust=0.5), legend.position="none",
      axis.text.x=element_text(angle=0,hjust=0.5,size=11,face="bold",color="black"))
  if(nrow(tf_pval)>0) p <- p + geom_segment(data=tf_pval, aes(x=x_start,xend=x_end,y=ypos,yend=ypos), inherit.aes=FALSE, linewidth=0.6, color="black") +
    geom_text(data=tf_pval, aes(x=(x_start+x_end)/2, y=ypos*1.02, label=label), inherit.aes=FALSE, size=3, fontface="bold")
  prefix <- if(target_tf %in% main_tfs) "Fig7b_" else "Fig7S_"
  ggsave(file.path(FIGDIR,sprintf("%s%s.png", prefix, target_tf)), p, w=4.5, h=4.5, dpi=300, bg="white")
}
cat("7b done\n")

# ── 7c: Module scores — PE groups ──
for(target_mod in c("ECM","Immune")) {
  mod_data <- data.frame(score=if(target_mod=="ECM") ecm_score else imm_score, disease=seu$disease) %>% na.omit()
  normal_vals <- mod_data$score[mod_data$disease=="Normal_Late"]
  mod_pval <- data.frame()
  # EP vs NL, EP vs LP, LP vs NL
  comp_pairs <- list(c("Early_PE","Normal_Late"), c("Early_PE","Late_PE"), c("Late_PE","Normal_Late"))
  for(cp in comp_pairs) {
    v1 <- mod_data$score[mod_data$disease==cp[1]]; v2 <- mod_data$score[mod_data$disease==cp[2]]
    if(length(v1)>10 && length(v2)>10) {
      pv <- tryCatch(wilcox.test(v1, v2, exact=FALSE)$p.value, error=function(e) NA)
      if(!is.na(pv)) mod_pval <- rbind(mod_pval, data.frame(d1=cp[1], d2=cp[2], pval=pv))
    }
  }
  if(nrow(mod_pval)>0) {
    mod_pval$label <- ifelse(mod_pval$pval<0.001,"***",ifelse(mod_pval$pval<0.01,"**",ifelse(mod_pval$pval<0.05,"*","ns")))
    mod_pval$ypos <- max(mod_data$score) * (1.05 + 0.10*(1:nrow(mod_pval)))
    mod_pval$x_start <- match(mod_pval$d1, pe_groups); mod_pval$x_end <- match(mod_pval$d2, pe_groups)
  }
  p <- ggplot(mod_data, aes(x=disease, y=score, fill=disease)) +
    geom_violin(scale="width", linewidth=0.3, alpha=0.7) +
    geom_boxplot(width=0.04, outlier.size=0.02, alpha=0.5, fill="white", linewidth=0.2) +
    scale_fill_manual(values=disease_colors) +
    labs(y=sprintf("%s score", target_mod), x="", title=sprintf("%s module", target_mod)) +
    theme_classic(11) + theme(plot.title=element_text(face="bold",size=14,hjust=0.5), legend.position="none",
      axis.text.x=element_text(angle=0,hjust=0.5,size=11,face="bold",color="black"))
  if(nrow(mod_pval)>0) p <- p + geom_segment(data=mod_pval, aes(x=x_start,xend=x_end,y=ypos,yend=ypos), inherit.aes=FALSE, linewidth=0.6, color="black") +
    geom_text(data=mod_pval, aes(x=(x_start+x_end)/2, y=ypos*1.02, label=label), inherit.aes=FALSE, size=3, fontface="bold")
  ggsave(file.path(FIGDIR,sprintf("Fig7c_%s.png", target_mod)), p, w=4.5, h=4.5, dpi=300, bg="white")
}
cat("7c done\n")

# ── 7d: ECM vs Immune scatter — PE groups ──
mod_df <- data.frame(ECM=ecm_score, Immune=imm_score, disease=seu$disease) %>% na.omit()
mod_sample <- do.call(rbind, lapply(split(mod_df, mod_df$disease), function(d) d[sample(nrow(d), min(800,nrow(d))),]))
p <- ggplot(mod_sample, aes(x=ECM, y=Immune, color=disease)) + geom_point(size=0.8, alpha=0.5) +
  scale_color_manual(values=disease_colors) + stat_ellipse(aes(group=disease), linewidth=0.9) +
  labs(x="ECM module score", y="Immune module score", title="ECM vs Immune — PE groups") +
  theme_classic(11) + theme(plot.title=element_text(face="bold",size=13,hjust=0.5), legend.position="right",
    axis.title.x=element_text(face="bold",size=12), axis.title.y=element_text(face="bold",size=12),
    axis.text=element_text(face="bold",size=10,color="black"))
ggsave(file.path(FIGDIR,"Fig7d_scatter.png"), p, w=7, h=6, dpi=300, bg="white")
cat("7d done\n")

# ── 7e: TF-subtype-disease heatmap ──
tfs_7e <- c("CEBPA","IRF1","STAT3","STAT1","KLF4")
tf_sub_dis <- FetchData(seu, vars=c(tfs_7e, "subtype_pred", "disease")) %>% na.omit()
colnames(tf_sub_dis)[colnames(tf_sub_dis)=="subtype_pred"] <- "subtype"
tf_long <- tf_sub_dis %>% pivot_longer(all_of(tfs_7e), names_to="TF", values_to="expr") %>%
  group_by(subtype, disease, TF) %>% summarise(mean_expr=mean(expr), .groups="drop")
tf_long$TF <- factor(tf_long$TF, levels=tfs_7e)

p7e <- ggplot(tf_long, aes(x=disease, y=subtype, fill=mean_expr)) +
  geom_tile(color="white", linewidth=0.5) + geom_text(aes(label=sprintf("%.2f",mean_expr)), size=2.5, color="black") +
  scale_fill_gradientn(colors=c("#4575B4","white","#D73027"), name="Mean expr") +
  facet_wrap(~TF, nrow=2, scales="free") + labs(title="TF-subtype-disease vulnerability map", x="", y="") +
  theme_minimal(10) + theme(axis.text.x=element_text(angle=30,hjust=1,size=9,face="bold",color="black"),
    axis.text.y=element_text(size=9,face="bold",color="black"), plot.title=element_text(face="bold",size=14,hjust=0.5,color="black"),
    strip.text=element_text(face="bold",size=10))
ggsave(file.path(FIGDIR,"Fig7e_heatmap.png"), p7e, w=14, h=8, dpi=300, bg="white")
cat("Fig7 all done\n")
