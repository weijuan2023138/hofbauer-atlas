#!/usr/bin/env Rscript
# Fig7: Disease disruption of the dual-track TF-communication model
setwd("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration")
suppressMessages(library(Seurat)); library(ggplot2); library(dplyr); library(patchwork); library(tidyr)

so <- readRDS("results/Hofbauer_Atlas_Final.rds")
DefaultAssay(so) <- "RNA"

# Define disease groups
disease_order <- c("Normal_Early","Normal_Late","PE","PTB","Miscarriage","Infection")
# Split into early/late groups
disease_early <- c("Normal_Early","Miscarriage","Infection")
disease_late  <- c("Normal_Late","PE","PTB")
disease_map <- list(
  "Normal_Early" = c("Normal 1st trimester","Normal 1st/2nd/Term"),
  "Normal_Late"  = c("Normal 3rd trimester / Preeclampsia"),
  "PE"           = c("Preeclampsia"),
  "PTB"          = c("Preterm Labor","Preterm No Labor","Term Labor"),
  "Miscarriage"  = c("Miscarriage / Normal"),
  "Infection"    = c("Infection")
)
so$disease_simple <- NA
for(d in names(disease_map)) {
  so$disease_simple[so$disease_group %in% disease_map[[d]]] <- d
}
so$disease_simple <- factor(so$disease_simple, levels=disease_order)

disease_colors <- c("Normal_Early"="#4575B4","Normal_Late"="#66C2A5","PE"="#FC8D59",
                     "PTB"="#E41A1C","Miscarriage"="#D73027","Infection"="#FDB462")
sty_title <- theme(plot.title=element_text(face="bold",size=13,hjust=0.5))
sty_void <- theme_void()

# Subtype colors
subtypes <- c("Vascular remodeling","MHCII+ Antigen-presenting","Pro-inflammatory",
              "Homeostatic","PRKN+ Autophagy","MKI67+ Proliferating")
subtype_colors <- c("Vascular remodeling"="#D73027","MHCII+ Antigen-presenting"="#4575B4",
                     "Pro-inflammatory"="#FC8D59","Homeostatic"="#66C2A5",
                     "PRKN+ Autophagy"="#8DA0CB","MKI67+ Proliferating"="grey70")

dir.create("figures/Fig7", showWarnings=FALSE)

# ===== 7a: Two-arm gene expression landscape (UMAP colored by ECM vs Immune) =====
ecm_genes  <- c("FN1","SPP1","COL1A1","MMP9")
imm_genes  <- c("IL1B","TNF","CXCL8","CD44","CD47")
ecm_score <- rowMeans(FetchData(so, vars=ecm_genes[ecm_genes %in% rownames(so)]))
imm_score <- rowMeans(FetchData(so, vars=imm_genes[imm_genes %in% rownames(so)]))

umap_df <- data.frame(UMAP_1=so$UMAP_1, UMAP_2=so$UMAP_2, ECM=ecm_score, Immune=imm_score)

p1 <- ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2, color=ECM)) + geom_point(size=0.15) +
  scale_color_gradientn(colors=c("grey90","#BDD7E7","#6BAED6","#2171B5","#08306B")) +
  labs(title="ECM module") + sty_void +
  theme(plot.title=element_text(face="bold",size=18,hjust=0.5), legend.position="right")
p2 <- ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2, color=Immune)) + geom_point(size=0.15) +
  scale_color_gradientn(colors=c("grey90","#FCBBA1","#FC9272","#CB181D","#67000D")) +
  labs(title="Immune module") + sty_void +
  theme(plot.title=element_text(face="bold",size=18,hjust=0.5), legend.position="right")

p7a <- p1 | p2
ggsave("figures/Fig7/Fig7a_two_arms.png", p7a, w=12, h=5.5, dpi=300, bg="white")

# ===== 7b: TF expression × disease — CEBPA/NFKB1 main, STAT3/STAT1/RELB supplementary =====
main_tfs <- c("CEBPA","NFKB1")
supp_tfs <- c("STAT3","STAT1","RELB")

for(target_tf in c(main_tfs, supp_tfs)) {
  tf_data <- FetchData(so, vars=c(target_tf, "disease_simple")) %>% na.omit()
  colnames(tf_data)[1] <- "expr"
  
  for(grp_name in c("early","late")) {
    grp_diseases <- if(grp_name=="early") disease_early else disease_late
    tf_sub <- tf_data %>% filter(disease_simple %in% grp_diseases)
    tf_sub$disease_simple <- factor(tf_sub$disease_simple, levels=grp_diseases)
    
    normal_vals <- tf_sub$expr[tf_sub$disease_simple==grp_diseases[1]]
    tf_pval <- data.frame()
    for(d in grp_diseases[-1]) {
      dis_vals <- tf_sub$expr[tf_sub$disease_simple==d]
      if(length(dis_vals) > 10 && length(normal_vals) > 10) {
        pv <- tryCatch(wilcox.test(normal_vals, dis_vals, exact=FALSE)$p.value, error=function(e) NA)
        if(!is.na(pv)) tf_pval <- rbind(tf_pval, data.frame(disease=d, pval=pv))
      }
    }
    if(nrow(tf_pval) > 0) {
      tf_pval$label <- ifelse(tf_pval$pval<0.001, "***", ifelse(tf_pval$pval<0.01, "**", ifelse(tf_pval$pval<0.05, "*", "ns")))
      tf_pval$txt_size <- ifelse(tf_pval$label=="ns", 2.8, 3.8)
      tf_pval$txt_face <- "bold"
      tf_pval$ypos <- max(tf_sub$expr) * (1.05 + 0.10 * (1:nrow(tf_pval)))
      # Horizontal significance lines from Normal to each disease
      norm_x <- 1  # Normal is first
      for(i in 1:nrow(tf_pval)) {
        tf_pval$x_start[i] <- norm_x
        tf_pval$x_end[i]  <- i + 1
      }
    }
    
    p <- ggplot(tf_sub, aes(x=disease_simple, y=expr, fill=disease_simple)) +
      geom_violin(scale="width", linewidth=0.3, alpha=0.7) +
      geom_boxplot(width=0.04, outlier.size=0.02, alpha=0.5, fill="white", linewidth=0.2, fatten=1.5) +
      
      scale_fill_manual(values=disease_colors) +
      labs(y=sprintf("%s expression", target_tf), x="",
           title=sprintf("%s — %s gestation", target_tf, grp_name)) +
      theme_classic(base_size=11) + sty_title +
      theme(legend.position="none", axis.text.x=element_text(angle=30,hjust=1,size=10,face="bold"))
    
    if(nrow(tf_pval) > 0) {
      p <- p +
        geom_segment(data=tf_pval, aes(x=x_start, xend=x_end, y=ypos, yend=ypos),
                     inherit.aes=FALSE, linewidth=0.6, color="black") +
        geom_text(data=tf_pval, aes(x=(x_start+x_end)/2, y=ypos*1.02, label=label),
                  inherit.aes=FALSE, size=tf_pval$txt_size, fontface=tf_pval$txt_face)
    }
    out_dir <- if(target_tf %in% main_tfs) "Fig7" else "FigS"
    ggsave(sprintf("figures/%s/Fig7b_%s_%s.png", out_dir, target_tf, grp_name), p, w=4.5, h=4.5, dpi=300, bg="white")
  }
}

# ===== 7c: Module scores × disease — grouped early/late =====
for(target_mod in c("ECM","Immune")) {
  mod_data <- data.frame(score=if(target_mod=="ECM") ecm_score else imm_score, disease=so$disease_simple) %>% na.omit()
  
  for(grp_name in c("early","late")) {
    grp_diseases <- if(grp_name=="early") disease_early else disease_late
    mod_sub <- mod_data %>% filter(disease %in% grp_diseases)
    mod_sub$disease <- factor(mod_sub$disease, levels=grp_diseases)
    
    normal_vals <- mod_sub$score[mod_sub$disease==grp_diseases[1]]
    mod_pval <- data.frame()
    for(d in grp_diseases[-1]) {
      dis_vals <- mod_sub$score[mod_sub$disease==d]
      if(length(dis_vals) > 10 && length(normal_vals) > 10) {
        pv <- tryCatch(wilcox.test(normal_vals, dis_vals, exact=FALSE)$p.value, error=function(e) NA)
        if(!is.na(pv)) mod_pval <- rbind(mod_pval, data.frame(disease=d, pval=pv))
      }
    }
    if(nrow(mod_pval) > 0) {
      mod_pval$label <- ifelse(mod_pval$pval<0.001, "***", ifelse(mod_pval$pval<0.01, "**", ifelse(mod_pval$pval<0.05, "*", "ns")))
      mod_pval$txt_size <- ifelse(mod_pval$label=="ns", 2.8, 3.8)
      mod_pval$txt_face <- "bold"
      mod_pval$ypos <- max(mod_sub$score) * (1.05 + 0.10 * (1:nrow(mod_pval)))
      norm_x <- 1
      for(i in 1:nrow(mod_pval)) { mod_pval$x_start[i] <- norm_x; mod_pval$x_end[i] <- i+1 }
    }
    
    p <- ggplot(mod_sub, aes(x=disease, y=score, fill=disease)) +
      geom_violin(scale="width", linewidth=0.3, alpha=0.7) +
      geom_boxplot(width=0.04, outlier.size=0.02, alpha=0.5, fill="white", linewidth=0.2, fatten=1.5) +
      
      scale_fill_manual(values=disease_colors) +
      labs(y=sprintf("%s score", target_mod), x="", title=sprintf("%s — %s gestation", target_mod, grp_name)) +
      theme_classic(base_size=11) + sty_title +
      theme(legend.position="none", axis.text.x=element_text(angle=30,hjust=1,size=10,face="bold"))
    
    if(nrow(mod_pval) > 0) {
      p <- p +
        geom_segment(data=mod_pval, aes(x=x_start, xend=x_end, y=ypos, yend=ypos),
                     inherit.aes=FALSE, linewidth=0.6, color="black") +
        geom_text(data=mod_pval, aes(x=(x_start+x_end)/2, y=ypos*1.02, label=label),
                  inherit.aes=FALSE, size=mod_pval$txt_size, fontface=mod_pval$txt_face)
    }
    ggsave(sprintf("figures/Fig7/Fig7c_%s_%s.png", target_mod, grp_name), p, w=4.5, h=4.5, dpi=300, bg="white")
  }
}

# ===== 7d: ECM vs Immune scatter — grouped early/late =====
mod_df <- data.frame(ECM=ecm_score, Immune=imm_score, disease=so$disease_simple) %>% na.omit()
set.seed(42)

for(grp_name in c("early","late")) {
  grp_diseases <- if(grp_name=="early") disease_early else disease_late
  mod_sub <- mod_df %>% filter(disease %in% grp_diseases)
  mod_sub$disease <- factor(mod_sub$disease, levels=grp_diseases)
  
  # Downsample
  mod_sample <- do.call(rbind, lapply(split(mod_sub, mod_sub$disease), function(d) {
    d[sample(nrow(d), min(800, nrow(d))), ]
  }))
  
  p <- ggplot(mod_sample, aes(x=ECM, y=Immune, color=disease)) +
    geom_point(size=0.8, alpha=0.5) +
    scale_color_manual(values=disease_colors) +
    stat_ellipse(aes(group=disease), linewidth=0.9) +
    labs(x="ECM module score", y="Immune module score",
         title=sprintf("ECM vs Immune module landscape — %s gestation", grp_name)) +
    theme_classic(base_size=11) + sty_title +
    theme(legend.position="right",
          axis.title.x=element_text(face="bold",size=12),
          axis.title.y=element_text(face="bold",size=12),
          axis.text=element_text(face="bold",size=10,color="black"))
  
  ggsave(sprintf("figures/Fig7/Fig7d_%s.png", grp_name), p, w=7, h=6, dpi=300, bg="white")
}

# ===== 7e: TF-subtype-disease — single fused heatmap =====
tf_sub_dis <- FetchData(so, vars=c("CEBPA","NFKB1","STAT3","STAT1","RELB","subtype","disease_simple")) %>%
  na.omit()

tf_long <- tf_sub_dis %>%
  pivot_longer(c(CEBPA,NFKB1,STAT3,STAT1,RELB), names_to="TF", values_to="expr") %>%
  group_by(subtype, disease_simple, TF) %>%
  summarise(mean_expr=mean(expr), .groups="drop")

tf_long$TF <- factor(tf_long$TF, levels=c("CEBPA","STAT3","NFKB1","STAT1","RELB"))

p7e <- ggplot(tf_long, aes(x=disease_simple, y=subtype, fill=mean_expr)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=sprintf("%.2f", mean_expr)), size=2.5) +
  scale_fill_gradientn(colors=c("#4575B4","white","#D73027"), name="Mean expr") +
  facet_wrap(~TF, nrow=2, scales="free") +
  labs(title="TF-subtype-disease vulnerability map", x="", y="") +
  theme_minimal(base_size=10) +
  theme(axis.text.x=element_text(angle=30,hjust=1,size=9,face="bold"),
        axis.text.y=element_text(size=9,face="bold"),
        plot.title=element_text(face="bold",size=14,hjust=0.5),
        strip.text=element_text(face="bold",size=10))

ggsave("figures/Fig7/Fig7e_heatmap.png", p7e, w=14, h=8, dpi=300, bg="white")

message("Fig7 done — disease disruption of dual-track model")
