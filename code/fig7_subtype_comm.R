#!/usr/bin/env Rscript
# Fig7: HB subtype communication — publication quality
setwd("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration")
suppressMessages(library(Seurat)); library(ggplot2); library(dplyr); library(patchwork); library(tidyr)

so <- readRDS("results/Hofbauer_Atlas_Final.rds")
DefaultAssay(so) <- "RNA"

comm_genes <- c("SPP1","FN1","COL1A1","COL1A2","TGFB1","BMP2","TNF","IL1B","IL6","IL10",
  "CXCL8","CCL2","CCL3","CCL4","CXCL2","CSF1","VEGFA","HGF","IGF1","MMP9","THBS1",
  "PTPRM","CD44","CD47","ITGAV","ITGB1","ITGA1","ITGA2","ITGB3","ITGB5")
comm_genes <- comm_genes[comm_genes %in% rownames(so)]

module_map <- list(
  "ECM"           = c("SPP1","FN1","COL1A1","COL1A2","MMP9","THBS1","ITGAV","ITGB1","ITGA1","ITGA2","ITGB3","ITGB5"),
  "TGFβ"          = c("TGFB1","BMP2"),
  "Growth Factor" = c("IGF1","VEGFA","HGF","CSF1"),
  "Immune"        = c("TNF","IL1B","IL6","IL10","CXCL8","CCL2","CCL3","CCL4","CXCL2","CD44","CD47","PTPRM")
)

subtypes <- c("Vascular remodeling","MHCII+ Antigen-presenting","Pro-inflammatory",
              "Homeostatic","PRKN+ Autophagy","MKI67+ Proliferating")
subtype_colors <- c("Vascular remodeling"="#D73027","MHCII+ Antigen-presenting"="#4575B4",
                     "Pro-inflammatory"="#FC8D59","Homeostatic"="#66C2A5",
                     "PRKN+ Autophagy"="#8DA0CB","MKI67+ Proliferating"="grey70")

so$subtype <- factor(so$subtype, levels=subtypes)
Idents(so) <- so$subtype

mod_col <- c("ECM"="#D73027","TGFβ"="#FC8D59","Growth Factor"="#4575B4","Immune"="#91BFDB")

sty_title <- theme(plot.title=element_text(face="bold",size=13,hjust=0.5))

# ====== Fig7a: Dotplot ======
gene_module <- bind_rows(lapply(names(module_map), function(m) {
  data.frame(gene=module_map[[m]][module_map[[m]] %in% comm_genes], module=m)
}))
gene_order <- gene_module$gene

p7a <- DotPlot(so, features=gene_order, group.by="subtype", cols=c("lightgrey","#D73027"), dot.scale=8) +
  RotatedAxis() +
  labs(x="", y="") +
  ggtitle("Communication gene expression across HB subtypes") +
  sty_title +
  theme(axis.text.x=element_text(size=7.5, angle=45, hjust=1),
        axis.text.y=element_text(size=10, face="bold"),
        legend.position="right")

ggsave("figures/Fig7/Fig7a_dotplot.png", p7a, w=16, h=5.5, dpi=300, bg="white")

# ====== Fig7b: Module scores with subtype colors ======
cell_mod <- data.frame(row.names=colnames(so))
for(m in names(module_map)) cell_mod[[m]] <- rowMeans(FetchData(so, vars=module_map[[m]][module_map[[m]] %in% rownames(so)]))
cell_mod$subtype <- so$subtype

p7b_list <- list()
for(m in names(module_map)) {
  df <- cell_mod[,c(m,"subtype")]; colnames(df)[1] <- "score"
  p <- ggplot(df, aes(x=subtype, y=score, fill=subtype)) +
    geom_violin(scale="width", linewidth=0.3, alpha=0.7) +
    geom_boxplot(width=0.08, outlier.size=0.05, alpha=0.5, fill="white", linewidth=0.3) +
    scale_fill_manual(values=subtype_colors) +
    labs(title=m, y="Expression", x="") +
    theme_classic(base_size=9) +
    theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,size=7,face="bold"),
          plot.title=element_text(face="bold",size=10,hjust=0.5,color=mod_col[m]))
  p7b_list[[m]] <- p
}
p7b <- wrap_plots(p7b_list, ncol=4) +
  plot_annotation(title="Communication module scores per subtype", theme=sty_title)
ggsave("figures/Fig7/Fig7b_modules.png", p7b, w=15, h=4, dpi=300, bg="white")

# ====== Fig7c: Communication role — with shapes ======
role_df <- cell_mod %>% group_by(subtype) %>%
  summarise(ECM_send=mean(!!sym("ECM")), Immune_recv=mean(!!sym("Immune")), .groups="drop")
med_ecm <- median(role_df$ECM_send); med_imm <- median(role_df$Immune_recv)
role_df$role <- case_when(
  role_df$ECM_send > med_ecm & role_df$Immune_recv > med_imm ~ "Dual",
  role_df$ECM_send > med_ecm ~ "Sender",
  role_df$Immune_recv > med_imm ~ "Receiver",
  TRUE ~ "Silent"
)

p7c <- ggplot(role_df, aes(x=ECM_send, y=Immune_recv, fill=subtype)) +
  geom_hline(yintercept=med_imm, linetype="dashed", color="grey70", linewidth=0.5) +
  geom_vline(xintercept=med_ecm, linetype="dashed", color="grey70", linewidth=0.5) +
  geom_point(size=8, shape=21, color="black", stroke=0.5) +
  geom_text(aes(label=subtype), vjust=-1.2, size=3.2, fontface="bold") +
  scale_fill_manual(values=subtype_colors) +
  annotate("text", x=max(role_df$ECM_send)*1.25, y=med_imm, label="Receiver", hjust=1, size=3, color="grey50") +
  annotate("text", x=max(role_df$ECM_send)*1.25, y=max(role_df$Immune_recv)*0.95, label="Dual", hjust=1, size=3, color="grey50") +
  annotate("text", x=med_ecm*0.3, y=max(role_df$Immune_recv)*0.95, label="Sender", size=3, color="grey50") +
  labs(x="ECM communication", y="Immune communication", title="Subtype communication roles") +
  xlim(0, max(role_df$ECM_send)*1.5) + ylim(0, max(role_df$Immune_recv)*1.1) +
  theme_classic(base_size=11) + sty_title +
  theme(legend.position="none")
ggsave("figures/Fig7/Fig7c_role.png", p7c, w=7.5, h=6, dpi=300, bg="white")

# ====== Fig7d: TF expression per subtype ======
tfs <- c("CEBPA","STAT3","STAT1","NFKB1","RELB")
tfs <- tfs[tfs %in% rownames(so)]
tf_data <- FetchData(so, vars=c(tfs, "subtype"))
tf_long <- pivot_longer(tf_data, -subtype, names_to="TF", values_to="expr")
tf_long$TF <- factor(tf_long$TF, levels=tfs)

p7d <- ggplot(tf_long, aes(x=subtype, y=expr, fill=subtype)) +
  geom_violin(scale="width", linewidth=0.3, alpha=0.7) +
  geom_boxplot(width=0.1, outlier.size=0.05, alpha=0.5, fill="white", linewidth=0.3) +
  facet_wrap(~TF, nrow=1, scales="free_y") +
  scale_fill_manual(values=subtype_colors) +
  labs(y="Expression", x="") +
  ggtitle("TF expression across HB subtypes") +
  theme_classic(base_size=10) +
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,size=8,face="bold"),
        plot.title=element_text(face="bold",size=13,hjust=0.5),
        strip.text=element_text(face="bold",size=10))
ggsave("figures/Fig7/Fig7d_tf_subtype.png", p7d, w=15, h=4, dpi=300, bg="white")

# ====== Fig7e: Disease proportions ======
prop <- so@meta.data %>% count(subtype, disease_group) %>%
  group_by(disease_group) %>% mutate(pct=n/sum(n)*100) %>% ungroup()
prop$disease_group <- factor(prop$disease_group, levels=c("Normal","PE","PTB"))

p7e <- ggplot(prop, aes(x=subtype, y=pct, fill=disease_group)) +
  geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.6,
           color="black", linewidth=0.3) +
  scale_fill_manual(values=c("Normal"="#4575B4","PE"="#D73027","PTB"="#FC8D59"),
                    name="Disease") +
  labs(y="Proportion (%)", x="", title="Subtype proportions by disease") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=9,face="bold"),
        plot.title=element_text(face="bold",size=13,hjust=0.5),
        legend.position="right")
ggsave("figures/Fig7/Fig7e_proportions.png", p7e, w=10, h=5, dpi=300, bg="white")

message("Fig7 panels saved")
