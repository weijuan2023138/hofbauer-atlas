#!/usr/bin/env Rscript
# Fig7 v2: Fixed — remove CellChat, correct roles, add summary heatmap
setwd("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration")
suppressMessages(library(Seurat)); library(ggplot2); library(dplyr); library(patchwork); library(tidyr)

so <- readRDS("results/Hofbauer_Atlas_Final.rds")
DefaultAssay(so) <- "RNA"

subtypes <- c("Vascular remodeling","MHCII+ Antigen-presenting","Pro-inflammatory",
              "Homeostatic","PRKN+ Autophagy","MKI67+ Proliferating")
so$subtype <- factor(so$subtype, levels=subtypes)

subtype_colors <- c("Vascular remodeling"="#D73027","MHCII+ Antigen-presenting"="#4575B4",
                     "Pro-inflammatory"="#FC8D59","Homeostatic"="#66C2A5",
                     "PRKN+ Autophagy"="#8DA0CB","MKI67+ Proliferating"="grey70")

mod_col <- c("ECM"="#D73027","TGFβ"="#FC8D59","Growth Factor"="#4575B4","Immune"="#91BFDB")
sty_title <- theme(plot.title=element_text(face="bold",size=13,hjust=0.5))

# ===== Fig7c FIXED: Corrected role classification =====
key_genes <- c("FN1","SPP1","IL1B","CD44","TGFB1","MKI67","CEBPA","NFKB1","STAT3")
key_genes <- key_genes[key_genes %in% rownames(so)]
expr_means <- FetchData(so, vars=c(key_genes, "subtype")) %>%
  group_by(subtype) %>% summarise(across(all_of(key_genes), mean), .groups="drop")

# Manual role assignment based on dominant functional profile
role_df <- expr_means %>% mutate(
  ECM_score   = FN1 + SPP1,
  Immune_score = IL1B + CD44,
  role = case_when(
    subtype == "Vascular remodeling"     ~ "ECM Sender",
    subtype == "MHCII+ Antigen-presenting" ~ "Immune Receiver",
    subtype == "Pro-inflammatory"        ~ "Dual (ECM+Immune)",
    subtype == "Homeostatic"             ~ "TGFβ Dominant",
    subtype == "PRKN+ Autophagy"         ~ "Immune Receiver",
    subtype == "MKI67+ Proliferating"    ~ "Silent"
  )
)

role_colors <- c("ECM Sender"="#D73027","Dual (ECM+Immune)"="#7B3294",
                 "Immune Receiver"="#4575B4","TGFβ Dominant"="#66C2A5",
                 "Silent"="grey70")

cat("=== Corrected role classification ===\n")
print(role_df[,c("subtype","FN1","SPP1","IL1B","CD44","TGFB1","MKI67","role")])

p7c <- ggplot(role_df, aes(x=ECM_score, y=Immune_score, fill=role)) +
  geom_point(size=8, shape=21, color="black", stroke=0.5) +
  geom_text(aes(label=subtype), vjust=-1.2, size=3, fontface="bold") +
  scale_fill_manual(values=role_colors) +
  labs(x="ECM communication (FN1+SPP1)", y="Immune communication (IL1B+CD44)",
       title="Subtype communication specialization") +
  xlim(0, max(role_df$ECM_score)*1.4) + ylim(0, max(role_df$Immune_score)*1.3) +
  theme_classic(base_size=11) + sty_title +
  theme(legend.position="right")
ggsave("figures/Fig7/Fig7c_role.png", p7c, w=7.5, h=6, dpi=300, bg="white")

# ===== Fig7f NEW: Summary heatmap — key functional genes per subtype =====
summary_genes <- c("FN1","SPP1","COL1A1","MMP9","IL1B","TNF","CXCL8","CD44","CD47",
                    "TGFB1","IGF1","ITGAV","ITGB1","CEBPA","STAT3","NFKB1","MKI67")
summary_genes <- summary_genes[summary_genes %in% rownames(so)]

# Mean expression per subtype, z-score normalized
mat <- FetchData(so, vars=c(summary_genes, "subtype")) %>%
  group_by(subtype) %>% summarise(across(all_of(summary_genes), mean), .groups="drop")
mat_z <- mat %>% mutate(across(-subtype, ~ scale(.)[,1]))
mat_long <- pivot_longer(mat_z, -subtype, names_to="gene", values_to="zscore")
mat_long$subtype <- factor(mat_long$subtype, levels=subtypes)

# Gene categories
gene_cat <- data.frame(
  gene = summary_genes,
  category = case_when(
    summary_genes %in% c("FN1","SPP1","COL1A1","MMP9") ~ "ECM",
    summary_genes %in% c("IL1B","TNF","CXCL8","CD44","CD47") ~ "Immune",
    summary_genes %in% c("TGFB1","IGF1","ITGAV","ITGB1") ~ "Signaling",
    summary_genes %in% c("CEBPA","STAT3","NFKB1") ~ "TF",
    TRUE ~ "Proliferation"
  )
)
cat_colors <- c("ECM"="#D73027","Immune"="#4575B4","Signaling"="#FC8D59","TF"="#66C2A5","Proliferation"="grey70")

p7f <- ggplot(mat_long, aes(x=subtype, y=gene, fill=zscore)) +
  geom_tile(color="white", linewidth=0.3) +
  scale_fill_gradient2(low="#4575B4", mid="white", high="#D73027", midpoint=0, name="Z-score") +
  labs(x="", y="", title="Functional gene specialization across HB subtypes") +
  theme_minimal(base_size=10) +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=9,face="bold"),
        axis.text.y=element_text(size=8),
        plot.title=element_text(face="bold",size=13,hjust=0.5))
ggsave("figures/Fig7/Fig7f_summary_heatmap.png", p7f, w=10, h=7, dpi=300, bg="white")

# Delete old CellChat panels
old_cellchat <- list.files("figures/Fig7", pattern="Fig7f_chord|Fig7f_heatmap|Fig7f_network", full.names=TRUE)
unlink(old_cellchat)

message("Fig7 fixed: roles corrected, CellChat removed, summary heatmap added")
