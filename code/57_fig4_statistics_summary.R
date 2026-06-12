#!/usr/bin/env Rscript
# 57_fig4_statistics_summary.R
# Summary statistics for CellChat: interactions per cell type + pathway activity

library(CellChat); library(dplyr); library(ggplot2); library(patchwork)

cellchat <- readRDS("results/cellchat_ucsf_mid.rds")

ct_levels <- c("CD14_M","dNK","FB","fEC","SCT","VCT","vEC")

# ── Panel A: Number of significant interactions (incoming + outgoing) ──
# Incoming to HB
inc_raw <- subsetCommunication(cellchat, targets.use="HB")
inc_count <- inc_raw %>%
  filter(source %in% ct_levels) %>%
  group_by(source) %>%
  summarise(incoming=n_distinct(interaction_name), .groups="drop") %>%
  rename(celltype=source)

# Outgoing from HB
out_raw <- subsetCommunication(cellchat, sources.use="HB")
out_count <- out_raw %>%
  filter(target %in% ct_levels) %>%
  group_by(target) %>%
  summarise(outgoing=n_distinct(interaction_name), .groups="drop") %>%
  rename(celltype=target)

# Combine
interaction_counts <- full_join(inc_count, out_count, by="celltype") %>%
  tidyr::pivot_longer(cols=c(incoming, outgoing),
                      names_to="direction", values_to="count") %>%
  mutate(celltype=factor(celltype, levels=ct_levels),
         direction=factor(direction, levels=c("incoming","outgoing")))

p_interactions <- ggplot(interaction_counts,
  aes(x=celltype, y=count, fill=direction)) +
  geom_col(position="dodge", width=0.7, color="black", linewidth=0.3) +
  scale_fill_manual(values=c("incoming"="#5B9BD5", "outgoing"="#ED7D31"),
                    labels=c("Incoming \u2192 HB", "Outgoing \u2190 HB")) +
  theme_minimal(base_size=11) +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_text(face="bold", size=11, color="black"),
    axis.text.x=element_text(angle=45, hjust=1, size=10, face="bold", color="black"),
    axis.text.y=element_text(size=10, color="black"),
    panel.grid.major.x=element_blank(),
    panel.grid.minor=element_blank(),
    panel.border=element_rect(fill=NA, color="black", linewidth=0.5),
    legend.position=c(0.7, 0.85),
    legend.title=element_blank(),
    legend.text=element_text(size=9, face="bold", color="black"),
    legend.background=element_rect(fill="white", color="black", linewidth=0.3),
    plot.title=element_text(face="bold", size=12, hjust=0.5, color="black")
  ) +
  labs(y="Number of significant L-R pairs", 
       title="Significant interactions with Hofbauer cells")

ggsave("figures/Fig4/Fig4B_cellchat_statistics_a_interactions.png",
  p_interactions, width=6, height=5, dpi=300, bg="white")

# ── Panel B: Top pathways by total communication probability ──
# Incoming pathways
inc_pathway <- inc_raw %>%
  filter(source %in% ct_levels) %>%
  group_by(pathway_name) %>%
  summarise(total_prob=sum(prob), n_pairs=n_distinct(interaction_name), .groups="drop") %>%
  arrange(desc(total_prob)) %>%
  head(15) %>%
  mutate(direction="Incoming")

# Outgoing pathways
out_pathway <- out_raw %>%
  filter(target %in% ct_levels) %>%
  group_by(pathway_name) %>%
  summarise(total_prob=sum(prob), n_pairs=n_distinct(interaction_name), .groups="drop") %>%
  arrange(desc(total_prob)) %>%
  head(15) %>%
  mutate(direction="Outgoing")

pathway_summary <- bind_rows(inc_pathway, out_pathway) %>%
  mutate(direction=factor(direction, levels=c("Incoming","Outgoing")))

p_pathways <- ggplot(pathway_summary,
  aes(x=reorder(pathway_name, total_prob), y=total_prob, fill=direction)) +
  geom_col(position="dodge", width=0.7, color="black", linewidth=0.3) +
  coord_flip() +
  scale_fill_manual(values=c("Incoming"="#5B9BD5", "Outgoing"="#ED7D31")) +
  theme_minimal(base_size=11) +
  theme(
    axis.title.y=element_blank(),
    axis.title.x=element_text(face="bold", size=11, color="black"),
    axis.text.y=element_text(size=9, face="bold", color="black"),
    axis.text.x=element_text(size=10, color="black"),
    panel.grid.minor=element_blank(),
    panel.border=element_rect(fill=NA, color="black", linewidth=0.5),
    legend.position=c(0.75, 0.2),
    legend.title=element_blank(),
    legend.text=element_text(size=9, face="bold", color="black"),
    legend.background=element_rect(fill="white", color="black", linewidth=0.3),
    plot.title=element_text(face="bold", size=12, hjust=0.5, color="black")
  ) +
  labs(x="Pathway", y="Total communication probability",
       title="Top 15 pathways by communication strength")

ggsave("figures/Fig4/补充图Fig4A_cellchat_statistics_b_pathways.png",
  p_pathways, width=8, height=6, dpi=300, bg="white")

# ── Panel C: Expanded pathway detail, grouped by functional category ──
func_cat <- c(
  "TGFb"="TGF-β/BMP", "BMP"="TGF-β/BMP", "GDF"="TGF-β/BMP",
  "FGF"="Growth Factor", "VEGF"="Growth Factor", "IGF"="Growth Factor",
  "PDGF"="Growth Factor", "ANGPT"="Growth Factor", "PROS"="Growth Factor",
  "GAS"="Growth Factor", "VISFATIN"="Growth Factor",
  "COLLAGEN"="ECM/Matrix", "LAMININ"="ECM/Matrix", "FN1"="ECM/Matrix",
  "THBS"="ECM/Matrix", "TENASCIN"="ECM/Matrix", "SPP1"="ECM/Matrix",
  "PERIOSTIN"="ECM/Matrix",
  "CD45"="Immune", "CD86"="Immune", "CD99"="Immune",
  "MHC-I"="Immune", "MHC-II"="Immune", "ADGRE5"="Immune",
  "MIF"="Immune", "CD46"="Immune", "NECTIN"="Immune",
  "VCAM"="Immune", "JAM"="Immune", "PECAM1"="Immune",
  "NCAM"="Immune", "APP"="Immune",
  "NOTCH"="Developmental", "SEMA3"="Developmental",
  "SEMA6"="Developmental", "PTPRM"="Developmental", "CHEMERIN"="Developmental",
  "DLK1"="Developmental", "GAS6"="Developmental"
)

cat_cols <- c(
  "TGF-β/BMP"      = "#E7298A",
  "Growth Factor"   = "#D95F02",
  "ECM/Matrix"      = "#7570B3",
  "Immune"          = "#1B9E77",
  "Developmental"   = "#66A61E"
)

all_inc <- inc_raw %>% filter(source %in% ct_levels) %>%
  mutate(category=func_cat[pathway_name], direction="Incoming")
all_out <- out_raw %>% filter(target %in% ct_levels) %>%
  mutate(category=func_cat[pathway_name], direction="Outgoing")

all_interactions <- bind_rows(all_inc, all_out) %>%
  mutate(category=ifelse(is.na(category), "Other", category))

# Count L-R pairs per pathway × direction
pathway_lr <- all_interactions %>%
  group_by(pathway_name, category, direction) %>%
  summarise(n_pairs=n_distinct(interaction_name), .groups="drop")

# Order pathways by total pairs within category
pathway_order <- pathway_lr %>%
  group_by(pathway_name, category) %>%
  summarise(total=sum(n_pairs), .groups="drop") %>%
  arrange(category, desc(total))

pathway_lr <- pathway_lr %>%
  mutate(pathway_name=factor(pathway_name, levels=pathway_order$pathway_name),
         category=factor(category, levels=names(cat_cols)),
         direction=factor(direction, levels=c("Incoming","Outgoing")))

p_categories <- ggplot(pathway_lr,
  aes(x=n_pairs, y=pathway_name, fill=direction)) +
  geom_col(position=position_dodge(width=0.72), width=0.58,
           alpha=0.85, color=NA) +
  geom_text(aes(label=n_pairs, group=direction),
            position=position_dodge(width=0.72),
            hjust=-0.25, size=2.5, color="grey40", fontface="plain") +
  facet_grid(category ~ ., scales="free_y", space="free_y", switch="y") +
  scale_fill_manual(values=c("Incoming"="#5B9BD5", "Outgoing"="#ED7D31"),
                    labels=c("Incoming \u2192 HB", "Outgoing \u2190 HB")) +
  scale_x_continuous(expand=expansion(mult=c(0, 0.2))) +
  theme_minimal(base_size=11) +
  theme(
    axis.title.y=element_blank(),
    axis.title.x=element_text(face="bold", size=11, color="black"),
    axis.text.y=element_text(size=8.5, face="bold", color="black"),
    axis.text.x=element_text(size=10, color="black"),
    axis.ticks.y=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor=element_blank(),
    strip.background=element_rect(fill="grey95", color="grey80", linewidth=0.4),
    strip.text.y.left=element_text(size=10, face="bold", angle=0,
                                   margin=margin(r=6)),
    strip.placement="outside",
    legend.position="top",
    legend.title=element_blank(),
    legend.text=element_text(size=9, face="bold", color="black"),
    legend.background=element_rect(fill="white", color="black", linewidth=0.3),
    legend.key.size=unit(0.5, "cm"),
    plot.title=element_text(face="bold", size=12, hjust=0.5, color="black", margin=margin(b=6)),
    plot.subtitle=element_text(size=8.5, hjust=0.5, color="black",
                               margin=margin(b=10)),
    panel.spacing=unit(0.6, "lines")
  ) +
  labs(x="Number of significant L-R pairs",
       title="Ligand-receptor interactions by pathway",
       subtitle=paste0(n_distinct(pathway_lr$pathway_name), " pathways across 5 functional categories  |  ",
                       sum(pathway_lr$n_pairs), " total L-R pairs"))

ggsave("figures/Fig4/补充图Fig4B_cellchat_statistics_c_categories.png",
  p_categories, width=7.5, height=9, dpi=300, bg="white")

cat("Done: Fig4_cellchat_statistics_a/b/c → 3 individual PNGs\n")

# Print summary
cat("\n=== Interaction counts ===\n")
print(interaction_counts %>% arrange(celltype, direction))

cat("\n=== Top incoming pathways ===\n")
print(inc_pathway %>% select(pathway_name, total_prob, n_pairs))

cat("\n=== Top outgoing pathways ===\n")
print(out_pathway %>% select(pathway_name, total_prob, n_pairs))
