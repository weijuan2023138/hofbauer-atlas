#!/usr/bin/env Rscript
# 54_fig4_bubble_category.R
# Recreate bubble plots with pathway category annotations on the right side
# Three categories: Cell-Cell Contact (blue), ECM-Receptor (orange), Secreted Signaling (red)

library(CellChat); library(dplyr); library(ggplot2); library(patchwork)

cellchat <- readRDS("results/cellchat_ucsf_mid.rds")
db <- cellchat@DB$interaction
pw_cat <- db %>% select(pathway_name, annotation) %>% distinct() %>% filter(!is.na(annotation))
cat_cols <- c("Cell-Cell Contact"="#4575B4", "ECM-Receptor"="#FDAE61", "Secreted Signaling"="#D73027")

# ── INCOMING (to HB) ──
inc_raw <- subsetCommunication(cellchat, targets.use="HB")
inc_agg <- inc_raw %>%
  group_by(pathway_name, source) %>%
  summarise(prob=sum(prob), .groups="drop") %>%
  rename(celltype=source) %>%
  left_join(pw_cat, by="pathway_name")
inc_agg$annotation[is.na(inc_agg$annotation)] <- "Other"
inc_agg$celltype <- factor(inc_agg$celltype, levels=c("FB","fEC","vEC","VCT","SCT","dNK","CD14_M"))

pw_inc <- inc_agg %>%
  group_by(pathway_name) %>%
  summarise(total=sum(prob), .groups="drop") %>%
  arrange(total) %>%
  pull(pathway_name)
inc_agg$pathway_name <- factor(inc_agg$pathway_name, levels=pw_inc)

p1 <- ggplot(inc_agg, aes(x=celltype, y=pathway_name)) +
  geom_point(aes(size=prob, color=prob), stroke=0) +
  scale_size_continuous(name="Probability", range=c(0.3, 8)) +
  scale_color_gradientn(colors=c("#FFFFCC","#FDAE61","#D73027","#67000D")) +
  theme_minimal(base_size=10) +
  theme(
    axis.title=element_blank(),
    axis.text.x=element_text(angle=45, hjust=1, size=10, face="bold"),
    axis.text.y=element_text(size=9, face="bold"),
    panel.grid=element_blank(),
    panel.border=element_rect(fill=NA, color="black", linewidth=0.5),
    plot.title=element_text(face="bold", size=13, hjust=0.5, margin=margin(b=10)),
    legend.position="right"
  ) +
  labs(title="Signaling to Hofbauer cells") +
  guides(size=guide_legend(order=2, override.aes=list(color="grey50")))

p1_cat <- inc_agg %>%
  distinct(pathway_name, annotation) %>%
  mutate(pathway_name=factor(pathway_name, levels=pw_inc)) %>%
  ggplot(aes(x="", y=pathway_name, fill=annotation)) +
  geom_tile(width=0.6) +
  scale_fill_manual(values=cat_cols) +
  theme_void() +
  theme(legend.position="none")

ggsave("figures/Fig4/Fig4_cellchat_bubble_incoming_v2.png",
  p1 + p1_cat + plot_layout(widths=c(6.5, 0.5)),
  width=8.5, height=8, dpi=300, bg="white")

# ── OUTGOING (from HB) ──
out_raw <- subsetCommunication(cellchat, sources.use="HB")
out_agg <- out_raw %>%
  group_by(pathway_name, target) %>%
  summarise(prob=sum(prob), .groups="drop") %>%
  rename(celltype=target) %>%
  left_join(pw_cat, by="pathway_name")
out_agg$annotation[is.na(out_agg$annotation)] <- "Other"
out_agg$celltype <- factor(out_agg$celltype, levels=c("FB","fEC","vEC","VCT","SCT","dNK","CD14_M"))

pw_out <- out_agg %>%
  group_by(pathway_name) %>%
  summarise(total=sum(prob), .groups="drop") %>%
  arrange(total) %>%
  pull(pathway_name)
out_agg$pathway_name <- factor(out_agg$pathway_name, levels=pw_out)

p2 <- ggplot(out_agg, aes(x=celltype, y=pathway_name)) +
  geom_point(aes(size=prob, color=prob), stroke=0) +
  scale_size_continuous(name="Probability", range=c(0.3, 8)) +
  scale_color_gradientn(colors=c("#FFFFCC","#FDAE61","#D73027","#67000D")) +
  theme_minimal(base_size=10) +
  theme(
    axis.title=element_blank(),
    axis.text.x=element_text(angle=45, hjust=1, size=10, face="bold"),
    axis.text.y=element_text(size=9, face="bold"),
    panel.grid=element_blank(),
    panel.border=element_rect(fill=NA, color="black", linewidth=0.5),
    plot.title=element_text(face="bold", size=13, hjust=0.5, margin=margin(b=10)),
    legend.position="right"
  ) +
  labs(title="Hofbauer signaling to neighbors") +
  guides(size=guide_legend(order=2, override.aes=list(color="grey50")))

p2_cat <- out_agg %>%
  distinct(pathway_name, annotation) %>%
  mutate(pathway_name=factor(pathway_name, levels=pw_out)) %>%
  ggplot(aes(x="", y=pathway_name, fill=annotation)) +
  geom_tile(width=0.6) +
  scale_fill_manual(values=cat_cols) +
  theme_void() +
  theme(legend.position="none")

ggsave("figures/Fig4/Fig4_cellchat_bubble_outgoing_v2.png",
  p2 + p2_cat + plot_layout(widths=c(6.5, 0.5)),
  width=8.5, height=7.5, dpi=300, bg="white")

cat("Done: Fig4_cellchat_bubble_incoming_v2.png + _outgoing_v2.png\n")
