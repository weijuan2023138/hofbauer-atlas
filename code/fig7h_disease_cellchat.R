#!/usr/bin/env Rscript
# Fig7h: CellChat comparison by disease — Normal vs PE vs PTB vs Miscarriage vs Infection
setwd("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration")
suppressMessages(library(Seurat)); library(CellChat); library(ggplot2); library(dplyr); library(patchwork)

so <- readRDS("results/Hofbauer_Atlas_Final.rds")
DefaultAssay(so) <- "RNA"

subtypes <- c("Vascular remodeling","MHCII+ Antigen-presenting","Pro-inflammatory",
              "Homeostatic","PRKN+ Autophagy","MKI67+ Proliferating")
so$subtype <- factor(so$subtype, levels=subtypes)

# Define disease groups (simplify)
disease_map <- list(
  "Normal" = c("Normal 1st/2nd/Term","Normal 1st trimester","Normal 3rd trimester / Preeclampsia"),
  "PE"     = c("Preeclampsia"),
  "PTB"    = c("Preterm Labor","Preterm No Labor","Term Labor"),
  "Miscarriage" = c("Miscarriage / Normal"),
  "Infection"   = c("Infection")
)

run_cellchat <- function(so_disease, label) {
  if(ncol(so_disease) < 50) return(NULL)
  # Downsample per subtype
  set.seed(42)
  cells_use <- c()
  for(st in subtypes) {
    cidx <- which(so_disease$subtype == st)
    if(length(cidx) > 200) cidx <- sample(cidx, 200)
    if(length(cidx) > 0) cells_use <- c(cells_use, cidx)
  }
  so_sub <- so_disease[, cells_use]
  if(ncol(so_sub) < 30) return(NULL)
  
  data_matrix <- as.matrix(GetAssayData(so_sub, assay="RNA", layer="data"))
  meta <- so_sub@meta.data
  cc <- createCellChat(object=data_matrix, meta=meta, group.by="subtype")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- projectData(cc, PPI.human)
  cc <- computeCommunProb(cc, type="triMean", trim=0.1)
  cc <- filterCommunication(cc, min.cells=5)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  
  cat(sprintf("  %s: %d cells, %d interactions, strength=%.2f\n",
    label, ncol(so_sub), sum(cc@net$count), sum(cc@net$weight)))
  return(cc)
}

results <- list()
for(disease_name in names(disease_map)) {
  cat(sprintf("\n=== %s ===\n", disease_name))
  disease_groups <- disease_map[[disease_name]]
  so_d <- subset(so, disease_group %in% disease_groups)
  cat(sprintf("  Total cells: %d\n", ncol(so_d)))
  print(table(so_d$subtype))
  
  cc <- run_cellchat(so_d, disease_name)
  if(!is.null(cc)) results[[disease_name]] <- cc
}

# Save
saveRDS(results, "results/cellchat_by_disease.rds")

# === Comparison heatmap ===
disease_names <- names(results)
if(length(disease_names) >= 2) {
  # Extract interaction matrices
  net_list <- lapply(results, function(cc) cc@net$count)
  
  # Compare Normal vs each disease
  if("Normal" %in% disease_names) {
    normal_net <- net_list[["Normal"]]
    for(dn in setdiff(disease_names, "Normal")) {
      if(!is.null(net_list[[dn]])) {
        diff_net <- net_list[[dn]] - normal_net
        
        df <- reshape2::melt(diff_net)
        colnames(df) <- c("Sender","Receiver","Diff")
        df$Sender <- factor(df$Sender, levels=subtypes)
        df$Receiver <- factor(df$Receiver, levels=subtypes)
        
        p <- ggplot(df, aes(x=Sender, y=Receiver, fill=Diff)) +
          geom_tile(color="white", linewidth=0.3) +
          scale_fill_gradient2(low="#4575B4", mid="white", high="#D73027", midpoint=0) +
          labs(title=sprintf("%s vs Normal: interaction change", dn), x="Sender", y="Receiver", fill="Δ") +
          theme_minimal(base_size=9) +
          theme(axis.text.x=element_text(angle=45,hjust=1,size=8),
                axis.text.y=element_text(size=8),
                plot.title=element_text(face="bold",size=11,hjust=0.5))
        
        ggsave(sprintf("figures/Fig7/Fig7h_diff_%s.png", dn), p, w=6, h=5, dpi=300, bg="white")
        cat(sprintf("Saved: Fig7h_diff_%s.png\n", dn))
      }
    }
  }
}

message("Fig7h done")
