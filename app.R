#!/usr/bin/env Rscript
# Hofbauer Atlas — Interactive Shiny App (disease groups + supplementary data)
# Deploy to shinyapps.io

library(shiny); library(ggplot2); library(dplyr); library(plotly); library(DT)
library(shinythemes); library(tidyr)

cat("Loading data...\n")
meta      <- readRDS("meta.rds")
prop      <- readRDS("proportions.rds")
# Exclude Preterm/Term Labor (batch effect)
prop      <- prop[!prop$disease_group %in% c("Preterm Labor", "Term Labor"), ]
prop$disease_group <- droplevels(prop$disease_group)
gene_expr <- readRDS("gene_expr.rds")
markers   <- readRDS("markers.rds")
mean_all  <- readRDS("mean_expr_all_genes.rds")
mean_ds   <- readRDS("mean_by_subtype_dataset.rds")
full_mat  <- readRDS("full_sparse_matrix.rds")
mean_sd   <- readRDS("mean_by_subtype_disease.rds")  # subtype x disease_group means (all cells)
cat(sprintf("Loaded: %d cells (%d main + %d supp)\n", nrow(meta),
  sum(meta$source=="Main"), sum(meta$source=="Supp")))

# ---- Colors ----
subtype_colors <- c(
  "Homeostatic"="#2E86AB",           # blue
  "CCL13+ Regulatory"="#50C878",     # green
  "C1Q+ Complement"="#F4A460",       # sand
  "PRKN+ Autophagy"="#9370DB",       # purple
  "Pro-inflammatory"="#DC143C",      # red
  "MKI67+ Proliferating"="#FFD700",  # gold
  "SPP1+ Remodeling"="#A0522D",      # brown
  "Not subtyped"="#BEBEBE"           # grey
)

disease_colors <- c(
  "Normal 1st trimester"="#4DBBD5","Normal 3rd trimester"="#00A087",
  "Preterm Labor"="#E18727","Term Labor"="#B09C85",
  "Preeclampsia"="#E64B35","Miscarriage"="#F39B7F",
  "Infection"="#DC0000"
)

all_genes <- setdiff(colnames(gene_expr), c("subtype_final","dataset"))
all_gene_names <- rownames(mean_all)

get_gene_expr <- function(gene) {
  # Only works for main dataset cells
  if(gene %in% all_genes) return(gene_expr[[gene]])
  if(gene %in% rownames(full_mat)) return(as.numeric(full_mat[gene, 1:nrow(gene_expr)]))
  NULL
}

# ---- UI ----
ui <- navbarPage(
  title="Hofbauer Macrophage Atlas",
  theme=shinythemes::shinytheme("flatly"),

  tabPanel("Overview",
    fluidRow(
      column(4,
        h4("Hofbauer Macrophage Atlas"),
        p(sprintf("%d cells from 4 main + 3 supplementary cohorts", nrow(meta))),
        p("7 subtypes across normal pregnancy, preeclampsia, miscarriage, and infection"),
        hr(),
        selectInput("ov_color", "Color by:",
          choices=c("Subtype"="subtype_final","Disease Group"="disease_group","Dataset"="dataset")),
        checkboxInput("ov_split", "Split by disease group", FALSE)
      ),
      column(8, plotlyOutput("ov_umap", height="550px"))
    ),
    hr(),
    h4("Subtype Composition by Disease Group"),
    plotlyOutput("ov_prop", height="400px")
  ),

  tabPanel("Subtype Atlas",
    sidebarLayout(
      sidebarPanel(
        selectInput("st_select", "Subtype:", choices=names(subtype_colors)[names(subtype_colors)!="Not subtyped"]),
        h5("Top Markers"), DTOutput("st_table"), width=3
      ),
      mainPanel(
        h4("Distribution Across Disease Groups"),
        plotlyOutput("st_prop", height="400px"),
        h4("Key Marker Expression (main cohorts)"),
        plotlyOutput("st_violin", height="400px")
      )
    )
  ),

  tabPanel("Gene Lookup",
    sidebarLayout(
      sidebarPanel(
        textInput("gene_search", "Gene:", value="FOLR2",
                  placeholder="Type a gene name, e.g. CXCL16..."),
        verbatimTextOutput("gene_stats"),
        p("Per-cell data from main cohorts. Mean heatmap includes all."), width=3
      ),
      mainPanel(
        fluidRow(
          column(6, h4("UMAP (main data)"), plotlyOutput("gene_umap", height="400px")),
          column(6, h4("By Subtype"), plotlyOutput("gene_violin", height="400px"))
        ),
        h4("Mean Expression by Subtype x Dataset"),
        plotlyOutput("gene_hm", height="300px")
      )
    )
  ),

  tabPanel("Disease Comparison",
    sidebarLayout(
      sidebarPanel(
        selectInput("comp_group_col","Group by:",
          choices=c("Disease Group"="disease_group","Dataset"="dataset","Raw Disease"="disease")),
        selectInput("comp_g1","Group A:",choices=NULL),
        selectInput("comp_g2","Group B:",choices=NULL),
        width=3
      ),
      mainPanel(
        h4("Subtype Proportions"),
        plotlyOutput("comp_prop", height="400px"),
        h4("UMAP"),
        plotlyOutput("comp_umap", height="500px")
      )
    )
  ),

  tabPanel("Download",
    h4("Data Tables"),
    DTOutput("dl_table"),
    br(),
    downloadButton("dl_meta","Cell Metadata (CSV)"),
    downloadButton("dl_markers","Subtype Markers (CSV)"),
    downloadButton("dl_prop","Proportions (CSV)")
  )
)

# ---- Server ----
server <- function(input, output, session) {

  output$ov_umap <- renderPlotly({
    col <- input$ov_color
    colors <- if(col=="subtype_final") subtype_colors else if(col=="disease_group") disease_colors else NULL
    df <- meta[!is.na(meta$UMAP_1),]  # 显示所有有 UMAP 的细胞
    p <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=.data[[col]],
               text=paste(subtype_final,"<br>",disease_group))) +
      geom_point(size=0.3, alpha=0.6) +
      {if(!is.null(colors)) scale_color_manual(values=colors, name="")} +
      theme_minimal(base_size=12) + labs(x="UMAP 1", y="UMAP 2") +
      theme(legend.position="bottom")
    if(input$ov_split) p <- p + facet_wrap(~disease_group, nrow=2)
    ggplotly(p, tooltip="text") %>% layout(legend=list(orientation="h", y=-0.2))
  })

  output$ov_prop <- renderPlotly({
    df <- prop
    p <- ggplot(df, aes(x=disease_group, y=pct, fill=subtype_final,
               text=sprintf("%s: %.1f%%", subtype_final, pct))) +
      geom_bar(stat="identity", position="fill", width=0.8, color="white", size=0.2) +
      scale_fill_manual(values=subtype_colors, name="") +
      labs(x="", y="Proportion") + theme_minimal(base_size=12) +
      theme(legend.position="bottom", axis.text.x=element_text(angle=30, hjust=1))
    ggplotly(p, tooltip="text") %>% layout(legend=list(orientation="h", y=-0.3))
  })

  output$st_table <- renderDT({
    st <- input$st_select
    cm <- unique(meta[meta$source=="Main", c("subtype_final","seurat_clusters")])
    cl <- cm$seurat_clusters[cm$subtype_final==st][1]
    if(!is.na(cl)) {
      top <- markers[markers$cluster==cl,]
      top <- top[order(-top$avg_log2FC),][1:20, c("gene","avg_log2FC","p_val_adj","pct.1")]
      top$pct.1 <- round(top$pct.1*100,1)
      top$p_val_adj <- format(top$p_val_adj, digits=2, scientific=TRUE)
      colnames(top) <- c("Gene","log2FC","p.adj","%")
      datatable(top, options=list(pageLength=10, dom="t"), rownames=FALSE)
    }
  })

  output$st_prop <- renderPlotly({
    st <- input$st_select
    df <- prop[prop$subtype_final==st,]
    p <- ggplot(df, aes(x=disease_group, y=pct, fill=subtype_final,
               text=sprintf("%s: %.1f%%", disease_group, pct))) +
      geom_bar(stat="identity", width=0.6) +
      scale_fill_manual(values=subtype_colors, guide="none") +
      labs(x="", y="%", title=paste(st)) +
      theme_minimal(base_size=12) + theme(axis.text.x=element_text(angle=30, hjust=1))
    ggplotly(p, tooltip="text")
  })

  output$st_violin <- renderPlotly({
    st <- input$st_select
    top_genes <- head(all_genes, 6)
    df <- gene_expr[gene_expr$subtype_final==st,]
    main_meta <- meta[meta$source=="Main",]
    df$disease_group <- main_meta$disease_group[match(rownames(df), main_meta$cell_id)]
    df_long <- tidyr::pivot_longer(df[,c(top_genes[1:6],"disease_group")], -disease_group, names_to="gene", values_to="expr")
    p <- ggplot(df_long, aes(x=disease_group, y=expr, fill=disease_group)) +
      geom_violin(scale="width", alpha=0.7) +
      facet_wrap(~gene, scales="free_y", nrow=2) +
      scale_fill_manual(values=disease_colors, guide="none") +
      labs(x="", y="Expression") + theme_minimal(base_size=11) +
      theme(axis.text.x=element_text(angle=30, hjust=1))
    ggplotly(p)
  })

  # Gene Lookup
  output$gene_umap <- renderPlotly({
    g <- input$gene_search; vals <- get_gene_expr(g)
    if(is.null(vals)) return(plotly_empty())
    df <- meta[meta$source=="Main",]; df$expr <- vals
    p <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=expr,
               text=paste(g, sprintf("%.2f", expr),"<br>",subtype_final))) +
      geom_point(size=0.4, alpha=0.7) +
      scale_color_viridis_c(option="magma", name=g) +
      theme_minimal(base_size=12) + labs(x="UMAP 1", y="UMAP 2")
    ggplotly(p, tooltip="text")
  })

  output$gene_violin <- renderPlotly({
    g <- input$gene_search; vals <- get_gene_expr(g)
    if(is.null(vals)) return(plotly_empty())
    df <- data.frame(expr=vals, subtype=meta$subtype_final[meta$source=="Main"])
    p <- ggplot(df, aes(x=subtype, y=expr, fill=subtype)) +
      geom_violin(scale="width", alpha=0.7) +
      scale_fill_manual(values=subtype_colors, guide="none") +
      labs(x="", y=g) + theme_minimal(base_size=12) +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    ggplotly(p)
  })

  output$gene_hm <- renderPlotly({
    g <- input$gene_search
    req(g)
    if(!g %in% rownames(mean_sd)) return(plotly_empty())
    vals <- mean_sd[g, , drop=FALSE]
    # Parse "Subtype | Disease" column names
    df <- data.frame(
      label = colnames(vals),
      expr  = as.numeric(vals[1, ])
    )
    df$subtype <- sub(" \\| .*", "", df$label)
    df$disease <- sub(".* \\| ", "", df$label)
    df <- df[!is.na(df$expr), ]
    # Set disease order: normal → disease
    disease_order <- c("Normal 1st trimester","Normal 3rd trimester",
                       "Preeclampsia","Miscarriage","Infection")
    df$disease <- factor(df$disease, levels=disease_order)
    p <- ggplot(df, aes(x=disease, y=subtype, fill=expr,
               text=sprintf("%s: %.2f", label, expr))) +
      geom_tile(color="white", size=0.3) +
      scale_fill_viridis_c(option="magma", name=paste(g,"mean")) +
      labs(x="", y="") + theme_minimal(base_size=11) +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    ggplotly(p, tooltip="text")
  })

  output$gene_stats <- renderText({
    g <- input$gene_search
    vals <- get_gene_expr(g)
    if(!is.null(vals)) {
      sprintf("Gene: %s\nMean: %.3f\nMedian: %.3f\nDetected: %.1f%%\nMax: %.3f",
        g, mean(vals), median(vals), mean(vals>0)*100, max(vals))
    } else if(g %in% rownames(full_mat)) {
      v <- as.numeric(full_mat[g, ])
      sprintf("Gene: %s\nMean: %.3f\nMedian: %.3f\nDetected: %.1f%%\nMax: %.3f",
        g, mean(v), median(v), mean(v>0)*100, max(v))
    } else {
      "Gene not found in expression data."
    }
  })

  # Disease Compare
  observe({
    col <- input$comp_group_col
    if(!is.null(col) && col %in% colnames(meta)) {
      vals <- sort(unique(meta[[col]])); vals <- vals[!is.na(vals)]
      # Exclude Preterm/Term Labor
      vals <- vals[!vals %in% c("Preterm Labor", "Term Labor")]
      updateSelectInput(session,"comp_g1", choices=vals, selected=vals[1])
      updateSelectInput(session,"comp_g2", choices=vals, selected=if(length(vals)>1) vals[2] else vals[1])
    }
  })

  get_comp <- reactive({
    list(col=input$comp_group_col, g1=input$comp_g1, g2=input$comp_g2,
         l1=input$comp_g1, l2=input$comp_g2)
  })

  output$comp_prop <- renderPlotly({
    cfg <- get_comp()
    if(is.null(cfg$col)||!cfg$col%in%colnames(meta)) return(NULL)
    df <- meta[meta[[cfg$col]] %in% c(cfg$g1, cfg$g2), ] %>%
      count(.data[[cfg$col]], subtype_final) %>% group_by(.data[[cfg$col]]) %>% mutate(pct=n/sum(n)*100)
    p <- ggplot(df, aes(x=subtype_final, y=pct, fill=.data[[cfg$col]])) +
      geom_bar(stat="identity", position="dodge", width=0.7, color="white", size=0.2) +
      scale_fill_manual(values=c("#E64B35","#4DBBD5")) +
      labs(x="", y="%", title=paste(cfg$l1, "vs", cfg$l2)) +
      theme_minimal(base_size=12) + theme(axis.text.x=element_text(angle=45, hjust=1))
    ggplotly(p)
  })

  output$comp_umap <- renderPlotly({
    cfg <- get_comp()
    if(is.null(cfg$col)||!cfg$col%in%colnames(meta)) return(NULL)
    df <- meta[meta[[cfg$col]] %in% c(cfg$g1, cfg$g2) & !is.na(meta$UMAP_1),]
    p <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=subtype_final,
               text=paste(subtype_final,"<br>",.data[[cfg$col]]))) +
      geom_point(size=0.5, alpha=0.6) +
      scale_color_manual(values=subtype_colors, name="") +
      facet_wrap(vars(.data[[cfg$col]])) + theme_minimal(base_size=12) +
      labs(x="UMAP 1", y="UMAP 2") +
      theme(legend.position="bottom")
    ggplotly(p, tooltip="text") %>% layout(legend=list(orientation="h", y=-0.2))
  })

  # Download
  output$dl_table <- renderDT(datatable(
    markers[markers$avg_log2FC>0.5, c("cluster","gene","avg_log2FC","p_val_adj")],
    options=list(pageLength=15)))
  output$dl_meta <- downloadHandler(
    filename="hofbauer_metadata.csv",
    content=function(file) write.csv(meta, file, row.names=FALSE))
  output$dl_markers <- downloadHandler(
    filename="subtype_markers.csv",
    content=function(file) write.csv(markers, file, row.names=FALSE))
  output$dl_prop <- downloadHandler(
    filename="subtype_proportions.csv",
    content=function(file) write.csv(prop, file, row.names=FALSE))
}

shinyApp(ui, server)
