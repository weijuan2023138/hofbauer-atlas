library(shiny)
library(shinythemes)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)
library(Matrix)

DATA_DIR <- "shiny_data"
cat("Loading data...\n")
umap_meta <- read.csv(file.path(DATA_DIR, "umap_meta.csv"))
expr      <- readRDS(file.path(DATA_DIR, "expr_full.rds"))
gene_list <- sort(rownames(expr))
deg_files  <- list.files(DATA_DIR, pattern="^deg_", full.names=TRUE)
deg_names  <- gsub("deg_|\\.csv","", basename(deg_files))
deg_list   <- setNames(deg_files, deg_names)

umap_meta$disease_short <- umap_meta$disease_group
umap_meta$disease_short[umap_meta$disease_group == "Normal 1st trimester"] <- "Normal 1st"
umap_meta$disease_short[umap_meta$disease_group == "Normal 1st/2nd/Term"] <- "Normal 1st/2nd"
umap_meta$disease_short[umap_meta$disease_group == "Normal 3rd trimester / Preeclampsia"] <- "Normal 3rd"
umap_meta$disease_short[umap_meta$disease_group == "Miscarriage / Normal"] <- "Miscarriage"
umap_meta$disease_short[umap_meta$disease_group == "Preeclampsia"] <- "PE"
umap_meta$disease_short[umap_meta$disease_group == "Preeclampsia / Control"] <- "PE"
umap_meta$disease_short[umap_meta$disease_group == "Severe Preeclampsia"] <- "PE"
umap_meta$disease_short[umap_meta$disease_group == "Infection"] <- "Infection"
umap_meta$disease_short[umap_meta$disease_group == "Preterm Labor / Term Labor"] <- "Preterm"
dis_levels <- c("Normal 1st","Miscarriage","Infection","Normal 3rd","PE","Preterm")
umap_meta$disease_short <- factor(umap_meta$disease_short, levels=dis_levels)

subtype_cols <- c(
  "Pro-inflammatory"="#C62828", "MHCII+ Antigen-presenting"="#E65100",
  "Homeostatic"="#1565C0", "PRKN+ Autophagy"="#6A1B9A",
  "Vascular remodeling"="#2E7D32", "MKI67+ Proliferating"="#455A64"
)
disease_cols <- c(
  "Normal 1st"="#4575B4","Miscarriage"="#D73027","Infection"="#FDB462",
  "Normal 3rd"="#66C2A5","PE"="#FC8D59","Preterm"="#E41A1C"
)
trimester_cols <- c("Early"="#4575B4","Mid"="#FDAE61","Late"="#D73027")

prop_tbl <- as.data.frame(prop.table(table(umap_meta$subtype, umap_meta$disease_short), margin=2))
colnames(prop_tbl) <- c("Subtype","Disease","Proportion")

# ── UI ──
ui <- navbarPage(
  title = "Hofbauer Cell Atlas",
  theme = shinytheme("flatly"),

  # ═══ OVERVIEW ═══
  tabPanel("Overview",
    fluidRow(
      column(4,
        h4("Hofbauer Cell Atlas"),
        p(sprintf("%d Hofbauer cells from 11 datasets", nrow(umap_meta))),
        p("6 subtypes across normal pregnancy, miscarriage, infection, preeclampsia, and preterm birth"),
        p("Integrated scRNA-seq + snATAC-seq + Stereo-seq spatial transcriptomics"),
        hr(),
        selectInput("ov_color", "Color by:",
          choices=c("Subtype"="subtype","Disease"="disease_short","Trimester"="trimester","Dataset"="dataset")),
        checkboxInput("ov_split", "Split by disease group", FALSE)
      ),
      column(8, plotlyOutput("ov_umap", height="550px"))
    ),
    hr(),
    h4("Subtype Composition by Disease Group"),
    plotlyOutput("ov_prop", height="400px")
  ),

  # ═══ GENE LOOKUP ═══
  tabPanel("Gene Lookup",
    sidebarLayout(
      sidebarPanel(
        selectizeInput("gene_search", "Gene (18,492 total):",
          choices = gene_list, selected = "FOLR2",
          options = list(placeholder = "Type any gene name...", create = TRUE, maxOptions = 5000)),
        actionButton("gene_go", "Search", class = "btn-primary", width = "100%"),
        tags$hr(),
        verbatimTextOutput("gene_stats"),
        width = 3
      ),
      mainPanel(
        fluidRow(
          column(6, h4("UMAP"), plotlyOutput("gene_umap", height="420px")),
          column(6, h4("By Subtype"), plotlyOutput("gene_violin", height="420px"))
        ),
        h4("Mean Expression by Disease Group"),
        plotlyOutput("gene_hm", height="420px"),
        hr(),
        fluidRow(
          column(6, h4("By Disease Group"), plotlyOutput("gene_dis_violin", height="420px")),
          column(6, h4("By Trimester"), plotlyOutput("gene_tri_violin", height="420px"))
        )
      )
    )
  ),

  # ═══ DISEASE ═══
  tabPanel("Disease Comparison",
    sidebarLayout(
      sidebarPanel(
        selectInput("comp_group", "Group by:",
          choices = c("Subtype"="subtype","Disease"="disease_short","Trimester"="trimester")),
        selectInput("comp_g1", "Group A:", choices = NULL),
        selectInput("comp_g2", "Group B:", choices = NULL),
        width = 3
      ),
      mainPanel(
        h4("UMAP Comparison"),
        plotlyOutput("comp_umap", height="500px"),
        h4("Subtype Proportions"),
        plotlyOutput("comp_prop", height="400px")
      )
    )
  ),

  # ═══ SPATIAL ═══
  tabPanel("Spatial (Stereo-seq)",
    sidebarLayout(
      sidebarPanel(
        h5("Stereo-seq — Placental Basal Plate"),
        p("16 mid-gestation sections at 500 nm resolution."),
        p("Hofbauer cells highlighted in red. Bin50 (~50 μm) resolution."),
        hr(),
        selectInput("spatial_sample", "Sample:",
          choices = c("Sample 001"="001","Sample 004"="004",
                      "Sample 010"="010","Sample 014"="014")),
        width = 3
      ),
      mainPanel(
        h4("Tissue Slice"),
        uiOutput("spatial_slice_img"),
        hr(),
        h4("Villus Zoom"),
        p("Hofbauer cells (red) near fetal capillaries, surrounded by fibroblasts (orange)."),
        tags$img(src="Fig3b_villus_zoom.png", style="max-width:100%"),
        hr(),
        h4("Neighborhood Enrichment"),
        p("Cell-type enrichment within 50 μm of Hofbauer cells."),
        tags$img(src="Fig3C_neighborhood.png", style="max-width:100%")
      )
    )
  ),

  # ═══ DOWNLOAD ═══
  tabPanel("Download",
    sidebarLayout(
      sidebarPanel(
        h5("Filter Cells"),
        checkboxGroupInput("dl_disease", "Disease Group:",
          choices = c("Normal 1st","Miscarriage","Infection","Normal 3rd","PE","Preterm"),
          selected = c("Normal 1st","Miscarriage","Infection","Normal 3rd","PE","Preterm")),
        checkboxGroupInput("dl_subtype", "Subtype:",
          choices = names(subtype_cols), selected = names(subtype_cols)),
        hr(),
        h5("Select Genes"),
        radioButtons("dl_gene_mode", "Genes:",
          choices = c("All 18,492 genes"="all","Top 5000 variable genes"="var",
                      "Custom gene list"="custom")),
        conditionalPanel("input.dl_gene_mode=='custom'",
          textAreaInput("dl_custom_genes", "Enter genes (one per line or comma-separated):",
            rows = 6, placeholder = "FOLR2\nCD163\nSPP1\n...")),
        hr(),
        downloadButton("dl_expr", "Download CSV", class="btn-primary", width="100%"),
        width = 3
      ),
      mainPanel(
        h4("Expression Matrix Download"),
        p("Filter by disease group, subtype, and genes. Downloads log-normalized expression values."),
        verbatimTextOutput("dl_preview")
      )
    )
  )
)

# ── Server ──
server <- function(input, output, session) {

  updateSelectizeInput(session, "gene_search", choices = gene_list, server = TRUE)

  current_gene <- reactiveVal("FOLR2")
  observeEvent(input$gene_go, { current_gene(input$gene_search) })

  # ── Overview ──
  output$ov_umap <- renderPlotly({
    col <- input$ov_color
    pal <- switch(col, subtype=subtype_cols, disease_short=disease_cols,
                  trimester=trimester_cols, dataset=NULL)
    df <- umap_meta; if(nrow(df) > 6000) df <- df[sample(nrow(df), 6000), ]
    p <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=.data[[col]],
               text=paste0(.data[[col]],"<br>",subtype))) +
      geom_point(size=0.4, alpha=0.7) +
      {if(!is.null(pal)) scale_color_manual(values=pal, name="")} +
      theme_minimal(base_size=12) + labs(x="UMAP 1", y="UMAP 2") +
      theme(legend.position="bottom")
    if(input$ov_split) p <- p + facet_wrap(~disease_short, nrow=2)
    ggplotly(p, tooltip="text") %>% layout(legend=list(orientation="h", y=-0.2))
  })

  output$ov_prop <- renderPlotly({
    p <- ggplot(prop_tbl, aes(x=Disease, y=Proportion, fill=Subtype,
               text=sprintf("%s: %.1f%%", Subtype, Proportion*100))) +
      geom_bar(stat="identity", position="fill", width=0.8, color="white", size=0.2) +
      scale_fill_manual(values=subtype_cols, name="") +
      labs(x="", y="Proportion") + theme_minimal(base_size=12) +
      theme(legend.position="bottom", axis.text.x=element_text(angle=30, hjust=1))
    ggplotly(p, tooltip="text") %>% layout(legend=list(orientation="h", y=-0.3))
  })

  # ── Gene Lookup ──
  output$gene_umap <- renderPlotly({
    g <- current_gene(); if(!g %in% rownames(expr)) return(plotly_empty())
    df <- umap_meta; df$Expression <- as.numeric(expr[g,])
    p <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=Expression,
               text=paste0(g,": ",round(Expression,3)))) +
      geom_point(size=0.4, alpha=0.7) +
      scale_color_gradientn(colors=c("grey90","#BDD7E7","#2171B5","#08306B"), name=g) +
      theme_minimal(base_size=12) + labs(x="UMAP 1", y="UMAP 2")
    ggplotly(p, tooltip="text")
  })

  output$gene_violin <- renderPlotly({
    g <- current_gene(); if(!g %in% rownames(expr)) return(plotly_empty())
    df <- data.frame(Expression=as.numeric(expr[g,]), Subtype=umap_meta$subtype)
    p <- ggplot(df, aes(x=Subtype, y=Expression, fill=Subtype)) +
      geom_violin(scale="width", alpha=0.7) +
      scale_fill_manual(values=subtype_cols, guide="none") +
      labs(x="", y=g) + theme_minimal(base_size=12) +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    ggplotly(p)
  })

  output$gene_hm <- renderPlotly({
    g <- current_gene(); if(!g %in% rownames(expr)) return(plotly_empty())
    df <- umap_meta; df$Expression <- as.numeric(expr[g,])
    keeps <- c("Normal 1st","Miscarriage","Infection","Normal 3rd","PE","Preterm")
    df <- df[df$disease_short %in% keeps, ]
    means <- df %>% group_by(subtype, disease_short) %>%
      summarise(Mean=mean(Expression), .groups="drop")
    p <- ggplot(means, aes(x=disease_short, y=subtype, fill=Mean,
               text=sprintf("%s: %.2f", disease_short, Mean))) +
      geom_tile(color="white", size=0.3) +
      scale_fill_gradientn(colors=c("grey90","#FDAE61","#D73027"), name=g) +
      labs(x="", y="") + theme_minimal(base_size=11) +
      theme(axis.text.x=element_text(angle=30, hjust=1))
    ggplotly(p, tooltip="text")
  })

  output$gene_stats <- renderText({
    g <- current_gene()
    if(!g %in% rownames(expr)) return("Gene not found.")
    vals <- as.numeric(expr[g,])
    sprintf("Gene: %s\nMean: %.3f\nMedian: %.3f\nDetected: %.1f%%\nMax: %.3f",
      g, mean(vals), median(vals), mean(vals>0)*100, max(vals))
  })

  output$gene_dis_violin <- renderPlotly({
    g <- current_gene(); if(!g %in% rownames(expr)) return(plotly_empty())
    dis_order <- c("Normal 1st","Miscarriage","Infection","Normal 3rd","PE","Preterm")
    df <- data.frame(Expression=as.numeric(expr[g,]), Disease=umap_meta$disease_short)
    df <- df[df$Disease %in% dis_order, ]
    df$Disease <- factor(df$Disease, levels=intersect(dis_order, unique(df$Disease)))
    p <- ggplot(df, aes(x=Disease, y=Expression, fill=Disease)) +
      geom_violin(scale="width", alpha=0.7) +
      scale_fill_manual(values=disease_cols, guide="none") +
      labs(x="", y=g) + theme_minimal(base_size=12) +
      theme(axis.text.x=element_text(angle=30, hjust=1))
    ggplotly(p)
  })

  output$gene_tri_violin <- renderPlotly({
    g <- current_gene(); if(!g %in% rownames(expr)) return(plotly_empty())
    df <- data.frame(Expression=as.numeric(expr[g,]), Trimester=umap_meta$trimester)
    df <- df[df$Trimester %in% c("Early","Mid","Late"), ]
    df$Trimester <- factor(df$Trimester, levels=c("Early","Mid","Late"))
    p <- ggplot(df, aes(x=Trimester, y=Expression, fill=Trimester)) +
      geom_violin(scale="width", alpha=0.7) +
      scale_fill_manual(values=trimester_cols, guide="none") +
      labs(x="", y=g) + theme_minimal(base_size=12)
    ggplotly(p)
  })

  # ── Disease Comparison ──
  observe({
    grp <- input$comp_group
    vals <- sort(unique(umap_meta[[grp]])); vals <- vals[!is.na(vals)]
    updateSelectInput(session, "comp_g1", choices=vals, selected=vals[1])
    updateSelectInput(session, "comp_g2", choices=vals,
      selected=if(length(vals)>1) vals[2] else vals[1])
  })

  get_comp <- reactive({
    list(grp=input$comp_group, g1=input$comp_g1, g2=input$comp_g2)
  })

  output$comp_umap <- renderPlotly({
    cfg <- get_comp()
    df <- umap_meta[umap_meta[[cfg$grp]] %in% c(cfg$g1, cfg$g2), ]
    if(nrow(df) > 4000) df <- df[sample(nrow(df), 4000), ]
    p <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=subtype,
               text=paste0(subtype,"<br>",.data[[cfg$grp]]))) +
      geom_point(size=0.5, alpha=0.6) +
      scale_color_manual(values=subtype_cols, name="") +
      facet_wrap(vars(.data[[cfg$grp]])) +
      theme_minimal(base_size=12) + labs(x="UMAP 1", y="UMAP 2") +
      theme(legend.position="bottom")
    ggplotly(p, tooltip="text") %>% layout(legend=list(orientation="h", y=-0.2))
  })

  output$comp_prop <- renderPlotly({
    cfg <- get_comp()
    df <- umap_meta[umap_meta[[cfg$grp]] %in% c(cfg$g1, cfg$g2), ] %>%
      count(.data[[cfg$grp]], subtype) %>%
      group_by(.data[[cfg$grp]]) %>% mutate(pct=n/sum(n)*100) %>% ungroup()
    p <- ggplot(df, aes(x=subtype, y=pct, fill=.data[[cfg$grp]])) +
      geom_bar(stat="identity", position="dodge", width=0.7, color="white", size=0.2) +
      scale_fill_manual(values=c("#E64B35","#4DBBD5")) +
      labs(x="", y="%", title=paste(cfg$g1, "vs", cfg$g2)) +
      theme_minimal(base_size=12) + theme(axis.text.x=element_text(angle=45, hjust=1))
    ggplotly(p)
  })

  output$dis_deg_table <- renderDT({
    req(input$dis_deg)
    deg <- read.csv(deg_list[[input$dis_deg]]); deg <- deg[!is.na(deg$p_val_adj), ]
    deg <- deg[order(deg$p_val_adj), ]; top <- head(deg, input$dis_topn)
    top <- top[, c("gene","avg_log2FC","p_val_adj","pct.1","pct.2")]
    top$avg_log2FC <- round(top$avg_log2FC, 3)
    top$p_val_adj  <- format.pval(top$p_val_adj, digits=2, eps=1e-300)
    top$pct.1 <- round(top$pct.1*100, 1); top$pct.2 <- round(top$pct.2*100, 1)
    colnames(top) <- c("Gene","log2FC","padj","% Disease","% Control")
    datatable(top, rownames=FALSE, options=list(pageLength=20, dom='tip', searching=FALSE),
      class='compact stripe hover') %>%
      formatStyle("log2FC", color=styleInterval(0, c("#4575B4","#D73027")))
  })

  # ── Spatial ──
  output$spatial_slice_img <- renderUI({
    req(input$spatial_sample)
    tags$div(
      tags$img(src=sprintf("Fig3a_spatial_%s.png", input$spatial_sample),
               style="max-width:100%;height:auto;border:1px solid #ddd;border-radius:4px;"))
  })

  # ── Download ──
  output$dl_preview <- renderText({
    cells_keep <- umap_meta$disease_short %in% input$dl_disease &
                  umap_meta$subtype %in% input$dl_subtype
    n_cells <- sum(cells_keep)
    
    if(input$dl_gene_mode == "all") {
      n_genes <- nrow(expr)
    } else if(input$dl_gene_mode == "var") {
      n_genes <- min(5000, nrow(expr))
    } else {
      genes <- strsplit(gsub("[,\\n]+", ",", input$dl_custom_genes), ",")[[1]]
      genes <- trimws(genes); genes <- genes[genes != ""]
      n_genes <- sum(genes %in% rownames(expr))
    }
    sprintf("Selection: %d cells × %d genes\nEstimated size: ~%.1f MB",
      n_cells, n_genes, n_cells * n_genes * 8 / 1e6)
  })

  output$dl_expr <- downloadHandler(
    filename = function() {
      sprintf("hofbauer_expr_%dcells_%s.csv.gz", 
        sum(umap_meta$disease_short %in% input$dl_disease & umap_meta$subtype %in% input$dl_subtype),
        format(Sys.time(), "%Y%m%d"))
    },
    content = function(file) {
      cells_keep <- umap_meta$disease_short %in% input$dl_disease &
                    umap_meta$subtype %in% input$dl_subtype
      
      if(input$dl_gene_mode == "all") {
        gene_idx <- 1:nrow(expr)
      } else if(input$dl_gene_mode == "var") {
        gene_idx <- 1:min(5000, nrow(expr))
      } else {
        genes <- strsplit(gsub("[,\\n]+", ",", input$dl_custom_genes), ",")[[1]]
        genes <- trimws(genes); genes <- genes[genes != ""]
        gene_idx <- which(rownames(expr) %in% genes)
      }
      
      mat <- as.matrix(expr[gene_idx, cells_keep, drop=FALSE])
      write.csv(mat, gzfile(file))
    }
  )
}

shinyApp(ui, server)
