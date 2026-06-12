library(shiny)
library(shinythemes)
library(ggplot2)
library(dplyr)
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

# Disease short labels
umap_meta$disease_short <- umap_meta$disease_group
umap_meta$disease_short[umap_meta$disease_group %in% c("Preterm Labor","Preterm No Labor","Term Labor")] <- "Preterm"

subtype_cols <- c(
  "Pro-inflammatory"="#C62828", "MHCII+ Antigen-presenting"="#E65100",
  "Homeostatic"="#1565C0", "PRKN+ Autophagy"="#6A1B9A",
  "Vascular remodeling"="#2E7D32", "MKI67+ Proliferating"="#455A64"
)
disease_cols <- c(
  "Normal 1st trimester"="#4DBBD5", "Normal 1st/2nd/Term"="#00A087",
  "Normal 3rd trimester / Preeclampsia"="#7E6148", "Preeclampsia"="#C62828",
  "Miscarriage / Normal"="#F39B7F", "Infection"="#DC0000", "Preterm"="#E18727"
)
trimester_cols <- c("Early"="#4575B4","Mid"="#FDAE61","Late"="#D73027")

# Prop data
prop_tbl <- as.data.frame(prop.table(table(umap_meta$subtype, umap_meta$disease_short), margin=2))
colnames(prop_tbl) <- c("Subtype","Disease","Proportion")

# ── UI ──
ui <- navbarPage(
  title = "Hofbauer Cell Atlas",
  theme = shinytheme("flatly"),

  # ═══ ATLAS ═══
  tabPanel("Atlas",
    sidebarLayout(
      sidebarPanel(
        selectInput("atlas_color","Color by:",
          choices=c("Subtype"="subtype","Trimester"="trimester","Disease"="disease_short")),
        checkboxInput("atlas_split","Split by disease group", FALSE),
        hr(),
        h5("Subtype Composition"),
        plotlyOutput("atlas_bar_small", height="220px"),
        width = 3
      ),
      mainPanel(
        plotlyOutput("atlas_umap", height="600px")
      )
    )
  ),

  # ═══ GENE LOOKUP ═══
  tabPanel("Gene Lookup",
    sidebarLayout(
      sidebarPanel(
        selectizeInput("gene_search","Gene (19,940 total):",
          choices = NULL, selected = "FOLR2",
          options = list(placeholder="Type gene name...")),
        verbatimTextOutput("gene_stats"),
        width = 3
      ),
      mainPanel(
        fluidRow(
          column(6, h4("UMAP"), plotlyOutput("gene_umap", height="420px")),
          column(6, h4("By Subtype"), plotlyOutput("gene_violin", height="420px"))
        )
      )
    )
  ),

  # ═══ DISEASE ═══
  tabPanel("Disease",
    sidebarLayout(
      sidebarPanel(
        selectInput("dis_select","Disease Group:",
          choices=c("All"="all","Miscarriage"="Miscarriage / Normal",
                    "Infection"="Infection","Preeclampsia"="Preeclampsia","Preterm"="Preterm")),
        hr(),
        h5("Cell Counts by Disease"),
        plotlyOutput("dis_bar_small", height="200px"),
        width = 3
      ),
      mainPanel(
        h4("UMAP by Subtype"),
        fluidRow(
          column(6, plotlyOutput("dis_umap_total", height="380px")),
          column(6, plotlyOutput("dis_umap_subset", height="380px"))
        ),
        hr(),
        h4("Subtype Proportions by Disease"),
        plotlyOutput("dis_prop", height="350px"),
        hr(),
        h4("Differential Expression"),
        fluidRow(
          column(4, selectInput("dis_deg","Comparison:",
            choices=setNames(deg_names, gsub("_vs_"," vs ",deg_names)))),
          column(2, numericInput("dis_topn","Top N",20,5,50,5))
        ),
        DTOutput("dis_deg_table")
      )
    )
  ),

  # ═══ SPATIAL ═══
  tabPanel("Spatial",
    sidebarLayout(
      sidebarPanel(
        h5("Stereo-seq — Placental Basal Plate"),
        p("16 mid-gestation sections at 500 nm resolution. Hofbauer cells highlighted in red."),
        selectInput("spatial_sample","Sample:",
          choices=c("Sample 001"="001","Sample 004"="004",
                    "Sample 010"="010","Sample 014"="014")),
        width = 3
      ),
      mainPanel(
        h4("Tissue Slice"),
        uiOutput("spatial_slice_img"),
        hr(),
        h4("Villus Zoom"),
        tags$img(src="Fig3b_villus_zoom.png", style="max-width:100%"),
        hr(),
        h4("Neighborhood Enrichment"),
        tags$img(src="Fig3C_neighborhood.png", style="max-width:100%")
      )
    )
  ),

  # ═══ DOWNLOAD ═══
  tabPanel("Download",
    h4("Data Tables"),
    DTOutput("dl_table"),
    br(),
    tags$a("UMAP + Metadata (CSV)", href="shiny_data/umap_meta.csv", download=NA),
    br(),
    tags$a("Subtype Proportions (CSV)", href="shiny_data/subtype_proportions.csv", download=NA),
    br(),
    lapply(deg_names, function(d) tagList(
      tags$a(sprintf("DEG: %s (CSV)", d), href=sprintf("shiny_data/deg_%s.csv",d), download=NA), tags$br()
    ))
  )
)

# ── Server ──
server <- function(input, output, session) {

  updateSelectizeInput(session, "gene_search", choices = gene_list, server = TRUE)

  # ── Atlas ──
  output$atlas_umap <- renderPlotly({
    col <- input$atlas_color
    pal <- switch(col, subtype=subtype_cols, disease_short=disease_cols, trimester=trimester_cols, NULL)
    df <- umap_meta; if(nrow(df) > 6000) df <- df[sample(nrow(df), 6000), ]
    p <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=.data[[col]],
               text=paste0(.data[[col]],"<br>",subtype))) +
      geom_point(size=0.4, alpha=0.7) +
      {if(!is.null(pal)) scale_color_manual(values=pal, name="")} +
      theme_minimal(base_size=12) + labs(x="UMAP 1", y="UMAP 2") +
      theme(legend.position="bottom")
    if(input$atlas_split) p <- p + facet_wrap(~disease_short, nrow=2)
    ggplotly(p, tooltip="text") %>% layout(legend=list(orientation="h", y=-0.2))
  })

  output$atlas_bar_small <- renderPlotly({
    cnt <- table(umap_meta$subtype)
    df <- data.frame(Subtype=factor(names(cnt), levels=names(sort(cnt))), Count=as.integer(cnt))
    p <- ggplot(df, aes(x=Count, y=Subtype, fill=Subtype)) +
      geom_bar(stat="identity") + scale_fill_manual(values=subtype_cols, guide="none") +
      theme_minimal(base_size=10) + labs(x="", y="")
    ggplotly(p, tooltip="text") %>% layout(showlegend=F, margin=list(l=130,r=5,t=5,b=5))
  })

  # ── Gene Lookup ──
  output$gene_umap <- renderPlotly({
    g <- input$gene_search
    if(!g %in% rownames(expr)) return(plotly_empty())
    df <- umap_meta; df$Expression <- as.numeric(expr[g,])
    p <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=Expression,
               text=paste0(g,": ",round(Expression,3)))) +
      geom_point(size=0.4, alpha=0.7) +
      scale_color_gradientn(colors=c("grey90","#BDD7E7","#2171B5","#08306B"), name=g) +
      theme_minimal(base_size=12) + labs(x="UMAP 1", y="UMAP 2")
    ggplotly(p, tooltip="text")
  })

  output$gene_violin <- renderPlotly({
    g <- input$gene_search
    if(!g %in% rownames(expr)) return(plotly_empty())
    df <- data.frame(Expression=as.numeric(expr[g,]), Subtype=umap_meta$subtype)
    p <- ggplot(df, aes(x=Subtype, y=Expression, fill=Subtype)) +
      geom_violin(scale="width", alpha=0.7) +
      scale_fill_manual(values=subtype_cols, guide="none") +
      labs(x="", y=g) + theme_minimal(base_size=12) +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    ggplotly(p)
  })

  output$gene_stats <- renderText({
    g <- input$gene_search
    if(!g %in% rownames(expr)) return("Gene not found.")
    vals <- as.numeric(expr[g,])
    sprintf("Gene: %s\nMean: %.3f\nMedian: %.3f\nDetected: %.1f%%\nMax: %.3f",
      g, mean(vals), median(vals), mean(vals>0)*100, max(vals))
  })

  # ── Disease ──
  output$dis_bar_small <- renderPlotly({
    dis_order <- c("Normal 1st trimester","Miscarriage / Normal","Infection",
                   "Normal 3rd trimester / Preeclampsia","Preeclampsia","Preterm")
    tbl <- table(umap_meta$disease_short)
    tbl <- tbl[intersect(dis_order, names(tbl))]
    df <- data.frame(Disease=factor(names(tbl), levels=names(tbl)), Count=as.integer(tbl))
    p <- ggplot(df, aes(x=Count, y=Disease, fill=Disease)) +
      geom_bar(stat="identity") + scale_fill_manual(values=disease_cols, guide="none") +
      theme_minimal(base_size=10) + labs(x="", y="")
    ggplotly(p, tooltip="text") %>% layout(showlegend=F, margin=list(l=120,r=5,t=5,b=5))
  })

  output$dis_umap_total <- renderPlotly({
    df <- umap_meta; if(nrow(df)>5000) df <- df[sample(nrow(df),5000),]
    p <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=subtype, text=subtype)) +
      geom_point(size=0.3, alpha=0.7) +
      scale_color_manual(values=subtype_cols, name="") +
      theme_minimal(base_size=11) + labs(x="UMAP 1", y="UMAP 2", title="All Cells") +
      theme(legend.position="bottom", legend.text=element_text(size=7))
    ggplotly(p, tooltip="text") %>% layout(legend=list(orientation="h", y=-0.2))
  })

  output$dis_umap_subset <- renderPlotly({
    req(input$dis_select)
    if(input$dis_select=="all") df <- umap_meta
    else df <- umap_meta[umap_meta$disease_short==input$dis_select,]
    if(nrow(df)==0) return(plotly_empty())
    if(nrow(df)>3000) df <- df[sample(nrow(df),3000),]
    p <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=subtype)) +
      geom_point(size=0.5, alpha=0.7) +
      scale_color_manual(values=subtype_cols, guide="none") +
      theme_minimal(base_size=11) + labs(x="UMAP 1", y="UMAP 2", title=input$dis_select)
    ggplotly(p)
  })

  output$dis_prop <- renderPlotly({
    p <- ggplot(prop_tbl, aes(x=Disease, y=Proportion, fill=Subtype,
               text=sprintf("%s: %.1f%%", Subtype, Proportion*100))) +
      geom_bar(stat="identity", color="white", size=0.2) +
      scale_fill_manual(values=subtype_cols, name="") +
      labs(x="", y="Proportion") + theme_minimal(base_size=12) +
      theme(legend.position="bottom", axis.text.x=element_text(angle=30, hjust=1))
    ggplotly(p, tooltip="text") %>% layout(legend=list(orientation="h", y=-0.3))
  })

  output$dis_deg_table <- renderDT({
    req(input$dis_deg)
    deg <- read.csv(deg_list[[input$dis_deg]])
    deg <- deg[!is.na(deg$p_val_adj), ]; deg <- deg[order(deg$p_val_adj), ]
    top <- head(deg, input$dis_topn)
    top <- top[, c("gene","avg_log2FC","p_val_adj","pct.1","pct.2")]
    top$avg_log2FC <- round(top$avg_log2FC, 3)
    top$p_val_adj  <- format.pval(top$p_val_adj, digits=2, eps=1e-300)
    top$pct.1 <- round(top$pct.1*100, 1); top$pct.2 <- round(top$pct.2*100, 1)
    colnames(top) <- c("Gene","log2FC","padj","% Disease","% Control")
    datatable(top, rownames=FALSE, options=list(pageLength=20, dom='tip', searching=FALSE),
      class='compact stripe hover') %>%
      formatStyle("log2FC", color=styleInterval(0, c("#4575B4","#D73027")))
  })

  # Download table
  output$dl_table <- renderDT({
    datatable(umap_meta[,c("subtype","disease_group","UMAP_1","UMAP_2")],
      options=list(pageLength=15))
  })

  # ── Spatial ──
  output$spatial_slice_img <- renderUI({
    req(input$spatial_sample)
    tags$div(
      tags$p(tags$b(sprintf("Stereo-seq — Sample %s", input$spatial_sample)),
        tags$br(), "Hofbauer cells in red. Bin50 resolution."),
      tags$img(src=sprintf("Fig3a_spatial_%s.png", input$spatial_sample),
               style="max-width:100%;height:auto;border:1px solid #ddd;border-radius:4px;"))
  })
}

shinyApp(ui, server)
