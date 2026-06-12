library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(plotly)
library(DT)
library(Matrix)

# ── Load data ──
DATA_DIR <- "shiny_data"
umap_meta <- read.csv(file.path(DATA_DIR, "umap_meta.csv"))
expr      <- readRDS(file.path(DATA_DIR, "expr_full.rds"))
tf_expr   <- readRDS(file.path(DATA_DIR, "tf_expr.rds"))
mod_scores <- read.csv(file.path(DATA_DIR, "module_scores.csv"))
prop_df    <- read.csv(file.path(DATA_DIR, "subtype_proportions.csv"))
colnames(prop_df) <- c("Subtype","Disease","Proportion")
gene_list  <- sort(rownames(expr))
atac_motif <- read.csv(file.path(DATA_DIR, "atac_motif_enrichment.csv"))
atac_peaks <- read.csv(file.path(DATA_DIR, "atac_differential_peaks.csv"))
deg_files  <- list.files(DATA_DIR, pattern="^deg_", full.names=TRUE)
deg_names  <- gsub("deg_|\\.csv","", basename(deg_files))
deg_list   <- setNames(deg_files, deg_names)

# ── Summary stats ──
n_cells   <- nrow(umap_meta)
n_subtypes <- length(unique(umap_meta$subtype))
n_datasets <- length(unique(umap_meta$dataset))
subtype_counts <- sort(table(umap_meta$subtype), decreasing=TRUE)

# ── Palettes ──
subtype_cols <- c(
  "Pro-inflammatory"="#C62828", "MHCII+ Antigen-presenting"="#E65100",
  "Homeostatic"="#1565C0", "PRKN+ Autophagy"="#6A1B9A",
  "Vascular remodeling"="#2E7D32", "MKI67+ Proliferating"="#455A64"
)
disease_cols <- c(
  "Normal 1st trimester"="#4DBBD5", "Normal 1st/2nd/Term"="#00A087",
  "Normal 3rd trimester / Preeclampsia"="#7E6148", "Preeclampsia"="#C62828",
  "Preterm Labor"="#E18727", "Preterm No Labor"="#BCAAA4", "Term Labor"="#3C5488",
  "Miscarriage / Normal"="#F39B7F", "Infection"="#DC0000"
)
trimester_cols <- c("Early"="#4575B4","Mid"="#FDAE61","Late"="#D73027")

# ── Helper: individual UMAP panel ──
make_umap_panel <- function(df, col, pal, title, n=3000) {
  if(nrow(df) > n) df <- df[sample(nrow(df), n), ]
  col_vals <- df[[col]]
  plot_ly(df, x=~UMAP_1, y=~UMAP_2, color=col_vals, colors=pal,
    type='scatter', mode='markers', marker=list(size=2.5, opacity=0.7),
    text=~paste(col,":", col_vals), hoverinfo='text') %>%
    layout(title=list(text=title, font=list(size=13)),
           xaxis=list(title="", showgrid=F, zeroline=F, showticklabels=F),
           yaxis=list(title="", showgrid=F, zeroline=F, showticklabels=F),
           showlegend=FALSE, margin=list(l=5,r=5,t=35,b=5)) %>%
    config(displayModeBar=FALSE)
}

# ── UI ──
ui <- page_navbar(
  title = "Hofbauer Cell Atlas",
  theme = bs_theme(version=5, bootswatch="flatly", primary="#1565C0"),

  # ═══ ATLAS ═══
  nav_panel("Atlas",
    layout_column_wrap(width=1/2, heights_equal="row",
      card(full_screen=FALSE,
        card_header("Dataset Summary"),
        layout_column_wrap(width=1/3,
          value_box("Total Cells", n_cells, theme="primary"),
          value_box("Subtypes", n_subtypes, theme="success"),
          value_box("Datasets", n_datasets, theme="warning")),
        tags$hr(),
        tags$p("Integrated single-cell atlas of human placental Hofbauer cells 
                (fetal tissue-resident macrophages) spanning GW4.5–38. 
                Complemented by snATAC-seq and Stereo-seq spatial transcriptomics.")),
      card(full_screen=FALSE,
        card_header("Subtype Composition"),
        plotlyOutput("atlas_bar", height="250px"))
    ),
    tags$hr(),
    card(full_screen=TRUE,
      card_header("UMAP — Subtype"),
      plotlyOutput("umap_subtype", height="300px")),
    card(full_screen=TRUE,
      card_header("UMAP — Trimester"),
      layout_column_wrap(width=1/4, heights_equal="row",
        card(plotlyOutput("umap_tri_all", height="260px")),
        card(plotlyOutput("umap_tri_early", height="260px")),
        card(plotlyOutput("umap_tri_mid", height="260px")),
        card(plotlyOutput("umap_tri_late", height="260px"))
      )),
    card(full_screen=TRUE,
      card_header("UMAP — Disease Group"),
      layout_column_wrap(width=1/5, heights_equal="row",
        card(plotlyOutput("umap_dis_all", height="230px")),
        card(plotlyOutput("umap_dis_normal1", height="230px")),
        card(plotlyOutput("umap_dis_pe", height="230px")),
        card(plotlyOutput("umap_dis_infection", height="230px")),
        card(plotlyOutput("umap_dis_miscarriage", height="230px"))
      ))),

  # ═══ GENE EXPRESSION ═══
  nav_panel("Gene Expression",
    layout_sidebar(
      sidebar=sidebar(
        selectizeInput("gene_search","Gene (19,940 available)",
          choices=NULL, selected="FOLR2",
          options=list(placeholder='Type gene name...',maxOptions=2000)),
        actionButton("gene_go","Search", class="btn-primary", width="100%"),
        tags$hr(),
        selectInput("gene_group","Group by",
          choices=c("Subtype"="subtype","Trimester"="trimester")), width=250),
      layout_column_wrap(width=1/2,
        card(full_screen=TRUE, card_header("UMAP"), plotlyOutput("gene_umap",height="500px")),
        card(full_screen=TRUE, card_header("Expression Distribution"), plotlyOutput("gene_violin",height="500px"))
      ))),

  # ═══ DISEASE ═══
  nav_panel("Disease",
    layout_sidebar(
      sidebar=sidebar(
        selectInput("disease_comp","Comparison",
          choices=setNames(deg_names, gsub("_vs_"," vs ",deg_names))),
        radioButtons("disease_view","View",
          choices=c("Volcano Plot"="volcano","Subtype Proportions"="prop")), width=250),
      card(full_screen=TRUE, plotlyOutput("disease_plot",height="620px"))
    )),

  # ═══ ATAC-seq ═══
  nav_panel("ATAC-seq",
    navset_card_underline(
      nav_panel("Motif Enrichment",
        plotlyOutput("atac_motif_plot",height="580px")),
      nav_panel("Differential Accessibility",
        plotlyOutput("atac_volcano_plot",height="620px")),
      nav_panel("Coverage Tracks",
        tags$div(style="padding:10px;",
          tags$img(src="Fig6b_ATAC_tracks.png",style="max-width:100%;height:auto"))))),

  # ═══ SPATIAL ═══
  nav_panel("Spatial",
    navset_card_underline(
      nav_panel("Tissue Slices",
        layout_sidebar(
          sidebar=sidebar(selectInput("spatial_sample","Sample",
            choices=c("Sample 001"="001","Sample 004"="004",
                      "Sample 010"="010","Sample 014"="014")), width=200, open=FALSE),
          uiOutput("spatial_slice_img"))),
      nav_panel("Villus Zoom",
        tags$div(style="padding:10px;",
          tags$img(src="Fig3b_villus_zoom.png",style="max-width:100%;height:auto"))),
      nav_panel("Neighborhood Enrichment",
        tags$div(style="padding:10px;",
          tags$img(src="Fig3C_neighborhood.png",style="max-width:100%;height:auto"))))),

  # ═══ DOWNLOADS ═══
  nav_panel("Downloads",
    layout_column_wrap(width=1/2,
      card(card_header("Metadata & Scores"),
        tags$ul(
          tags$li(tags$a("UMAP coordinates + metadata (CSV)", href="shiny_data/umap_meta.csv", download=NA)),
          tags$li(tags$a("Subtype proportions (CSV)", href="shiny_data/subtype_proportions.csv", download=NA)),
          tags$li(tags$a("Module scores (CSV)", href="shiny_data/module_scores.csv", download=NA)))),
      card(card_header("ATAC-seq"),
        tags$ul(
          tags$li(tags$a("Motif enrichment (CSV)", href="shiny_data/atac_motif_enrichment.csv", download=NA)),
          tags$li(tags$a("Differential peaks (CSV)", href="shiny_data/atac_differential_peaks.csv", download=NA)))),
      card(card_header("Differential Expression"),
        tags$ul(lapply(deg_names, function(d)
          tags$li(tags$a(sprintf("%s (CSV)", d), href=sprintf("shiny_data/deg_%s.csv",d), download=NA))))),
      card(card_header("Code & Data"),
        tags$p("Analysis code available on GitHub. Full Seurat object (229 MB) 
                available upon request from corresponding author."))))
)

# ── Server ──
server <- function(input, output, session) {
  updateSelectizeInput(session,"gene_search",choices=gene_list,server=TRUE)
  current_gene <- reactiveVal("FOLR2")
  observeEvent(input$gene_go, { current_gene(input$gene_search) })

  # ── Atlas bar ──
  output$atlas_bar <- renderPlotly({
    df <- data.frame(Subtype=names(subtype_counts), Count=as.integer(subtype_counts))
    df$Subtype <- factor(df$Subtype, levels=names(subtype_counts))
    plot_ly(df, x=~Count, y=~Subtype, type='bar', orientation='h',
      marker=list(color=subtype_cols[as.character(df$Subtype)]),
      text=~Subtype, hoverinfo='text') %>%
      layout(xaxis=list(title="Number of cells"), yaxis=list(title="", categoryorder="total ascending"),
             margin=list(l=180), showlegend=F)
  })

  # ── UMAP: Subtype ──
  output$umap_subtype <- renderPlotly({
    df <- umap_meta[sample(nrow(umap_meta), 5000), ]
    plot_ly(df, x=~UMAP_1, y=~UMAP_2, color=~subtype, colors=subtype_cols,
      type='scatter', mode='markers', marker=list(size=3, opacity=0.7),
      text=~paste("Subtype:",subtype), hoverinfo='text') %>%
      layout(xaxis=list(title="UMAP 1", showgrid=F, zeroline=F),
             yaxis=list(title="UMAP 2", showgrid=F, zeroline=F),
             legend=list(orientation='v', y=0.5, font=list(size=10)))
  })

  # ── UMAP: Trimester (4 panels) ──
  output$umap_tri_all <- renderPlotly({
    make_umap_panel(umap_meta, "trimester", trimester_cols, "All", 4000)
  })
  for(tr in c("Early","Mid","Late")) {
    local({
      t <- tr
      output[[paste0("umap_tri_",tolower(t))]] <- renderPlotly({
        df <- umap_meta[umap_meta$trimester == t, ]
        if(nrow(df)==0) return(plotly_empty())
        make_umap_panel(df, "trimester", trimester_cols, t, 2000)
      })
    })
  }

  # ── UMAP: Disease (5 panels) ──
  dis_groups <- c("Normal 1st trimester","Preeclampsia","Infection","Miscarriage / Normal")
  dis_labels <- c("normal1","pe","infection","miscarriage")
  output$umap_dis_all <- renderPlotly({
    make_umap_panel(umap_meta, "disease_group", disease_cols, "All", 4000)
  })
  for(i in seq_along(dis_groups)) {
    local({
      dg <- dis_groups[i]; lb <- dis_labels[i]
      output[[paste0("umap_dis_",lb)]] <- renderPlotly({
        df <- umap_meta[umap_meta$disease_group == dg, ]
        if(nrow(df)==0) return(plotly_empty())
        make_umap_panel(df, "disease_group", disease_cols, dg, 2000)
      })
    })
  }

  # ── Gene ──
  output$gene_umap <- renderPlotly({
    gene <- current_gene()
    if(!gene %in% rownames(expr)) return(NULL)
    df <- umap_meta; df$Expression <- as.numeric(expr[gene,])
    plot_ly(df, x=~UMAP_1, y=~UMAP_2, color=~Expression,
      colors=colorRamp(c("grey90","#BDD7E7","#2171B5","#08306B")),
      type='scatter', mode='markers', marker=list(size=3, opacity=0.7),
      text=~paste(gene,":",round(Expression,3)), hoverinfo='text') %>%
      layout(title=list(text=gene,font=list(size=16)),
             xaxis=list(title="UMAP 1"), yaxis=list(title="UMAP 2"))
  })
  output$gene_violin <- renderPlotly({
    gene <- current_gene()
    if(!gene %in% rownames(expr)) return(NULL)
    grp <- input$gene_group
    df <- data.frame(Group=umap_meta[[grp]], Expression=as.numeric(expr[gene,]))
    pal <- if(grp=="subtype") subtype_cols else trimester_cols
    plot_ly(df, x=~Group, y=~Expression, color=~Group, colors=pal,
      type='violin', box=list(visible=T, width=0.1),
      points='all', jitter=0.3, pointpos=-0.5, marker=list(size=1, opacity=0.3)) %>%
      layout(title=list(text=gene,font=list(size=16)),
             xaxis=list(title=""), yaxis=list(title="Log-normalized expression"), showlegend=F)
  })

  # ── Disease ──
  output$disease_plot <- renderPlotly({
    req(input$disease_comp, input$disease_view)
    if(input$disease_view=="volcano"){
      deg <- read.csv(deg_list[[input$disease_comp]]); deg <- deg[!is.na(deg$p_val_adj),]
      deg$logP <- -log10(deg$p_val_adj+1e-300)
      deg$sig <- "NS"
      deg$sig[deg$p_val_adj<0.05 & abs(deg$avg_log2FC)>0.5 & deg$avg_log2FC>0] <- "Up in Disease"
      deg$sig[deg$p_val_adj<0.05 & abs(deg$avg_log2FC)>0.5 & deg$avg_log2FC<0] <- "Up in Control"
      top <- rbind(head(deg[deg$sig=="Up in Disease",][order(-deg[deg$sig=="Up in Disease",]$avg_log2FC),],10),
                   head(deg[deg$sig=="Up in Control",][order(deg[deg$sig=="Up in Control",]$avg_log2FC),],10))
      deg$label <- ifelse(deg$gene %in% top$gene, deg$gene, "")
      plot_ly(deg,x=~avg_log2FC,y=~logP,color=~sig,
        colors=c("Up in Disease"="#D73027","Up in Control"="#4575B4","NS"="grey80"),
        type='scatter',mode='markers',marker=list(size=3,opacity=0.5),
        text=~paste(gene,"<br>log2FC:",round(avg_log2FC,3),"<br>padj:",format.pval(p_val_adj,digits=2)),
        hoverinfo='text') %>%
        add_annotations(x=deg$avg_log2FC[deg$label!=""],y=deg$logP[deg$label!=""],
          text=deg$label[deg$label!=""],showarrow=F,font=list(size=10)) %>%
        layout(title=list(text=input$disease_comp,font=list(size=16)),
               xaxis=list(title="log2 Fold Change"),yaxis=list(title="-log10(adjusted P)"))
    } else {
      plot_ly(prop_df,x=~Disease,y=~Proportion,color=~Subtype,colors=subtype_cols,
        type='bar',text=~paste(Subtype,":",round(Proportion*100,1),"%"),hoverinfo='text') %>%
        layout(title=list(text="Subtype Proportions by Disease Group",font=list(size=16)),
               xaxis=list(title=""),yaxis=list(title="Proportion"),barmode='stack')
    }
  })

  # ── ATAC ──
  output$atac_motif_plot <- renderPlotly({
    df <- atac_motif[atac_motif$sig==TRUE,]; df <- df[order(-df$delta),]; df <- head(df,12)
    df$name <- factor(df$name, levels=rev(df$name))
    plot_ly(df, x=~delta, y=~name, type='bar', orientation='h',
      marker=list(color=~delta, colorscale=list(c(0,"#FDAE61"),c(1,"#D73027"))),
      text=~paste(name,"<br>Δ:",round(delta,1),"%<br>p:",format.pval(pval,digits=2)),hoverinfo='text') %>%
      layout(title=list(text="TF Motif Enrichment in Term-Specific Open Peaks",font=list(size=16)),
             xaxis=list(title="Δ (% Term − Mid)"), yaxis=list(title=""))
  })
  output$atac_volcano_plot <- renderPlotly({
    df <- atac_peaks[!is.na(atac_peaks$p_val_adj),]; df <- df[order(df$p_val_adj),]
    if(nrow(df)>8000) df <- df[1:8000,]
    df$logP <- -log10(df$p_val_adj+1e-300); df$dir <- "NS"
    df$dir[df$sig=="Sig" & df$avg_log2FC>0.25] <- "Term"
    df$dir[df$sig=="Sig" & df$avg_log2FC< -0.25] <- "Mid"
    plot_ly(df,x=~avg_log2FC,y=~logP,color=~dir,
      colors=c("Term"="#D73027","Mid"="#4575B4","NS"="grey80"),
      type='scatter',mode='markers',marker=list(size=2,opacity=0.5),
      text=~paste("log2FC:",round(avg_log2FC,3),"<br>padj:",format.pval(p_val_adj,digits=2)),
      hoverinfo='text') %>%
      layout(title=list(text="Differential Accessibility: Term vs Mid",font=list(size=16)),
             xaxis=list(title="log2 Fold Change"),yaxis=list(title="-log10(adjusted P)"))
  })

  # ── Spatial ──
  output$spatial_slice_img <- renderUI({
    req(input$spatial_sample)
    tags$div(style="padding:10px;",
      tags$p(tags$b(sprintf("Stereo-seq — Sample %s",input$spatial_sample)),
        tags$br(),"Hofbauer cells (red). Mid-gestation placental basal plate, bin50 resolution."),
      tags$img(src=sprintf("Fig3a_spatial_%s.png",input$spatial_sample),
               style="max-width:100%;height:auto;border:1px solid #ddd;border-radius:4px;"))
  })
}

shinyApp(ui, server)
