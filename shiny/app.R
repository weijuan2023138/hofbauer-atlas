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

# ── UI ──
ui <- page_navbar(
  title = "Hofbauer Cell Atlas",
  theme = bs_theme(version=5, bootswatch="flatly", primary="#1565C0",
    "nav-link-font-size" = "0.9rem"),
  
  # ═══ HOME ═══
  nav_panel("Home",
    layout_column_wrap(width=1/2, heights_equal="row",
      card(full_screen=FALSE,
        card_header("About this Atlas"),
        tags$p("An integrated single-cell atlas of human placental Hofbauer cells 
                (fetal tissue-resident macrophages) spanning early to late gestation. 
                We profiled 17,896 high-confidence Hofbauer cells from 9 independent 
                scRNA-seq datasets, complemented by snATAC-seq and Stereo-seq spatial 
                transcriptomics from the UCSF Li 2026 cohort.")),
      card(full_screen=FALSE,
        card_header("Key Numbers"),
        layout_column_wrap(width=1/3,
          value_box("Cells", "17,896", theme="primary"),
          value_box("Subtypes", "6", theme="success"),
          value_box("Datasets", "9", theme="warning"))),
      card(full_screen=FALSE,
        card_header("Data Modalities"),
        tags$ul(
          tags$li(tags$b("scRNA-seq:"), "17,896 Hofbauer cells × 19,940 genes"),
          tags$li(tags$b("snATAC-seq:"), "1,230 mid-gestation + 925 term Hofbauer cells (UCSF Li 2026)"),
          tags$li(tags$b("Stereo-seq:"), "16 mid-gestation placental sections at 500 nm resolution"))),
      card(full_screen=FALSE,
        card_header("Six Transcriptional Subtypes"),
        tags$div(style="display:flex;flex-wrap:wrap;gap:8px;",
          lapply(names(subtype_cols), function(s)
            tags$span(style=sprintf("background:%s;color:white;padding:4px 10px;border-radius:4px;font-size:0.85em",subtype_cols[s]), s)))
      ))
  ),
  
  # ═══ ATLAS ═══
  nav_panel("Atlas",
    layout_sidebar(
      sidebar=sidebar(
        selectInput("atlas_color", "Color by",
          choices=c("Subtype"="subtype","Trimester"="trimester",
                    "Disease"="disease_group","Dataset"="dataset")),
        checkboxInput("atlas_downsample","Subsample to 5,000 cells",TRUE), width=250),
      card(full_screen=TRUE, plotlyOutput("atlas_umap",height="680px"))
    )),
  
  # ═══ GENE EXPRESSION ═══
  nav_panel("Gene Expression",
    layout_sidebar(
      sidebar=sidebar(
        selectizeInput("gene_search","Search gene (19,940 available)",
          choices=NULL, selected="FOLR2",
          options=list(placeholder='Type gene name...',maxOptions=2000)),
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
  
  # ═══ TF REGULATION ═══
  nav_panel("TF Regulation",
    layout_sidebar(
      sidebar=sidebar(
        selectInput("tf_select","Transcription Factor",
          choices=sort(rownames(tf_expr)), selected="CEBPA"),
        selectInput("tf_group","Group by",
          choices=c("Subtype"="subtype","Disease"="disease_group")), width=250),
      layout_column_wrap(width=1/2,
        card(full_screen=TRUE, card_header("UMAP"), plotlyOutput("tf_umap",height="500px")),
        card(full_screen=TRUE, card_header("Expression Distribution"), plotlyOutput("tf_violin",height="500px"))
      ))),
  
  # ═══ ECM-IMMUNE ═══
  nav_panel("ECM-Immune",
    layout_sidebar(
      sidebar=sidebar(
        selectInput("mod_color","Color by",
          choices=c("Disease"="disease_group","Subtype"="subtype")),
        checkboxInput("mod_ellipse","95% confidence ellipses",TRUE), width=250),
      card(full_screen=TRUE, plotlyOutput("mod_scatter",height="680px"))
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
          tags$p(tags$b("ATAC-seq fragment coverage at key gene loci (Term vs Mid)")),
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
          tags$p(tags$b("Magnified view of a representative placental villus. 
            Hofbauer cells (red) localize near fetal capillaries, 
            surrounded by fibroblasts (orange).")),
          tags$img(src="Fig3b_villus_zoom.png",style="max-width:100%;height:auto"))),
      nav_panel("Neighborhood Enrichment",
        tags$div(style="padding:10px;",
          tags$p(tags$b("Cell-type enrichment within 50 μm of Hofbauer cells. 
            FB (2.77×) and fVEC/fEC (1.91×) are the most significantly enriched neighbors.")),
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
        tags$p("All analysis code is available on GitHub. 
                The full Seurat object (Hofbauer_Atlas_Final.rds, 229 MB) 
                is available upon request from the corresponding author."))))
)

# ── Server ──
server <- function(input, output, session) {
  updateSelectizeInput(session,"gene_search",choices=gene_list,server=TRUE)
  
  # ── Atlas ──
  output$atlas_umap <- renderPlotly({
    n <- if(input$atlas_downsample) 5000 else nrow(umap_meta)
    df <- umap_meta; if(n < nrow(df)) df <- df[sample(nrow(df),n),]
    col <- input$atlas_color
    pal <- switch(col, subtype=subtype_cols, disease_group=disease_cols,
                  trimester=trimester_cols, dataset=NULL)
    plot_ly(df, x=~UMAP_1, y=~UMAP_2, color=~get(col), colors=pal,
      type='scatter', mode='markers', marker=list(size=3, opacity=0.7),
      text=~paste(col,":",get(col)), hoverinfo='text') %>%
      layout(title=list(text=paste("Hofbauer Atlas —",col),font=list(size=16)),
             xaxis=list(title="UMAP 1", showgrid=F, zeroline=F),
             yaxis=list(title="UMAP 2", showgrid=F, zeroline=F),
             legend=list(orientation='v', y=0.5, font=list(size=11)))
  })
  
  # ── Gene ──
  output$gene_umap <- renderPlotly({
    req(input$gene_search); gene <- input$gene_search
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
    req(input$gene_search); gene <- input$gene_search
    if(!gene %in% rownames(expr)) return(NULL)
    grp <- input$gene_group
    df <- data.frame(Group=umap_meta[[grp]], Expression=as.numeric(expr[gene,]))
    pal <- if(grp=="subtype") subtype_cols else trimester_cols
    plot_ly(df, x=~Group, y=~Expression, color=~Group, colors=pal,
      type='violin', box=list(visible=T, width=0.1),
      points='all', jitter=0.3, pointpos=-0.5, marker=list(size=1, opacity=0.3)) %>%
      layout(title=list(text=gene,font=list(size=16)),
             xaxis=list(title=""), yaxis=list(title="Log-normalized expression"),
             showlegend=F)
  })
  
  # ── Disease ──
  output$disease_plot <- renderPlotly({
    req(input$disease_comp, input$disease_view)
    if(input$disease_view == "volcano") {
      deg <- read.csv(deg_list[[input$disease_comp]])
      deg <- deg[!is.na(deg$p_val_adj),]
      deg$logP <- -log10(deg$p_val_adj + 1e-300)
      deg$sig <- "NS"
      deg$sig[deg$p_val_adj < 0.05 & abs(deg$avg_log2FC) > 0.5 & deg$avg_log2FC > 0] <- "Up in Disease"
      deg$sig[deg$p_val_adj < 0.05 & abs(deg$avg_log2FC) > 0.5 & deg$avg_log2FC < 0] <- "Up in Control"
      top <- rbind(head(deg[deg$sig=="Up in Disease",][order(-deg[deg$sig=="Up in Disease",]$avg_log2FC),],10),
                   head(deg[deg$sig=="Up in Control",][order(deg[deg$sig=="Up in Control",]$avg_log2FC),],10))
      deg$label <- ifelse(deg$gene %in% top$gene, deg$gene, "")
      plot_ly(deg, x=~avg_log2FC, y=~logP, color=~sig,
        colors=c("Up in Disease"="#D73027","Up in Control"="#4575B4","NS"="grey80"),
        type='scatter', mode='markers', marker=list(size=3, opacity=0.5),
        text=~paste(gene,"<br>log2FC:",round(avg_log2FC,3),"<br>padj:",format.pval(p_val_adj,digits=2)),
        hoverinfo='text') %>%
        add_annotations(x=deg$avg_log2FC[deg$label!=""], y=deg$logP[deg$label!=""],
          text=deg$label[deg$label!=""], showarrow=F, font=list(size=10)) %>%
        layout(title=list(text=input$disease_comp,font=list(size=16)),
               xaxis=list(title="log2 Fold Change"), yaxis=list(title="-log10(adjusted P)"))
    } else {
      plot_ly(prop_df, x=~Disease, y=~Proportion, color=~Subtype, colors=subtype_cols,
        type='bar', text=~paste(Subtype,":",round(Proportion*100,1),"%"), hoverinfo='text') %>%
        layout(title=list(text="Subtype Proportions by Disease Group",font=list(size=16)),
               xaxis=list(title=""), yaxis=list(title="Proportion"), barmode='stack')
    }
  })
  
  # ── TF ──
  output$tf_umap <- renderPlotly({
    req(input$tf_select); tf <- input$tf_select
    df <- umap_meta; df$TF <- as.numeric(tf_expr[tf,])
    plot_ly(df, x=~UMAP_1, y=~UMAP_2, color=~TF,
      colors=colorRamp(c("grey90","#FDAE61","#D73027")),
      type='scatter', mode='markers', marker=list(size=3, opacity=0.7),
      text=~paste(tf,":",round(TF,3)), hoverinfo='text') %>%
      layout(title=list(text=tf,font=list(size=16)),
             xaxis=list(title="UMAP 1"), yaxis=list(title="UMAP 2"))
  })
  output$tf_violin <- renderPlotly({
    req(input$tf_select); tf <- input$tf_select
    grp <- input$tf_group
    df <- data.frame(Group=umap_meta[[grp]], Expression=as.numeric(tf_expr[tf,]))
    plot_ly(df, x=~Group, y=~Expression, color=~Group,
      type='violin', box=list(visible=T, width=0.1),
      points='all', jitter=0.3, pointpos=-0.5, marker=list(size=1, opacity=0.3)) %>%
      layout(title=list(text=tf,font=list(size=16)),
             xaxis=list(title=""), yaxis=list(title="Expression"), showlegend=F)
  })
  
  # ── ECM-Immune ──
  output$mod_scatter <- renderPlotly({
    df <- mod_scores; df <- df[complete.cases(df),]
    col <- input$mod_color; set.seed(42)
    if(nrow(df)>4000) df <- df[sample(nrow(df),4000),]
    plot_ly(df, x=~ECM_score, y=~Immune_score, color=~get(col),
      type='scatter', mode='markers', marker=list(size=4, opacity=0.6),
      text=~paste(col,":",get(col)), hoverinfo='text') %>%
      layout(title=list(text="ECM vs Immune Module Landscape",font=list(size=16)),
             xaxis=list(title="ECM module score"), yaxis=list(title="Immune module score"))
  })
  
  # ── ATAC ──
  output$atac_motif_plot <- renderPlotly({
    df <- atac_motif[atac_motif$sig==TRUE,]; df <- df[order(-df$delta),]; df <- head(df,12)
    df$name <- factor(df$name, levels=rev(df$name))
    plot_ly(df, x=~delta, y=~name, type='bar', orientation='h',
      marker=list(color=~delta, colorscale=list(c(0,"#FDAE61"),c(1,"#D73027"))),
      text=~paste(name,"<br>Δ:",round(delta,1),"%<br>p:",format.pval(pval,digits=2)),
      hoverinfo='text') %>%
      layout(title=list(text="TF Motif Enrichment in Term-Specific Open Peaks",font=list(size=16)),
             xaxis=list(title="Δ (% Term − Mid)"), yaxis=list(title=""))
  })
  output$atac_volcano_plot <- renderPlotly({
    df <- atac_peaks[!is.na(atac_peaks$p_val_adj),]; df <- df[order(df$p_val_adj),]
    if(nrow(df)>8000) df <- df[1:8000,]
    df$logP <- -log10(df$p_val_adj + 1e-300)
    df$dir <- "NS"
    df$dir[df$sig=="Sig" & df$avg_log2FC > 0.25] <- "Term"
    df$dir[df$sig=="Sig" & df$avg_log2FC < -0.25] <- "Mid"
    plot_ly(df, x=~avg_log2FC, y=~logP, color=~dir,
      colors=c("Term"="#D73027","Mid"="#4575B4","NS"="grey80"),
      type='scatter', mode='markers', marker=list(size=2, opacity=0.5),
      text=~paste("log2FC:",round(avg_log2FC,3),"<br>padj:",format.pval(p_val_adj,digits=2)),
      hoverinfo='text') %>%
      layout(title=list(text="Differential Accessibility: Term vs Mid",font=list(size=16)),
             xaxis=list(title="log2 Fold Change (Term / Mid)"),
             yaxis=list(title="-log10(adjusted P)"))
  })
  
  # ── Spatial ──
  output$spatial_slice_img <- renderUI({
    req(input$spatial_sample)
    tags$div(style="padding:10px;",
      tags$p(tags$b(sprintf("Stereo-seq spatial transcriptomics — Sample %s", input$spatial_sample)),
        tags$br(), "Hofbauer cells highlighted in red. Mid-gestation placental basal plate, bin50 resolution."),
      tags$img(src=sprintf("Fig3a_spatial_%s.png", input$spatial_sample),
               style="max-width:100%;height:auto;border:1px solid #ddd;border-radius:4px;"))
  })
}

shinyApp(ui, server)
