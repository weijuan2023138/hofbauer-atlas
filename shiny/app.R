library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(plotly)
library(DT)
library(Matrix)

# ── Load pre-computed data ──
DATA_DIR <- "shiny_data"
umap_meta <- read.csv(file.path(DATA_DIR, "umap_meta.csv"))
expr      <- readRDS(file.path(DATA_DIR, "expr_full.rds"))
tf_expr   <- readRDS(file.path(DATA_DIR, "tf_expr.rds"))
mod_scores <- read.csv(file.path(DATA_DIR, "module_scores.csv"))
prop_df   <- read.csv(file.path(DATA_DIR, "subtype_proportions.csv"))
colnames(prop_df) <- c("Subtype","Disease","Proportion")

gene_list <- sort(rownames(expr))

# ATAC data
atac_motif <- read.csv(file.path(DATA_DIR, "atac_motif_enrichment.csv"))
atac_peaks <- read.csv(file.path(DATA_DIR, "atac_differential_peaks.csv"))

# Disease DEG files
deg_files <- list.files(DATA_DIR, pattern="^deg_", full.names=TRUE)
deg_names <- gsub("deg_|\\.csv","", basename(deg_files))
deg_list <- setNames(deg_files, deg_names)

# ── Color palettes ──
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
  theme = bs_theme(version=5, bootswatch="flatly", primary="#1565C0"),
  
  # ── Tab 1: Atlas ──
  nav_panel("发育图谱",
    layout_sidebar(
      sidebar=sidebar(
        selectInput("atlas_color","着色方式",
          choices=c("亚型 (Subtype)"="subtype","孕期 (Trimester)"="trimester",
                    "疾病分组 (Disease)"="disease_group","数据集 (Dataset)"="dataset")),
        checkboxInput("atlas_downsample","下采样至5000细胞",TRUE), width=250),
      card(full_screen=TRUE, plotlyOutput("atlas_umap",height="650px")))),

  # ── Tab 2: Gene ──
  nav_panel("基因表达",
    layout_sidebar(
      sidebar=sidebar(
        selectizeInput("gene_search","搜索基因",choices=NULL,selected="FOLR2",
          options=list(placeholder='输入基因名...',maxOptions=2000)),
        selectInput("gene_group","分组方式",
          choices=c("亚型 (Subtype)"="subtype","孕期 (Trimester)"="trimester")),width=250),
      layout_column_wrap(width=1/2,
        card(full_screen=TRUE,plotlyOutput("gene_umap",height="500px")),
        card(full_screen=TRUE,plotlyOutput("gene_violin",height="500px"))))),

  # ── Tab 3: Disease ──
  nav_panel("疾病对比",
    layout_sidebar(
      sidebar=sidebar(
        selectInput("disease_comp","疾病比较",
          choices=setNames(deg_names,gsub("_vs_"," vs ",deg_names))),
        selectInput("disease_view","视图",
          choices=c("火山图 (Volcano)"="volcano","亚型比例 (Proportions)"="prop")),width=250),
      card(full_screen=TRUE,plotlyOutput("disease_plot",height="600px")))),

  # ── Tab 4: TF ──
  nav_panel("TF调控",
    layout_sidebar(
      sidebar=sidebar(
        selectInput("tf_select","转录因子",choices=sort(rownames(tf_expr)),selected="CEBPA"),
        selectInput("tf_group","分组方式",
          choices=c("亚型 (Subtype)"="subtype","疾病分组 (Disease)"="disease_group")),width=250),
      layout_column_wrap(width=1/2,
        card(full_screen=TRUE,plotlyOutput("tf_umap",height="500px")),
        card(full_screen=TRUE,plotlyOutput("tf_violin",height="500px"))))),

  # ── Tab 5: ECM-Immune ──
  nav_panel("ECM-免疫双维空间",
    layout_sidebar(
      sidebar=sidebar(
        selectInput("mod_color","着色方式",
          choices=c("疾病分组 (Disease)"="disease_group","亚型 (Subtype)"="subtype")),
        checkboxInput("mod_ellipse","显示95%置信椭圆",TRUE),width=250),
      card(full_screen=TRUE,plotlyOutput("mod_scatter",height="650px")))),

  # ── Tab 6: ATAC-seq ──
  nav_panel("ATAC-seq",
    layout_sidebar(
      sidebar=sidebar(
        radioButtons("atac_view","视图",
          choices=c("Motif富集 (Motif Enrichment)"="motif",
                    "差异可及性火山图 (DA Volcano)"="volcano",
                    "基因座覆盖度轨迹 (Coverage Tracks)"="tracks")),width=250),
      card(full_screen=TRUE,
        conditionalPanel("input.atac_view=='motif'", plotlyOutput("atac_motif_plot",height="550px")),
        conditionalPanel("input.atac_view=='volcano'", plotlyOutput("atac_volcano_plot",height="600px")),
        conditionalPanel("input.atac_view=='tracks'",
          tags$img(src="Fig6b_ATAC_tracks.png",style="max-width:100%;height:auto"))))),

  # ── Tab 7: Spatial ──
  nav_panel("空间转录组",
    layout_sidebar(
      sidebar=sidebar(
        selectInput("spatial_view","视图",
          choices=c("空间切片 (Spatial Slices)"="slices",
                    "绒毛放大 (Villus Zoom)"="villus",
                    "邻域富集 (Neighborhood)"="neighbor")),
        conditionalPanel("input.spatial_view=='slices'",
          selectInput("spatial_sample","样本",
            choices=c("Sample 001"="001","Sample 004"="004","Sample 010"="010","Sample 014"="014"))),
        width=250),
      card(full_screen=TRUE,
        conditionalPanel("input.spatial_view=='slices'",
          uiOutput("spatial_slice_img")),
        conditionalPanel("input.spatial_view=='villus'",
          tags$img(src="Fig3b_villus_zoom.png",style="max-width:100%;height:auto")),
        conditionalPanel("input.spatial_view=='neighbor'",
          tags$img(src="Fig3C_neighborhood.png",style="max-width:100%;height:auto"))))),

  # ── Tab 8: Downloads ──
  nav_panel("数据下载",
    card(
      card_header("预计算数据"),
      tags$ul(
        tags$li(tags$a("UMAP坐标与元数据 (CSV)",href="shiny_data/umap_meta.csv",download=NA)),
        tags$li(tags$a("亚型比例 (CSV)",href="shiny_data/subtype_proportions.csv",download=NA)),
        tags$li(tags$a("模块评分 (CSV)",href="shiny_data/module_scores.csv",download=NA)),
        tags$li(tags$a("ATAC Motif富集 (CSV)",href="shiny_data/atac_motif_enrichment.csv",download=NA)),
        tags$li(tags$a("ATAC差异Peak (CSV)",href="shiny_data/atac_differential_peaks.csv",download=NA))),
      card_header("差异表达基因"),
      tags$ul(lapply(deg_names,function(d)
        tags$li(tags$a(sprintf("%s (CSV)",d),href=sprintf("shiny_data/deg_%s.csv",d),download=NA)))),
      card_header("完整分析代码"),
      tags$p("所有分析代码可在 GitHub 获取。RDS对象（229MB）可向通讯作者请求。"))))

# ── Server ──
server <- function(input, output, session) {
  updateSelectizeInput(session,"gene_search",choices=gene_list,server=TRUE)

  # ── Tab 1: Atlas ──
  output$atlas_umap <- renderPlotly({
    n <- if(input$atlas_downsample) 5000 else nrow(umap_meta)
    df <- umap_meta; if(n<nrow(df)) df <- df[sample(nrow(df),n),]
    col <- input$atlas_color
    pal <- switch(col, "subtype"=subtype_cols, "disease_group"=disease_cols,
                  "trimester"=trimester_cols, "dataset"=NULL)
    plot_ly(df,x=~UMAP_1,y=~UMAP_2,color=~get(col),colors=pal,
      type='scatter',mode='markers',marker=list(size=3,opacity=0.7),
      text=~paste(col,":",get(col)),hoverinfo='text') %>%
      layout(title=paste("Hofbauer Atlas —",col),
             xaxis=list(title="UMAP 1",showgrid=F,zeroline=F),
             yaxis=list(title="UMAP 2",showgrid=F,zeroline=F),
             legend=list(orientation='v',y=0.5))
  })

  # ── Tab 2: Gene ──
  output$gene_umap <- renderPlotly({
    req(input$gene_search); gene <- input$gene_search
    if(!gene %in% rownames(expr)) return(NULL)
    df <- umap_meta; df$Expression <- as.numeric(expr[gene,])
    plot_ly(df,x=~UMAP_1,y=~UMAP_2,color=~Expression,
      colors=colorRamp(c("grey90","#BDD7E7","#2171B5","#08306B")),
      type='scatter',mode='markers',marker=list(size=3,opacity=0.7),
      text=~paste(gene,":",round(Expression,3)),hoverinfo='text') %>%
      layout(title=paste(gene,"— UMAP"),xaxis=list(title="UMAP 1"),yaxis=list(title="UMAP 2"))
  })
  output$gene_violin <- renderPlotly({
    req(input$gene_search); gene <- input$gene_search
    if(!gene %in% rownames(expr)) return(NULL)
    grp <- input$gene_group
    df <- data.frame(Group=umap_meta[[grp]],Expression=as.numeric(expr[gene,]))
    pal <- if(grp=="subtype") subtype_cols else trimester_cols
    plot_ly(df,x=~Group,y=~Expression,color=~Group,colors=pal,
      type='violin',box=list(visible=T,width=0.1),
      points='all',jitter=0.3,pointpos=-0.5,marker=list(size=1,opacity=0.3)) %>%
      layout(title=paste(gene,"— by",grp),xaxis=list(title=""),
             yaxis=list(title="Log-normalized expression"),showlegend=F)
  })

  # ── Tab 3: Disease ──
  output$disease_plot <- renderPlotly({
    req(input$disease_comp,input$disease_view)
    if(input$disease_view=="volcano"){
      deg <- read.csv(deg_list[[input$disease_comp]])
      deg <- deg[!is.na(deg$p_val_adj),]
      deg$logP <- -log10(deg$p_val_adj)
      deg$sig <- "NS"
      deg$sig[deg$p_val_adj<0.05 & abs(deg$avg_log2FC)>0.5 & deg$avg_log2FC>0] <- "Up"
      deg$sig[deg$p_val_adj<0.05 & abs(deg$avg_log2FC)>0.5 & deg$avg_log2FC<0] <- "Down"
      top <- rbind(head(deg[deg$sig=="Up",][order(-deg[deg$sig=="Up",]$avg_log2FC),],10),
                   head(deg[deg$sig=="Down",][order(deg[deg$sig=="Down",]$avg_log2FC),],10))
      deg$label <- ifelse(deg$gene %in% top$gene, deg$gene, "")
      plot_ly(deg,x=~avg_log2FC,y=~logP,color=~sig,
        colors=c("Up"="#D73027","Down"="#4575B4","NS"="grey80"),
        type='scatter',mode='markers',marker=list(size=3,opacity=0.5),
        text=~paste(gene,"<br>log2FC:",round(avg_log2FC,3),"<br>padj:",format.pval(p_val_adj,digits=2)),
        hoverinfo='text') %>%
        add_annotations(x=deg$avg_log2FC[deg$label!=""],y=deg$logP[deg$label!=""],
          text=deg$label[deg$label!=""],showarrow=F,font=list(size=10)) %>%
        layout(title=input$disease_comp,xaxis=list(title="log2 Fold Change"),
               yaxis=list(title="-log10(adjusted P)"))
    } else {
      plot_ly(prop_df,x=~Disease,y=~Proportion,color=~Subtype,colors=subtype_cols,
        type='bar',text=~paste(Subtype,":",round(Proportion*100,1),"%"),hoverinfo='text') %>%
        layout(title="Subtype proportions by disease group",
               xaxis=list(title=""),yaxis=list(title="Proportion"),barmode='stack')
    }
  })

  # ── Tab 4: TF ──
  output$tf_umap <- renderPlotly({
    req(input$tf_select); tf <- input$tf_select
    df <- umap_meta; df$TF <- as.numeric(tf_expr[tf,])
    plot_ly(df,x=~UMAP_1,y=~UMAP_2,color=~TF,
      colors=colorRamp(c("grey90","#FDAE61","#D73027")),
      type='scatter',mode='markers',marker=list(size=3,opacity=0.7),
      text=~paste(tf,":",round(TF,3)),hoverinfo='text') %>%
      layout(title=paste(tf,"— UMAP"),xaxis=list(title="UMAP 1"),yaxis=list(title="UMAP 2"))
  })
  output$tf_violin <- renderPlotly({
    req(input$tf_select); tf <- input$tf_select
    grp <- input$tf_group
    df <- data.frame(Group=umap_meta[[grp]],Expression=as.numeric(tf_expr[tf,]))
    plot_ly(df,x=~Group,y=~Expression,color=~Group,
      type='violin',box=list(visible=T,width=0.1),
      points='all',jitter=0.3,pointpos=-0.5,marker=list(size=1,opacity=0.3)) %>%
      layout(title=paste(tf,"— by",grp),xaxis=list(title=""),yaxis=list(title="Expression"),showlegend=F)
  })

  # ── Tab 5: ECM-Immune ──
  output$mod_scatter <- renderPlotly({
    df <- mod_scores; df <- df[complete.cases(df),]
    col <- input$mod_color; set.seed(42)
    if(nrow(df)>4000) df <- df[sample(nrow(df),4000),]
    plot_ly(df,x=~ECM_score,y=~Immune_score,color=~get(col),
      type='scatter',mode='markers',marker=list(size=4,opacity=0.6),
      text=~paste(col,":",get(col)),hoverinfo='text') %>%
      layout(title="ECM vs Immune module landscape",
             xaxis=list(title="ECM module score"),yaxis=list(title="Immune module score"))
  })

  # ── Tab 6: ATAC ──
  output$atac_motif_plot <- renderPlotly({
    df <- atac_motif[atac_motif$sig==TRUE,]
    df <- df[order(-df$delta),]
    df <- head(df, 12)
    df$name <- factor(df$name, levels=rev(df$name))
    plot_ly(df, x=~delta, y=~name, type='bar', orientation='h',
      marker=list(color=~delta, colorscale=list(c(0,"#FDAE61"),c(1,"#D73027"))),
      text=~paste(name,"<br>Δ:",round(delta,1),"%<br>p:",format.pval(pval,digits=2)),
      hoverinfo='text') %>%
      layout(title="ATAC Motif Enrichment — Term vs Mid specific peaks",
             xaxis=list(title="Δ (% Term − Mid)"), yaxis=list(title=""))
  })

  output$atac_volcano_plot <- renderPlotly({
    df <- atac_peaks[!is.na(atac_peaks$p_val_adj),]
    df <- df[order(df$p_val_adj),]
    if(nrow(df)>8000) df <- df[1:8000,]
    df$logP <- -log10(df$p_val_adj + 1e-300)
    df$dir <- "NS"
    df$dir[df$sig=="Sig" & df$avg_log2FC > 0.25] <- "Term-up"
    df$dir[df$sig=="Sig" & df$avg_log2FC < -0.25] <- "Mid-up"
    plot_ly(df, x=~avg_log2FC, y=~logP, color=~dir,
      colors=c("Term-up"="#D73027","Mid-up"="#4575B4","NS"="grey80"),
      type='scatter',mode='markers',marker=list(size=2,opacity=0.5),
      text=~paste("log2FC:",round(avg_log2FC,3),"<br>padj:",format.pval(p_val_adj,digits=2)),
      hoverinfo='text') %>%
      layout(title="Differential Accessibility: Term vs Mid",
             xaxis=list(title="log2 Fold Change (Term/Mid)"),
             yaxis=list(title="-log10(adjusted P)"))
  })

  # ── Tab 7: Spatial ──
  output$spatial_slice_img <- renderUI({
    req(input$spatial_sample)
    tags$img(src=sprintf("Fig3a_spatial_%s.png",input$spatial_sample),
             style="max-width:100%;height:auto")
  })
}

shinyApp(ui, server)
