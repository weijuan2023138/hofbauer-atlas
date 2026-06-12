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
expr      <- readRDS(file.path(DATA_DIR, "expr_top2000.rds"))
tf_expr   <- readRDS(file.path(DATA_DIR, "tf_expr.rds"))
mod_scores <- read.csv(file.path(DATA_DIR, "module_scores.csv"))
prop_df   <- read.csv(file.path(DATA_DIR, "subtype_proportions.csv"))
colnames(prop_df) <- c("Subtype","Disease","Proportion")

gene_list <- sort(rownames(expr))

# Disease DEG files
deg_files <- list.files(DATA_DIR, pattern="^deg_", full.names=TRUE)
deg_names <- gsub("deg_|\\.csv","", basename(deg_files))
deg_list <- setNames(deg_files, deg_names)

# ── Color palettes ──
subtype_cols <- c(
  "Pro-inflammatory"="#C62828",
  "MHCII+ Antigen-presenting"="#E65100",
  "Homeostatic"="#1565C0",
  "PRKN+ Autophagy"="#6A1B9A",
  "Vascular remodeling"="#2E7D32",
  "MKI67+ Proliferating"="#455A64"
)
disease_cols <- c(
  "Normal 1st trimester"="#4DBBD5",
  "Normal 1st/2nd/Term"="#00A087",
  "Normal 3rd trimester / Preeclampsia"="#7E6148",
  "Preeclampsia"="#C62828",
  "Preterm Labor"="#E18727",
  "Preterm No Labor"="#BCAAA4",
  "Term Labor"="#3C5488",
  "Miscarriage / Normal"="#F39B7F",
  "Infection"="#DC0000"
)
trimester_cols <- c("Early"="#4575B4","Mid"="#FDAE61","Late"="#D73027")

# ── UI ──
ui <- page_navbar(
  title = "Hofbauer Cell Atlas",
  theme = bs_theme(version=5, bootswatch="flatly", primary="#1565C0"),
  
  # ── Tab 1: Atlas ──
  nav_panel(
    "发育图谱",
    layout_sidebar(
      sidebar = sidebar(
        selectInput("atlas_color", "着色方式",
          choices = c("亚型 (Subtype)"="subtype",
                      "孕期 (Trimester)"="trimester",
                      "疾病分组 (Disease)"="disease_group",
                      "数据集 (Dataset)"="dataset")),
        checkboxInput("atlas_downsample", "下采样至5000细胞 (加快渲染)", TRUE),
        width = 250
      ),
      card(full_screen=TRUE, plotlyOutput("atlas_umap", height="650px"))
    )
  ),
  
  # ── Tab 2: Gene Expression ──
  nav_panel(
    "基因表达",
    layout_sidebar(
      sidebar = sidebar(
        selectizeInput("gene_search", "搜索基因",
          choices = NULL, selected = "FOLR2",
          options = list(placeholder='输入基因名...', maxOptions=2000)),
        selectInput("gene_group", "分组方式",
          choices = c("亚型 (Subtype)"="subtype",
                      "孕期 (Trimester)"="trimester")),
        width = 250
      ),
      layout_column_wrap(
        width = 1/2,
        card(full_screen=TRUE, plotlyOutput("gene_umap", height="500px")),
        card(full_screen=TRUE, plotlyOutput("gene_violin", height="500px"))
      )
    )
  ),
  
  # ── Tab 3: Disease ──
  nav_panel(
    "疾病对比",
    layout_sidebar(
      sidebar = sidebar(
        selectInput("disease_comp", "疾病比较",
          choices = setNames(deg_names, gsub("_vs_"," vs ",deg_names))),
        selectInput("disease_view", "视图",
          choices = c("火山图 (Volcano)"="volcano",
                      "亚型比例 (Subtype Proportions)"="prop")),
        width = 250
      ),
      card(full_screen=TRUE, plotlyOutput("disease_plot", height="600px"))
    )
  ),
  
  # ── Tab 4: TF Regulation ──
  nav_panel(
    "TF调控",
    layout_sidebar(
      sidebar = sidebar(
        selectInput("tf_select", "转录因子",
          choices = sort(rownames(tf_expr)),
          selected = "CEBPA"),
        selectInput("tf_group", "分组方式",
          choices = c("亚型 (Subtype)"="subtype",
                      "疾病分组 (Disease)"="disease_group")),
        width = 250
      ),
      layout_column_wrap(
        width = 1/2,
        card(full_screen=TRUE, plotlyOutput("tf_umap", height="500px")),
        card(full_screen=TRUE, plotlyOutput("tf_violin", height="500px"))
      )
    )
  ),
  
  # ── Tab 5: ECM-Immune ──
  nav_panel(
    "ECM-免疫双维空间",
    layout_sidebar(
      sidebar = sidebar(
        selectInput("mod_color", "着色方式",
          choices = c("疾病分组 (Disease)"="disease_group",
                      "亚型 (Subtype)"="subtype")),
        checkboxInput("mod_ellipse", "显示95%置信椭圆", TRUE),
        width = 250
      ),
      card(full_screen=TRUE, plotlyOutput("mod_scatter", height="650px"))
    )
  ),
  
  # ── Tab 6: Downloads ──
  nav_panel(
    "数据下载",
    card(
      card_header("下载预计算数据"),
      tags$ul(
        tags$li(tags$a("UMAP坐标与元数据 (CSV)", href="shiny_data/umap_meta.csv", download=NA)),
        tags$li(tags$a("亚型比例 (CSV)", href="shiny_data/subtype_proportions.csv", download=NA)),
        tags$li(tags$a("模块评分 (CSV)", href="shiny_data/module_scores.csv", download=NA))
      ),
      card_header("差异表达基因"),
      tags$ul(lapply(deg_names, function(d) {
        tags$li(tags$a(sprintf("%s (CSV)", d), href=sprintf("shiny_data/deg_%s.csv", d), download=NA))
      })),
      card_header("完整分析代码"),
      tags$p("所有分析代码可在 GitHub 获取。RDS对象（Hofbauer_Atlas_Final.rds, 229MB）可向通讯作者请求。")
    )
  )
)

# ── Server ──
server <- function(input, output, session) {
  
  updateSelectizeInput(session, "gene_search", choices = gene_list, server = TRUE)
  
  # Helper: downsample UMAP data
  get_umap <- function(n=5000) {
    df <- umap_meta
    if(n < nrow(df)) df <- df[sample(nrow(df), n), ]
    df
  }
  
  # Helper: build UMAP plotly
  build_umap <- function(df, color_col, palette, title="") {
    p <- plot_ly(df, x=~UMAP_1, y=~UMAP_2, color=~get(color_col),
          colors=palette, type='scatter', mode='markers',
          marker=list(size=3, opacity=0.7),
          text=~paste(color_col,":", get(color_col)),
          hoverinfo='text') %>%
      layout(title=title, xaxis=list(title="UMAP 1", showgrid=FALSE, zeroline=FALSE),
             yaxis=list(title="UMAP 2", showgrid=FALSE, zeroline=FALSE),
             legend=list(orientation='v', y=0.5))
    p
  }
  
  # ── Tab 1: Atlas UMAP ──
  output$atlas_umap <- renderPlotly({
    n <- if(input$atlas_downsample) 5000 else nrow(umap_meta)
    df <- get_umap(n)
    col <- input$atlas_color
    pal <- switch(col,
      "subtype"=subtype_cols,
      "disease_group"=disease_cols,
      "trimester"=trimester_cols,
      "dataset"=NULL)
    build_umap(df, col, pal, paste("Hofbauer Atlas —", col))
  })
  
  # ── Tab 2: Gene expression ──
  output$gene_umap <- renderPlotly({
    req(input$gene_search)
    gene <- input$gene_search
    if(!gene %in% rownames(expr)) return(NULL)
    vals <- as.numeric(expr[gene, ])
    df <- umap_meta; df$Expression <- vals
    p <- plot_ly(df, x=~UMAP_1, y=~UMAP_2, color=~Expression,
          colors=colorRamp(c("grey90","#BDD7E7","#2171B5","#08306B")),
          type='scatter', mode='markers',
          marker=list(size=3, opacity=0.7),
          text=~paste(gene,":", round(Expression,3)),
          hoverinfo='text') %>%
      layout(title=paste(gene,"— UMAP"),
             xaxis=list(title="UMAP 1"), yaxis=list(title="UMAP 2"))
    p
  })
  
  output$gene_violin <- renderPlotly({
    req(input$gene_search)
    gene <- input$gene_search
    if(!gene %in% rownames(expr)) return(NULL)
    vals <- as.numeric(expr[gene, ])
    grp <- input$gene_group
    df <- data.frame(Group=umap_meta[[grp]], Expression=vals)
    pal <- if(grp=="subtype") subtype_cols else trimester_cols
    
    p <- plot_ly(df, x=~Group, y=~Expression, color=~Group, colors=pal,
          type='violin', box=list(visible=TRUE, width=0.1),
          points='all', jitter=0.3, pointpos=-0.5,
          marker=list(size=1, opacity=0.3)) %>%
      layout(title=paste(gene,"— expression by", grp),
             xaxis=list(title=""), yaxis=list(title="Log-normalized expression"),
             showlegend=FALSE)
    p
  })
  
  # ── Tab 3: Disease ──
  output$disease_plot <- renderPlotly({
    req(input$disease_comp, input$disease_view)
    
    if(input$disease_view == "volcano") {
      deg <- read.csv(deg_list[[input$disease_comp]])
      deg <- deg[!is.na(deg$p_val_adj), ]
      deg$logP <- -log10(deg$p_val_adj)
      deg$sig <- "NS"
      deg$sig[deg$p_val_adj < 0.05 & abs(deg$avg_log2FC) > 0.5 & deg$avg_log2FC > 0] <- "Up"
      deg$sig[deg$p_val_adj < 0.05 & abs(deg$avg_log2FC) > 0.5 & deg$avg_log2FC < 0] <- "Down"
      top_genes <- rbind(
        head(deg[deg$sig=="Up", ][order(deg[deg$sig=="Up",]$avg_log2FC, decreasing=TRUE), ], 10),
        head(deg[deg$sig=="Down",][order(deg[deg$sig=="Down",]$avg_log2FC), ], 10)
      )
      deg$label <- ifelse(deg$gene %in% top_genes$gene, deg$gene, "")
      
      p <- plot_ly(deg, x=~avg_log2FC, y=~logP, color=~sig,
            colors=c("Up"="#D73027","Down"="#4575B4","NS"="grey80"),
            type='scatter', mode='markers',
            marker=list(size=3, opacity=0.5),
            text=~paste(gene,"<br>log2FC:",round(avg_log2FC,3),"<br>padj:",format.pval(p_val_adj,digits=2)),
            hoverinfo='text') %>%
        add_annotations(x=deg$avg_log2FC[deg$label!=""], y=deg$logP[deg$label!=""],
          text=deg$label[deg$label!=""], showarrow=FALSE, font=list(size=10)) %>%
        layout(title=input$disease_comp, xaxis=list(title="log2 Fold Change"),
               yaxis=list(title="-log10(adjusted P)"))
      p
    } else {
      p <- plot_ly(prop_df, x=~Disease, y=~Proportion, color=~Subtype,
            colors=subtype_cols, type='bar',
            text=~paste(Subtype,":",round(Proportion*100,1),"%"),
            hoverinfo='text') %>%
        layout(title="Subtype proportions by disease group",
               xaxis=list(title=""), yaxis=list(title="Proportion"),
               barmode='stack')
      p
    }
  })
  
  # ── Tab 4: TF regulation ──
  output$tf_umap <- renderPlotly({
    req(input$tf_select)
    tf <- input$tf_select
    vals <- as.numeric(tf_expr[tf, ])
    df <- umap_meta; df$TF <- vals
    plot_ly(df, x=~UMAP_1, y=~UMAP_2, color=~TF,
      colors=colorRamp(c("grey90","#FDAE61","#D73027")),
      type='scatter', mode='markers',
      marker=list(size=3, opacity=0.7),
      text=~paste(tf,":", round(TF,3)),
      hoverinfo='text') %>%
      layout(title=paste(tf,"— UMAP"), xaxis=list(title="UMAP 1"), yaxis=list(title="UMAP 2"))
  })
  
  output$tf_violin <- renderPlotly({
    req(input$tf_select)
    tf <- input$tf_select
    vals <- as.numeric(tf_expr[tf, ])
    grp <- input$tf_group
    df <- data.frame(Group=umap_meta[[grp]], Expression=vals)
    plot_ly(df, x=~Group, y=~Expression, color=~Group,
      type='violin', box=list(visible=TRUE, width=0.1),
      points='all', jitter=0.3, pointpos=-0.5,
      marker=list(size=1, opacity=0.3)) %>%
      layout(title=paste(tf,"— by", grp), xaxis=list(title=""), yaxis=list(title="Expression"),
             showlegend=FALSE)
  })
  
  # ── Tab 5: ECM-Immune ──
  output$mod_scatter <- renderPlotly({
    df <- mod_scores; df <- df[complete.cases(df), ]
    col <- input$mod_color
    set.seed(42)
    if(nrow(df) > 4000) df <- df[sample(nrow(df), 4000), ]
    
    p <- plot_ly(df, x=~ECM_score, y=~Immune_score, color=~get(col),
          type='scatter', mode='markers',
          marker=list(size=4, opacity=0.6),
          text=~paste(col,":", get(col)),
          hoverinfo='text') %>%
      layout(title="ECM vs Immune module landscape",
             xaxis=list(title="ECM module score"), yaxis=list(title="Immune module score"))
    p
  })
}

shinyApp(ui, server)
