# Hofbauer Macrophage Atlas — Shiny App

Interactive single-cell atlas of Hofbauer placental macrophages across pregnancy complications.

**Live:** https://sc-macrophage.shinyapps.io/hofbauer-atlas/

## Features

- **Overview** — UMAP + subtype composition across disease groups
- **Subtype Atlas** — Per-subtype distribution and marker violin plots
- **Gene Lookup** — Search 20K+ genes, UMAP/violin/heatmap expression
- **Disease Comparison** — Compare subtype proportions between two groups
- **Download** — Export metadata, markers, and proportions as CSV

## Data

15,526 Hofbauer macrophages (10,903 main + 4,623 supplementary) across 7 datasets:

| Dataset | Disease Groups |
|---------|---------------|
| Arutyunyan 2023 | Normal 1st trimester |
| GSE290578 | Normal 3rd trimester, Preeclampsia |
| gse214607 | Normal 1st trimester, Miscarriage |
| hoo_2024 | Normal 1st trimester, Infection |
| gse173193 | Preeclampsia |
| gse298119 | Preeclampsia |
| gse183338 | Preeclampsia |

## Deployment

Hosted on [shinyapps.io](https://www.shinyapps.io/) (Posit). To deploy:

```r
rsconnect::setAccountInfo(name='<name>', token='<token>', secret='<secret>')
rsconnect::deployApp('hofbauer-atlas')
```

## Local Run

```bash
cd hofbauer-atlas
Rscript app.R
# → http://localhost:3838
```

Requires: R ≥ 4.3, packages: shiny, ggplot2, dplyr, plotly, DT, shinythemes, tidyr, Matrix

## Files

- `app.R` — Shiny application
- `*.rds` — Pre-computed data (not in repo, generated from Seurat objects)
