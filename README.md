# Hofbauer Cell Atlas — Multi-Omic Developmental Atlas Across Gestation and Disease

Integration of 11 scRNA-seq datasets (21,681 Hofbauer cells) + snATAC-seq + Stereo-seq spatial transcriptomics + CellChat analysis, constructing a multi-omic atlas of human placental Hofbauer cells.

**Manuscript**: "The architecture of Hofbauer cells" (submitted to *Communications Biology*)  
**Interactive Data Portal**: https://sc-macrophage.shinyapps.io/hofbauer-atlas/  
**Manuscript Files**: `../文章手稿/Manuscript_Full_CommsBio.md` (markdown) / `.docx` (Word)

## Repository Structure

```
ucsf_integration/
├── code/                  # Analysis scripts (R + Python) — main pipeline
├── vento_integration/     # Vento-Tormo integration sub-project
│   └── code/              # Sub-pipeline scripts (R + Python, 64 scripts)
├── results/               # Output data (RDS, CSV, h5ad)
│   └── classification/    # Classifier output
├── figures/               # Publication figures (300dpi PNG/PDF)
│   ├── Fig1/              # Developmental atlas + protein validation
│   ├── Fig2/              # TF switch + ATAC-seq
│   ├── Fig3/              # Spatial transcriptomics niche
│   ├── Fig4/              # CellChat communication network
│   ├── Fig5/              # Dual-track TF regulatory model
│   ├── Fig6/              # Disease perturbation analysis
│   ├── Fig7/              # Disease disruption of dual-track model
│   └── FigS/              # Supplementary figures
├── write/                 # Draft results text (Chinese)
├── ref/                   # Reference files (GMT gene sets)
├── deprecated/            # Deprecated intermediate files
│   ├── code/
│   ├── figures/
│   └── results/
├── shiny/                 # Shiny app for interactive data exploration
│   ├── app.R
│   ├── shiny_data/
│   └── www/
├── logs/                  # Run logs
├── app.R                  # Standalone Shiny app
├── 工作记录.md             # Analysis decision log (Chinese)
└── README.md
```

## Figure-to-Script Mapping

| Figure | Key Scripts | Description |
|--------|-------------|-------------|
| Fig1A | `code/42_fig1a_model.R` | Study design schematic |
| Fig1B | `code/17_final_comprehensive_figures.R` | Subtype UMAP |
| Fig1C | `code/fig5a_umap.R` | Trimester UMAP |
| Fig1D | `code/fig5b_subtype_proportions.R` | Subtype proportion line plots |
| Fig1E | `code/fig5c_gsea.R` | GSEA functional transition |
| Fig1F | `code/29_module_volcano.R` | Module score trajectories |
| Fig1G | `code/fig5d_gene_dotplot.R` | Developmental gene dot plot |
| Fig2A | `code/30_fig2a_final.R` | TF dot plot |
| Fig2B | `code/34_atac_analysis.R` | ATAC volcano plot |
| Fig2C | `code/37_atac_rna_joint.R` | ATAC-RNA scatter plot |
| Fig2D | `code/fig6c_motif.R` / `code/56_motif_enrich.py` | Motif enrichment |
| Fig2E | `code/41_fig2b_gsea_network.R` | GSEA network |
| Fig2F | `code/20_GSEA_trimester.R` | Expression heatmap |
| Fig3 | `code/32_spatial_plot.py`, `code/33_neighborhood.py` | Spatial transcriptomics |
| Fig4 | `code/45_fig4_cellchat.R`, `code/63_spatial_LR_validation.py` | CellChat + spatial validation |
| Fig5 | `code/fig6a_tf_comm_corr.R` ~ `code/fig6e_model.R` | Dual-track regulatory model |
| Fig6 | `code/40_disease_analysis.R`, `code/fig5*.R` | Disease perturbation |
| Fig7 | `code/fig7_disease_dual.R`, `code/fig7g_stomics_subtypes.py` | Disease disruption of dual-track |

Additional analysis scripts are in `vento_integration/code/` (64 scripts covering data extraction, classification, Harmony integration, DEG, GSEA, spatial analysis).

## Key Data Files

| File | Content |
|------|---------|
| `results/hofbauer_final_clean.rds` | Final Seurat object (21,681 cells) |
| `results/Hofbauer_Atlas_Final.rds` | Annotated atlas object |
| `results/Hofbauer_ATAC_mid_term.rds` | ATAC Mid vs Term object |
| `results/cellchat_ucsf_mid.rds` | Mid-gestation CellChat results |
| `results/cellchat_hoo2024_*.rds` | Infection vs control CellChat |
| `results/classification/all_hofbauer_final_corrected.h5ad` | Classifier output (21,681 cells) |

## Datasets Integrated

| Dataset | Gestation | Condition | n (HBCs) |
|---------|-----------|-----------|----------|
| E-MTAB-6701 | GW4.5–10 | Normal (elective termination) | 5,002 |
| GSE290578 | GW32–38 | Normal + PE | 6,704 |
| GSE214607 | GW6–10 | Normal + Miscarriage | 874 |
| E-MTAB-12421 | GW5–8 | Normal + Infection | 1,975 |
| GSE173193 | GW29–40 | Normal + PE | 2,731 |
| GSE298119 | GW29–34 | Normal + PE | 3,012 |
| GSE298602 | GW37–40 | Normal + PE | TBD |
| GSE329173 | GW30–38 | Severe PE | 337 |
| GSE333257 | GW32–38 | Normal + PTB | 3,716 |
| E-MTAB-12795 | GW6–14 | Normal | TBD |
| UCSF Li 2026 | GW11–24 | Normal (mid-gestation) | TBD |
| **Total** | **GW4.5–38** | | **21,681** |

## Environment

- R >= 4.3, Python >= 3.10
- R packages: Seurat v5, Signac, harmony, CellChat v2, fgsea, clusterProfiler, slingshot, monocle3, ComplexHeatmap, ggplot2
- Python packages: scanpy, squidpy, pandas, numpy, macs3
- MACS3: `~/.local/bin/macs3`

## Known Limitations

1. Hi-C/ABC model data unavailable — speculative statements in manuscript are noted
2. Stereo-seq data GEO accession pending
3. CellChat only for mid-gestation normal; disease-state CellChat coverage incomplete
4. Infection phenotype from single dataset (hoo_2024, n=1,975); generalizability requires validation
5. snATAC-seq limited to mid-gestation vs term (first-trimester chromatin landscape missing)
