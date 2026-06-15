#!/usr/bin/env python3
"""Merge old Atlas + new datasets in Python, save as h5ad for R."""
import scanpy as sc
import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

BASE = '/home/weijuan/文档/胎盘单细胞数据'
OUT = f'{BASE}/ucsf_integration/vento_integration'

# ── 1. Load old Atlas (from R .rds via Seurat-to-h5ad conversion is hard)
# Instead, load the original merged h5ad that was the INPUT to Harmony
# The old pipeline used hofbauer_corrected_clustered.rds → but we need h5ad
# Alternative: load each old dataset individually from their classified h5ad files

CLASS_DIR = f'{BASE}/ucsf_integration/results/classification'
PROCESSED = f'{BASE}/processed'

# Old datasets with their paths and extraction methods
old_datasets = {
    'E-MTAB-12421': {'path': f'{PROCESSED}/arutyunyan_processed.h5ad', 'col': 'cell_type_fine', 'label': 'Hofbauer'},
    'GSE290578':    {'path': f'{BASE}/results/phase3_gse290578/gse290578_reclassified.h5ad', 'col': 'cell_type_new', 'label': 'Hofbauer'},
    'GSE214607':    {'path': f'{PROCESSED}/gse214607_reclassified.h5ad', 'col': 'cell_type_new', 'label': 'Hofbauer'},
    'E-MTAB-12795': {'path': f'{PROCESSED}/hoo_2024_reclassified.h5ad', 'col': 'cell_type_new', 'label': 'Hofbauer'},
    'GSE173193':    {'path': f'{PROCESSED}/gse173193_reclassified.h5ad', 'col': 'cell_type_new', 'label': 'Hofbauer'},
    'GSE298119':    {'path': f'{PROCESSED}/gse298119_reclassified.h5ad', 'col': 'cell_type_new', 'label': 'Hofbauer'},
    'GSE333257':    {'path': f'{PROCESSED}/my_cohort_processed.h5ad', 'col': 'cell_type_fine', 'label': 'Hofbauer'},
    'UCSF Li 2026': {'path': f'{CLASS_DIR}/UCSF_Li_2026_hofbauer.h5ad', 'col': None, 'label': None},
}

# New datasets (already classified as Hofbauer)
new_datasets = {
    'E-MTAB-6701': f'{CLASS_DIR}/vento_tormo_reclassified_hofbauer.h5ad',
    'GSE298602':   f'{CLASS_DIR}/gse298602_all_hofbauer.h5ad',
}

all_ads = []

# Load old datasets
for ds_name, info in old_datasets.items():
    print(f"Loading {ds_name}...")
    ad = sc.read_h5ad(info['path'])
    
    if info['col'] is not None:
        mask = ad.obs[info['col']] == info['label']
        ad = ad[mask].copy()
    
    ad.obs['dataset'] = ds_name
    all_ads.append(ad)
    print(f"  {ad.n_obs} Hofbauer cells")

# Load new datasets
for ds_name, path in new_datasets.items():
    print(f"Loading {ds_name}...")
    ad = sc.read_h5ad(path)
    ad.obs['dataset'] = ds_name
    all_ads.append(ad)
    print(f"  {ad.n_obs} Hofbauer cells")

# Find common genes
common_genes = set(all_ads[0].var_names)
for ad in all_ads[1:]:
    common_genes &= set(ad.var_names)
print(f"\nCommon genes: {len(common_genes)}")

# Subset to common genes and concat
all_sub = [ad[:, list(common_genes)] for ad in all_ads]
combined = sc.concat(all_sub, join='inner')
print(f"Combined: {combined.n_obs} cells, {combined.n_vars} genes")

# Save
out_path = f'{OUT}/final_10datasets_merged.h5ad'
combined.write_h5ad(out_path)
print(f"Saved: {out_path}")

# Per-dataset count
for ds in sorted(combined.obs['dataset'].unique()):
    print(f"  {ds}: {sum(combined.obs['dataset']==ds)} cells")
print(f"  TOTAL: {combined.n_obs}")
