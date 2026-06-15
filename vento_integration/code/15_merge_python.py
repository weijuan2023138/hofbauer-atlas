#!/usr/bin/env python3
"""Merge old Atlas h5ad + new datasets in Python → single h5ad for R Harmony."""
import scanpy as sc
import numpy as np
import os

BASE = '/home/weijuan/文档/胎盘单细胞数据'
OUT = f'{BASE}/ucsf_integration/vento_integration'
CLASS = f'{BASE}/ucsf_integration/results/classification'

# Load old base
print("Loading old Atlas base...")
old = sc.read_h5ad(f'{OUT}/old_atlas_base.h5ad')
print(f"  {old.n_obs} cells, {old.n_vars} genes")
print(f"  Datasets: {old.obs['dataset'].value_counts().to_dict()}")

# Load new
new_paths = {
    'E-MTAB-6701': f'{CLASS}/vento_tormo_reclassified_hofbauer.h5ad',
    'GSE298602':   f'{CLASS}/gse298602_all_hofbauer.h5ad',
}

new_ads = []
for ds, path in new_paths.items():
    print(f"\nLoading {ds}...")
    ad = sc.read_h5ad(path)
    # Keep only gene expression, drop extra obs columns
    ad.obs['dataset'] = ds
    # Keep only necessary obs columns
    for col in list(ad.obs.columns):
        if col not in ['dataset']:
            del ad.obs[col]
    new_ads.append(ad)
    print(f"  {ad.n_obs} cells, {ad.n_vars} genes")

# Find common genes  
print("\nFinding common genes...")
all_genes = set(old.var_names)
for ad in new_ads:
    all_genes &= set(ad.var_names)
print(f"Common genes: {len(all_genes)}")

# Subset and concat
print("Concatenating...")
ads = [old[:, list(all_genes)]] + [ad[:, list(all_genes)] for ad in new_ads]
combined = sc.concat(ads, join='inner')
print(f"Combined: {combined.n_obs} cells, {combined.n_vars} genes")

# Save
out_path = f'{OUT}/final_10datasets_merged.h5ad'
combined.write_h5ad(out_path)
print(f"\nSaved: {out_path}")
print(f"Total: {combined.n_obs} cells")

# Per dataset
for ds in sorted(combined.obs['dataset'].unique()):
    n = sum(combined.obs['dataset'] == ds)
    print(f"  {ds}: {n}")
