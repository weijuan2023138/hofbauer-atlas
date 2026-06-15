#!/usr/bin/env python3
"""Extract Vento-Tormo 2018 Hofbauer cells and save as h5ad."""
import scanpy as sc
import os

INPUT = "/home/weijuan/文档/胎盘单细胞数据/processed/vento_tormo_2018_processed.h5ad"
OUTPUT_DIR = "/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration"

os.makedirs(OUTPUT_DIR, exist_ok=True)

adata = sc.read_h5ad(INPUT)
print(f"Total cells: {adata.n_obs}")

# Use cell_type_fine annotation
hb = adata[adata.obs['cell_type_fine'] == 'Hofbauer'].copy()
print(f"Hofbauer cells: {hb.n_obs}")

# Save
out = os.path.join(OUTPUT_DIR, "vento_tormo_hofbauer.h5ad")
hb.write(out)
print(f"Saved: {out}")

# Quick stats
print(f"\nGestational weeks range: {hb.obs['Fetus'].unique()[:10]}")
print(f"Conditions: {hb.obs['condition'].value_counts().to_dict()}")
print(f"Locations: {hb.obs['location'].value_counts().to_dict()}")
