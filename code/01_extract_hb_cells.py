#!/usr/bin/env python3
"""
Step 1: Extract Hofbauer (HB) cells from UCSF Li 2026 dataset
Input: scPlacenta_host.h5ad (191,735 cells, 32,981 genes)
Output: ucsf_hb.h5ad (2,950 HB cells)
"""

import scanpy as sc
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

# Paths
INPUT = '/home/weijuan/文档/胎盘单细胞数据/raw_data/UCSF_Li_2026/scPlacenta_host.h5ad'
OUTPUT = '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/ucsf_hb.h5ad'
META_OUTPUT = '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/results/ucsf_hb_meta.csv'

print("=" * 60)
print("Step 1: Extract Hofbauer cells from UCSF Li 2026")
print("=" * 60)

# 1. Load data in backed mode (file is 3GB)
print("\n[1/5] Loading scPlacenta_host.h5ad in backed mode...")
adata = sc.read_h5ad(INPUT, backed='r')
print(f"  Total cells: {adata.shape[0]}")
print(f"  Total genes: {adata.shape[1]}")

# 2. Check major_class distribution
print("\n[2/5] Cell type distribution (major_class):")
major_counts = adata.obs['major_class'].value_counts()
for ct, count in major_counts.items():
    print(f"  {ct}: {count}")

# 3. Filter HB cells
print("\n[3/5] Filtering Hofbauer (HB) cells...")
mask = adata.obs['major_class'] == 'HB'
n_hb = mask.sum()
print(f"  HB cells found: {n_hb}")

# 4. Convert to memory and save
print("\n[4/5] Converting to memory and saving...")
hb = adata[mask].to_memory()

# Add metadata
print("\n  Metadata columns:")
for col in hb.obs.columns:
    print(f"    {col}: {hb.obs[col].nunique()} unique values")

# Save h5ad
hb.write_h5ad(OUTPUT)
print(f"\n  Saved: {OUTPUT}")

# Save metadata as CSV for R
meta = hb.obs.copy()
meta.to_csv(META_OUTPUT)
print(f"  Saved: {META_OUTPUT}")

# 5. Summary
print("\n[5/5] Summary:")
print(f"  HB cells: {n_hb}")
print(f"  Genes: {hb.shape[1]}")

# Gestational age distribution
print("\n  Gestational age group:")
ga_counts = hb.obs['gestational_age_group'].value_counts()
for ga, count in ga_counts.items():
    print(f"    {ga}: {count}")

# Sample distribution
print("\n  Sample distribution:")
sid_counts = hb.obs['sample_id'].value_counts().sort_index()
for sid, count in sid_counts.items():
    print(f"    {sid}: {count}")

# Gestational week distribution
print("\n  Gestational week:")
gw_counts = hb.obs['gestational_week'].value_counts().sort_index()
for gw, count in gw_counts.items():
    print(f"    {gw}: {count}")

print("\n" + "=" * 60)
print("Done! Next step: 02_merge_and_integrate.R")
print("=" * 60)
