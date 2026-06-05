#!/usr/bin/env python3
"""
Fix my_preterm_cohort metadata and re-integrate all datasets
Correct conditions: PTL, PTNL, TL (not just PTL/TL)
"""

import scanpy as sc
import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

# Paths
BASE_DIR = '/home/weijuan/文档/胎盘单细胞数据'
OUTPUT_DIR = os.path.join(BASE_DIR, 'ucsf_integration/results/classification')
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("Fix my_preterm_cohort Metadata + Re-integrate")
print("=" * 60)

# 1. Load all Hofbauer cells (8 datasets without Arutyunyan)
print("\n[1/10] Loading existing combined data (8 datasets)...")
adata_8 = sc.read_h5ad(os.path.join(OUTPUT_DIR, 'all_hofbauer_complete.h5ad'))
print(f"  8 datasets: {adata_8.shape[0]} cells")

# 2. Load Arutyunyan Hofbauer cells
print("\n[2/10] Loading Arutyunyan Hofbauer cells...")
adata_arut = sc.read_h5ad(os.path.join(BASE_DIR, 'results/phase2_arutyunyan/arutyunyan_hofbauer.h5ad'))
print(f"  Arutyunyan: {adata_arut.shape[0]} cells")

# 3. Add metadata to Arutyunyan
print("\n[3/10] Adding metadata to Arutyunyan...")
adata_arut.obs['dataset'] = 'Arutyunyan'
adata_arut.obs['disease'] = 'Normal_1st'
adata_arut.obs['disease_group'] = 'Normal 1st trimester'

# 4. Fix my_preterm_cohort metadata
print("\n[4/10] Fixing my_preterm_cohort metadata...")

# Load my_preterm_cohort Hofbauer cells
adata_my = sc.read_h5ad(os.path.join(BASE_DIR, 'processed/my_cohort_processed.h5ad'))
hb_mask = adata_my.obs['cell_type_fine'] == 'Hofbauer'
adata_my_hb = adata_my[hb_mask].copy()

# Add correct disease_group based on condition
def get_disease_group(condition):
    if condition == 'PTL':
        return 'Preterm Labor'
    elif condition == 'PTNL':
        return 'Preterm No Labor'
    elif condition == 'TL':
        return 'Term Labor'
    else:
        return 'Unknown'

adata_my_hb.obs['disease_group'] = adata_my_hb.obs['condition'].apply(get_disease_group)
adata_my_hb.obs['dataset'] = 'my_preterm_cohort'
adata_my_hb.obs['disease'] = adata_my_hb.obs['condition']

print(f"  my_preterm_cohort Hofbauer: {adata_my_hb.shape[0]} cells")
print(f"  disease_group distribution:")
print(adata_my_hb.obs['disease_group'].value_counts())

# 5. Remove old my_preterm_cohort from 8-dataset combined
print("\n[5/10] Removing old my_preterm_cohort from combined data...")
adata_other = adata_8[adata_8.obs['dataset'] != 'my_preterm_cohort'].copy()
print(f"  Other datasets: {adata_other.shape[0]} cells")

# 6. Find common genes
print("\n[6/10] Finding common genes...")
common_genes = set(adata_other.var_names).intersection(set(adata_arut.var_names)).intersection(set(adata_my_hb.var_names))
print(f"  Common genes: {len(common_genes)}")

# 7. Subset to common genes
print("\n[7/10] Subsetting to common genes...")
adata_other_sub = adata_other[:, list(common_genes)]
adata_arut_sub = adata_arut[:, list(common_genes)]
adata_my_hb_sub = adata_my_hb[:, list(common_genes)]

# 8. Merge all datasets
print("\n[8/10] Merging all 9 datasets...")
adata_all = sc.concat([adata_other_sub, adata_arut_sub, adata_my_hb_sub], join='inner')
print(f"  Total cells: {adata_all.shape[0]}")
print(f"  Dataset distribution:")
print(adata_all.obs['dataset'].value_counts())
print(f"\n  Disease group distribution:")
print(adata_all.obs['disease_group'].value_counts())

# 9. Save merged data
print("\n[9/10] Saving merged data...")
merged_path = os.path.join(OUTPUT_DIR, 'all_hofbauer_final_fixed.h5ad')
adata_all.write_h5ad(merged_path)
print(f"  Saved: {merged_path}")

# 10. Summary
print("\n[10/10] Summary")
print("=" * 60)
print(f"  Total Hofbauer cells: {adata_all.shape[0]}")
print(f"  Datasets: {adata_all.obs['dataset'].nunique()}")
print(f"  Disease groups: {adata_all.obs['disease_group'].nunique()}")

print("\n" + "=" * 60)
print("Done! All datasets integrated with correct metadata.")
print("=" * 60)
