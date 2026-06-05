#!/usr/bin/env python3
"""
Re-classify my_preterm_cohort with proper classifier and re-integrate all datasets
Proper criteria: Mac_score > 0 AND DIFF > 0.32
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
print("Re-classify my_preterm_cohort + Re-integrate All Datasets")
print("=" * 60)

# 1. Load my_preterm_cohort
print("\n[1/10] Loading my_preterm_cohort...")
adata_my = sc.read_h5ad(os.path.join(BASE_DIR, 'processed/my_cohort_processed.h5ad'))
print(f"  Total cells: {adata_my.shape[0]}")

# 2. Calculate DIFF and apply proper classifier
print("\n[2/10] Applying proper classifier to my_preterm_cohort...")
adata_my.obs['DIFF'] = adata_my.obs['HBC_score'] - adata_my.obs['Maternal_score']

# Proper criteria: Mac_score > 0 AND DIFF > 0.32
mac_mask = adata_my.obs['Mac_score'] > 0
diff_mask = adata_my.obs['DIFF'] > 0.32
proper_hofbauer = mac_mask & diff_mask

print(f"  Mac_score > 0: {mac_mask.sum()}")
print(f"  DIFF > 0.32: {diff_mask.sum()}")
print(f"  Proper Hofbauer (Mac>0 AND DIFF>0.32): {proper_hofbauer.sum()}")

# 3. Extract proper Hofbauer cells
print("\n[3/10] Extracting proper Hofbauer cells...")
adata_my_hb = adata_my[proper_hofbauer].copy()

# Add metadata
adata_my_hb.obs['dataset'] = 'my_preterm_cohort'
adata_my_hb.obs['disease'] = adata_my_hb.obs['condition']
adata_my_hb.obs['disease_group'] = adata_my_hb.obs['condition'].map({
    'PTL': 'Preterm Labor',
    'PTNL': 'Preterm No Labor',
    'TL': 'Term Labor'
})

print(f"  Proper Hofbauer cells: {adata_my_hb.shape[0]}")
print(f"  Disease group distribution:")
print(adata_my_hb.obs['disease_group'].value_counts())

# 4. Load other datasets (already properly classified)
print("\n[4/10] Loading other datasets...")

# Arutyunyan
adata_arut = sc.read_h5ad(os.path.join(BASE_DIR, 'results/phase2_arutyunyan/arutyunyan_hofbauer.h5ad'))
adata_arut.obs['dataset'] = 'Arutyunyan'
adata_arut.obs['disease'] = 'Normal_1st'
adata_arut.obs['disease_group'] = 'Normal 1st trimester'
print(f"  Arutyunyan: {adata_arut.shape[0]} cells")

# Load combined 8 datasets (without my_preterm_cohort)
adata_8 = sc.read_h5ad(os.path.join(OUTPUT_DIR, 'all_hofbauer_complete.h5ad'))
adata_other = adata_8[adata_8.obs['dataset'] != 'my_preterm_cohort'].copy()
print(f"  Other datasets: {adata_other.shape[0]} cells")

# 5. Find common genes
print("\n[5/10] Finding common genes...")
common_genes = set(adata_other.var_names).intersection(set(adata_arut.var_names)).intersection(set(adata_my_hb.var_names))
print(f"  Common genes: {len(common_genes)}")

# 6. Subset to common genes
print("\n[6/10] Subsetting to common genes...")
adata_other_sub = adata_other[:, list(common_genes)]
adata_arut_sub = adata_arut[:, list(common_genes)]
adata_my_hb_sub = adata_my_hb[:, list(common_genes)]

# 7. Merge all datasets
print("\n[7/10] Merging all datasets...")
adata_all = sc.concat([adata_other_sub, adata_arut_sub, adata_my_hb_sub], join='inner')
print(f"  Total cells: {adata_all.shape[0]}")
print(f"  Dataset distribution:")
print(adata_all.obs['dataset'].value_counts())
print(f"\n  Disease group distribution:")
print(adata_all.obs['disease_group'].value_counts())

# 8. Save merged data
print("\n[8/10] Saving merged data...")
merged_path = os.path.join(OUTPUT_DIR, 'all_hofbauer_final_corrected.h5ad')
adata_all.write_h5ad(merged_path)
print(f"  Saved: {merged_path}")

# 9. Summary
print("\n[9/10] Summary")
print("=" * 60)
print(f"  Total Hofbauer cells: {adata_all.shape[0]}")
print(f"  Datasets: {adata_all.obs['dataset'].nunique()}")
print(f"  Disease groups: {adata_all.obs['disease_group'].nunique()}")

print("\n" + "=" * 60)
print("Done! Ready for clustering.")
print("=" * 60)
