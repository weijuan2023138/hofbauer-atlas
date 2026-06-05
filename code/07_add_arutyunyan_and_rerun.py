#!/usr/bin/env python3
"""
Add Arutyunyan to the integration and re-run Harmony
Input: all_hofbauer_complete.h5ad + arutyunyan_hofbauer.h5ad
Output: all_hofbauer_final.h5ad (all 9 datasets)
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

# Paths
BASE_DIR = '/home/weijuan/文档/胎盘单细胞数据'
OUTPUT_DIR = os.path.join(BASE_DIR, 'ucsf_integration/results/classification')
FIGURES_DIR = os.path.join(BASE_DIR, 'ucsf_integration/figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

print("=" * 60)
print("Add Arutyunyan and Re-run Harmony (All 9 Datasets)")
print("=" * 60)

# 1. Load existing combined data (8 datasets)
print("\n[1/8] Loading existing combined data (8 datasets)...")
adata_8 = sc.read_h5ad(os.path.join(OUTPUT_DIR, 'all_hofbauer_complete.h5ad'))
print(f"  8 datasets: {adata_8.shape[0]} cells")

# 2. Load Arutyunyan Hofbauer cells
print("\n[2/8] Loading Arutyunyan Hofbauer cells...")
adata_arut = sc.read_h5ad(os.path.join(BASE_DIR, 'results/phase2_arutyunyan/arutyunyan_hofbauer.h5ad'))
print(f"  Arutyunyan: {adata_arut.shape[0]} cells")

# 3. Add metadata to Arutyunyan
print("\n[3/8] Adding metadata to Arutyunyan...")
adata_arut.obs['dataset'] = 'Arutyunyan'
adata_arut.obs['disease'] = 'Normal_1st'
adata_arut.obs['disease_group'] = 'Normal 1st trimester'

# 4. Find common genes
print("\n[4/8] Finding common genes...")
common_genes = set(adata_8.var_names).intersection(set(adata_arut.var_names))
print(f"  Common genes: {len(common_genes)}")

# 5. Subset to common genes
print("\n[5/8] Subsetting to common genes...")
adata_8_sub = adata_8[:, list(common_genes)]
adata_arut_sub = adata_arut[:, list(common_genes)]

# 6. Merge all datasets
print("\n[6/8] Merging all 9 datasets...")
adata_all = sc.concat([adata_8_sub, adata_arut_sub], join='inner')
print(f"  Total cells: {adata_all.shape[0]}")
print(f"  Dataset distribution:")
print(adata_all.obs['dataset'].value_counts())

# 7. Save merged data
print("\n[7/8] Saving merged data...")
merged_path = os.path.join(OUTPUT_DIR, 'all_hofbauer_final.h5ad')
adata_all.write_h5ad(merged_path)
print(f"  Saved: {merged_path}")

# 8. Run Harmony correction
print("\n[8/8] Running Harmony correction...")
import harmonypy as hm

# Normalize and find variable features
sc.pp.normalize_total(adata_all, target_sum=1e4)
sc.pp.log1p(adata_all)
sc.pp.highly_variable_genes(adata_all, n_top_genes=2000, batch_key='dataset')

# PCA
sc.tl.pca(adata_all, n_comps=50)

# Harmony
ho = hm.run_harmony(adata_all.obsm['X_pca'][:, :30], adata_all.obs, 'dataset')
adata_all.obsm['X_pca_harmony'] = ho.Z_corr

# UMAP
sc.pp.neighbors(adata_all, n_neighbors=15, use_rep='X_pca_harmony')
sc.tl.umap(adata_all)

# 9. Visualize
print("\n[Visualization] Creating batch effect plots...")

# Plot 1: UMAP colored by dataset
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# By dataset
sc.pl.umap(adata_all, color='dataset', ax=axes[0], show=False, title='Hofbauer Cells by Dataset (All 9 Datasets)')
axes[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

# By disease_group
sc.pl.umap(adata_all, color='disease_group', ax=axes[1], show=False, title='Hofbauer Cells by Disease Group')
axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'batch_effect_umap_final.png'), dpi=300, bbox_inches='tight')
print(f"  Saved: {FIGURES_DIR}/batch_effect_umap_final.png")

# 10. Statistical assessment
print("\n" + "=" * 60)
print("Batch Effect Statistics (All 9 Datasets)")
print("=" * 60)

from sklearn.metrics import silhouette_score

# Sample a subset for speed
n_sample = min(5000, adata_all.shape[0])
idx = np.random.choice(adata_all.shape[0], n_sample, replace=False)
adata_sub = adata_all[idx]

# Silhouette score for dataset
sil_dataset = silhouette_score(adata_sub.obsm['X_pca_harmony'][:, :30], 
                                adata_sub.obs['dataset'].astype('category').cat.codes)
print(f"\nSilhouette score (dataset): {sil_dataset:.3f}")
print(f"  Interpretation: {'Low batch effect' if sil_dataset < 0.3 else 'Moderate batch effect' if sil_dataset < 0.5 else 'High batch effect'}")

# 11. Summary
print("\n" + "=" * 60)
print("Final Summary")
print("=" * 60)
print(f"  Total Hofbauer cells: {adata_all.shape[0]}")
print(f"  Datasets: {adata_all.obs['dataset'].nunique()}")
print(f"  Common genes: {adata_all.shape[1]}")
print(f"  Batch effect after Harmony (silhouette): {sil_dataset:.3f}")

print("\n" + "=" * 60)
print("Done! All 9 datasets integrated.")
print("=" * 60)
