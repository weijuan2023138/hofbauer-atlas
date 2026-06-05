#!/usr/bin/env python3
"""
Assess batch effects across all Hofbauer cells
Input: all_hofbauer_combined.h5ad
Output: Batch effect visualization + statistics
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
INPUT_DIR = os.path.join(BASE_DIR, 'ucsf_integration/results/classification')
OUTPUT_DIR = os.path.join(BASE_DIR, 'ucsf_integration/figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("Batch Effect Assessment: Hofbauer Cells")
print("=" * 60)

# 1. Load combined data
print("\n[1/5] Loading combined Hofbauer cells...")
adata = sc.read_h5ad(os.path.join(INPUT_DIR, 'all_hofbauer_combined.h5ad'))
print(f"  Total cells: {adata.shape[0]}")
print(f"  Genes: {adata.shape[1]}")

# 2. Dataset distribution
print("\n[2/5] Dataset distribution:")
print(adata.obs['dataset'].value_counts())

# 3. Normalize and find variable features
print("\n[3/5] Normalizing and finding variable features...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='dataset')
print(f"  Highly variable genes: {adata.var['highly_variable'].sum()}")

# 4. PCA
print("\n[4/5] Running PCA...")
sc.tl.pca(adata, n_comps=50)

# 5. Compute neighbors and UMAP
print("\n[5/5] Computing neighbors and UMAP...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)

# 6. Visualize batch effects
print("\n[Visualization] Creating batch effect plots...")

# Plot 1: UMAP colored by dataset
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# By dataset
sc.pl.umap(adata, color='dataset', ax=axes[0], show=False, title='Hofbauer Cells by Dataset')
axes[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

# By disease_group
sc.pl.umap(adata, color='disease_group', ax=axes[1], show=False, title='Hofbauer Cells by Disease Group')
axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'batch_effect_umap.png'), dpi=300, bbox_inches='tight')
print(f"  Saved: {OUTPUT_DIR}/batch_effect_umap.png")

# Plot 2: Key marker expression
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

markers = ['FOLR2', 'CD163', 'DAB2', 'HLA-DRA', 'C1QA', 'SPP1']
for i, gene in enumerate(markers):
    ax = axes[i // 3, i % 3]
    if gene in adata.var_names:
        sc.pl.umap(adata, color=gene, ax=ax, show=False, title=gene, vmax='p99')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'batch_effect_markers.png'), dpi=300, bbox_inches='tight')
print(f"  Saved: {OUTPUT_DIR}/batch_effect_markers.png")

# 7. Statistical assessment
print("\n" + "=" * 60)
print("Batch Effect Statistics")
print("=" * 60)

# Calculate silhouette score (dataset vs biological)
from sklearn.metrics import silhouette_score

# Sample a subset for speed
n_sample = min(5000, adata.shape[0])
idx = np.random.choice(adata.shape[0], n_sample, replace=False)
adata_sub = adata[idx]

# Silhouette score for dataset (should be low if batch effect is small)
sil_dataset = silhouette_score(adata_sub.obsm['X_pca'][:, :30], 
                                adata_sub.obs['dataset'].astype('category').cat.codes)
print(f"\nSilhouette score (dataset): {sil_dataset:.3f}")
print(f"  Interpretation: {'Low batch effect' if sil_dataset < 0.3 else 'Moderate batch effect' if sil_dataset < 0.5 else 'High batch effect'}")

# Check if datasets cluster separately
print("\nDataset clustering analysis:")
for ds in adata.obs['dataset'].unique():
    mask = adata.obs['dataset'] == ds
    ds_pca = adata[mask].obsm['X_pca'][:, :30].mean(axis=0)
    print(f"  {ds}: PCA center = [{ds_pca[0]:.2f}, {ds_pca[1]:.2f}, {ds_pca[2]:.2f}]")

# 8. Summary
print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print(f"  Total Hofbauer cells: {adata.shape[0]}")
print(f"  Datasets: {adata.obs['dataset'].nunique()}")
print(f"  Common genes: {adata.shape[1]}")
print(f"  Batch effect (silhouette): {sil_dataset:.3f}")

print("\n" + "=" * 60)
print("Done!")
print("=" * 60)
