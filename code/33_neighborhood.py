#!/usr/bin/env python3
"""Neighborhood analysis: what cell types surround Hofbauer cells"""
import anndata
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

print("Loading STOMICS...")
a = anndata.read_h5ad('/home/weijuan/文档/胎盘单细胞数据/raw_data/UCSF_Li_2026/STOMICS.h5ad')

# Use sample 007 (largest, 180K spots, 15K HB)
sub = a[a.obs['sample_id']=='007'].copy()
print(f"Sample 007: {sub.shape[0]} cells")

coords = sub.obsm['spatial']
celltypes = sub.obs['celltype'].values
hb_idx = np.where(celltypes == 'HB')[0]

# Build KD-tree and find neighbors
tree = cKDTree(coords)
radius = 50  # 50 pixels ~ 25 µm
neighbor_counts = {}
hb_neighbors_all = []

for i in hb_idx:
    neighbors = tree.query_ball_point(coords[i], radius)
    neighbor_types = celltypes[neighbors]
    # Exclude self (HB)
    neighbor_types = neighbor_types[neighbor_types != 'HB']
    for ct in neighbor_types:
        neighbor_counts[ct] = neighbor_counts.get(ct, 0) + 1

# Also compute background expectation
bg_counts = {}
for ct in celltypes:
    if ct != 'HB':
        bg_counts[ct] = np.sum(celltypes == ct)
total_bg = sum(bg_counts.values())

# Enrichment: observed/expected
total_obs = sum(neighbor_counts.values())
enrichment = {}
for ct in neighbor_counts:
    obs_freq = neighbor_counts[ct] / total_obs
    exp_freq = bg_counts[ct] / total_bg
    enrichment[ct] = obs_freq / exp_freq

# Sort by enrichment
sorted_cts = sorted(enrichment.items(), key=lambda x: -x[1])

print("\nCell types enriched around Hofbauer:")
for ct, enr in sorted_cts:
    print(f"  {ct:12s}  enrichment={enr:.2f}  observed={neighbor_counts.get(ct,0)}")

# Plot
fig, ax = plt.subplots(figsize=(6, 4))
cts = [x[0] for x in sorted_cts]
enrs = [x[1] for x in sorted_cts]
colors = ['#D73027' if e > 1.2 else '#4575B4' if e < 0.8 else 'grey' for e in enrs]
ax.barh(range(len(cts)), enrs, color=colors, height=0.6)
ax.set_yticks(range(len(cts)))
ax.set_yticklabels(cts, fontsize=8)
ax.axvline(1, color='black', linewidth=0.5, linestyle='--')
ax.set_xlabel('Enrichment (observed/expected)')
ax.set_title('Cell types surrounding Hofbauer cells')
plt.tight_layout()

OUTDIR = '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'
plt.savefig(f'{OUTDIR}/Fig3b_neighborhood.png', dpi=300, bbox_inches='tight')
print(f"\nSaved Fig3b_neighborhood.png")
