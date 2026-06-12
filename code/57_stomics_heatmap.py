#!/usr/bin/env python3
"""Cell type proportion heatmap — HB row highlighted in red"""
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

print("Loading STOMICS...")
adata = sc.read_h5ad('/home/weijuan/文档/胎盘单细胞数据/raw_data/UCSF_Li_2026/STOMICS.h5ad')

prop = pd.crosstab(adata.obs['sample_id'], adata.obs['celltype'], normalize='index') * 100
top_cts = prop.sum().sort_values(ascending=False).head(12).index
prop = prop[top_cts]
sample_order = sorted(prop.index.tolist())
prop = prop.loc[sample_order]

# Annotation text colors: red for HB row, black for others
annot_colors = pd.DataFrame('black', index=prop.index, columns=prop.columns)
annot_colors.loc[:, 'HB'] = '#C62828'

# Plot with per-cell text colors
fig, ax = plt.subplots(figsize=(6, 5))
sns.heatmap(prop, annot=True, fmt='.1f', cmap='YlOrRd', linewidths=0.5,
            cbar_kws={'label': '% of cells'}, ax=ax,
            annot_kws={'fontsize': 7})
ax.set_title('Cell type composition across STOMICS samples', fontweight='bold', fontsize=12)
ax.set_ylabel('Sample')
ax.set_xlabel('Cell type')

# Overlay red text for HB column
for i, sample in enumerate(prop.index):
    val = prop.loc[sample, 'HB']
    ax.text(list(prop.columns).index('HB') + 0.5, i + 0.5, f'{val:.1f}',
            ha='center', va='center', fontsize=7, fontweight='bold', color='#C62828')

plt.tight_layout()
plt.savefig('figures/Fig3/补充图 Fig3b_stomics_composition.png', dpi=300, bbox_inches='tight')
print("Saved with HB highlighted in red")
