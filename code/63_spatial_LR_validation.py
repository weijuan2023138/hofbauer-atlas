#!/usr/bin/env python3
"""
空间验证：Fig4 CellChat鉴定的关键L-R对在空间转录组中的共定位验证
输出：Fig3/Fig4补充图 → supplement/FigS_spatial_LR_validation.png
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.spatial import KDTree
from matplotlib.colors import LinearSegmentedColormap

# ============================================================
# 1. 加载数据
# ============================================================
ad = sc.read_h5ad('/home/weijuan/文档/胎盘单细胞数据/raw_data/UCSF_Li_2026/STOMICS.h5ad')

# 选3个代表性样本
samples = ['001', '006', '015']
ad_sub = ad[ad.obs['sample_id'].isin(samples)].copy()

# 关键细胞类型
cell_types = ['HB', 'FB', 'fVEC']
ct_colors = {'HB': '#D73027', 'FB': '#FDAE61', 'fVEC': '#4575B4'}

# 关键L-R对
lr_pairs = [
    ('SPP1',  'ITGAV', 'SPP1 → ITGAV',  ['HB','FB'],      ['HB','fVEC','FB']),
    ('FN1',   'ITGB1', 'FN1 → ITGB1',   ['FB','HB','fVEC'], ['HB','fVEC','FB']),
    ('COL1A2','ITGA1', 'COL1A2 → ITGA1',['FB'],            ['HB']),
    ('PTPRM', 'PTPRM', 'PTPRM ↔ PTPRM', ['fVEC'],          ['HB']),
]

figdir = '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures/Fig3/supplement'
import os; os.makedirs(figdir, exist_ok=True)

# ============================================================
# 2. 绘制函数
# ============================================================

def plot_spatial_celltypes(ad_samp, sample_id, ax):
    """Panel a: 三种细胞的空间位置"""
    for ct in cell_types:
        mask = ad_samp.obs['celltype'] == ct
        xy = ad_samp.obsm['spatial'][mask]
        if len(xy) > 5000:
            idx = np.random.choice(len(xy), 5000, replace=False)
            xy = xy[idx]
        ax.scatter(xy[:,0], xy[:,1], c=ct_colors[ct], s=0.5, alpha=0.6, label=ct, rasterized=True)
    ax.set_title(f'Sample {sample_id}\nHB + FB + fVEC', fontsize=9, fontweight='bold')
    ax.legend(markerscale=8, fontsize=7, loc='upper right', frameon=True)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect('equal')

def plot_lr_coexpression(ad_samp, lig_gene, rec_gene, lig_cts, rec_cts, title, ax):
    """Panel b-e: 配体×受体的空间共表达"""
    lig_mask = ad_samp.obs['celltype'].isin(lig_cts)
    rec_mask = ad_samp.obs['celltype'].isin(rec_cts)
    
    lig_expr = ad_samp[lig_mask, lig_gene].X.toarray().flatten() if hasattr(ad_samp.X, 'toarray') else ad_samp[lig_mask, lig_gene].X.flatten()
    rec_expr = ad_samp[rec_mask, rec_gene].X.toarray().flatten() if hasattr(ad_samp.X, 'toarray') else ad_samp[rec_mask, rec_gene].X.flatten()
    
    lig_xy = ad_samp.obsm['spatial'][lig_mask]
    rec_xy = ad_samp.obsm['spatial'][rec_mask]
    
    # Ligand — red intensity
    lig_colors = np.clip(lig_expr / (np.percentile(lig_expr, 95) + 1e-8), 0, 1)
    rgba_lig = np.zeros((len(lig_xy), 4))
    rgba_lig[:, 0] = 1.0  # R
    rgba_lig[:, 3] = lig_colors * 0.8
    
    # Receptor — blue intensity
    rec_colors = np.clip(rec_expr / (np.percentile(rec_expr, 95) + 1e-8), 0, 1)
    rgba_rec = np.zeros((len(rec_xy), 4))
    rgba_rec[:, 2] = 1.0  # B
    rgba_rec[:, 3] = rec_colors * 0.8
    
    # Plot receptor first (background), ligand on top
    ax.scatter(rec_xy[:,0], rec_xy[:,1], c=rgba_rec, s=0.3, rasterized=True)
    ax.scatter(lig_xy[:,0], lig_xy[:,1], c=rgba_lig, s=0.5, rasterized=True)
    
    ax.set_title(title, fontsize=8, fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([])
    
    # Legend
    rp = mpatches.Patch(color='red', alpha=0.6, label=f'{lig_gene} (lig)')
    bp = mpatches.Patch(color='blue', alpha=0.6, label=f'{rec_gene} (rec)')
    ax.legend(handles=[rp, bp], fontsize=6, loc='upper right', frameon=True)

def compute_neighborhood_enrichment(ad_samp, lig_gene, rec_gene, lig_cts, rec_cts, radius=50, n_perm=100):
    """Panel f: 邻域共表达富集分析"""
    lig_mask = ad_samp.obs['celltype'].isin(lig_cts)
    rec_mask = ad_samp.obs['celltype'].isin(rec_cts)
    
    lig_expr = ad_samp[lig_mask, lig_gene].X.toarray().flatten() if hasattr(ad_samp.X, 'toarray') else ad_samp[lig_mask, lig_gene].X.flatten()
    rec_expr = ad_samp[rec_mask, rec_gene].X.toarray().flatten() if hasattr(ad_samp.X, 'toarray') else ad_samp[rec_mask, rec_gene].X.flatten()
    
    lig_xy = ad_samp.obsm['spatial'][lig_mask]
    rec_xy = ad_samp.obsm['spatial'][rec_mask]
    
    lig_high = lig_expr > np.percentile(lig_expr, 50)
    rec_high = rec_expr > np.percentile(rec_expr, 50)
    
    lig_tree = KDTree(lig_xy[lig_high])
    rec_tree = KDTree(rec_xy[rec_high])
    
    # Observed: count lig_high bins within radius of rec_high bins
    counts = lig_tree.query_ball_tree(rec_tree, r=radius)
    observed = np.mean([len(c) for c in counts])
    
    # Permutation
    null_vals = []
    all_xy = ad_samp.obsm['spatial']
    for _ in range(n_perm):
        rand_idx = np.random.choice(len(all_xy), sum(lig_high), replace=False)
        rand_tree = KDTree(all_xy[rand_idx])
        rand_counts = rand_tree.query_ball_tree(rec_tree, r=radius)
        null_vals.append(np.mean([len(c) for c in rand_counts]))
    
    null_vals = np.array(null_vals)
    enrichment = observed / np.mean(null_vals)
    pval = (np.sum(null_vals >= observed) + 1) / (n_perm + 1)
    return enrichment, pval

# ============================================================
# 3. 生成主图
# ============================================================
n_pairs = len(lr_pairs)
fig, axes = plt.subplots(3, 1 + n_pairs, figsize=(3.5 * (1 + n_pairs), 10), 
                          dpi=150, facecolor='white')

for i, (sample_id, ax_row) in enumerate(zip(samples, axes)):
    ad_samp = ad_sub[ad_sub.obs['sample_id'] == sample_id].copy()
    
    # Panel a: cell type map
    plot_spatial_celltypes(ad_samp, sample_id, ax_row[0])
    
    # Panels b-e: L-R co-expression
    for j, (lig, rec, title, lig_cts, rec_cts) in enumerate(lr_pairs):
        plot_lr_coexpression(ad_samp, lig, rec, lig_cts, rec_cts, title, ax_row[1+j])

plt.tight_layout(pad=0.5)
plt.savefig(f'{figdir}/FigS_spatial_LR_maps.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')
plt.close()
print(f'Saved: {figdir}/FigS_spatial_LR_maps.png')

# ============================================================
# 4. 邻域富集量化
# ============================================================
enrichment_results = []
for lig, rec, title, lig_cts, rec_cts in lr_pairs:
    for sample_id in samples:
        ad_samp = ad_sub[ad_sub.obs['sample_id'] == sample_id]
        enrich, pval = compute_neighborhood_enrichment(ad_samp, lig, rec, lig_cts, rec_cts)
        enrichment_results.append({
            'L-R_pair': title, 'Sample': sample_id,
            'Enrichment': enrich, 'p_value': pval
        })

df_enrich = pd.DataFrame(enrichment_results)
df_summary = df_enrich.groupby('L-R_pair')['Enrichment'].agg(['mean','std']).reset_index()

fig2, ax2 = plt.subplots(figsize=(4.5, 3.5), facecolor='white')
colors_bar = ['#D73027', '#4575B4', '#FDAE61', '#66A61E']
bars = ax2.bar(range(len(df_summary)), df_summary['mean'], 
               yerr=df_summary['std'], color=colors_bar[:len(df_summary)],
               edgecolor='black', linewidth=0.5, capsize=4)
ax2.axhline(y=1.0, color='#808080', linestyle='--', linewidth=0.8, label='Random expectation')
ax2.set_xticks(range(len(df_summary)))
ax2.set_xticklabels(df_summary['L-R_pair'], rotation=25, ha='right', fontsize=9, fontweight='bold')
ax2.set_ylabel('Neighborhood enrichment\n(observed / random)', fontsize=11, fontweight='bold')
ax2.set_title('Spatial co-localization of L-R pairs\n(radius=50μm, HB-centric)', 
              fontsize=12, fontweight='bold')
ax2.legend(fontsize=8)

# Add significance stars
for i, (_, row) in enumerate(df_summary.iterrows()):
    pair_pvals = df_enrich[df_enrich['L-R_pair'] == row['L-R_pair']]['p_value']
    sig = sum(pv < 0.05 for pv in pair_pvals)
    if sig >= 2:
        ax2.text(i, row['mean'] + row['std'] + 0.15, '**', ha='center', fontsize=14, fontweight='bold')

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures/Fig4/补充图Fig4E_spatial_LR_enrichment.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()
print(f'Saved: {figdir}/FigS_spatial_LR_enrichment.png')

# Print summary
print('\n=== Neighborhood enrichment summary ===')
print(df_summary.to_string(index=False))
print('\nDone.')
