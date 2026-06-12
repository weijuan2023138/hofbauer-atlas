#!/usr/bin/env python3
"""Spatial validation: 2 pairs × 2 samples × 3 panels — tight layout"""
import scanpy as sc, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde
from pathlib import Path

OUT = Path("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures/FigS")
OUT.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad("/home/weijuan/文档/胎盘单细胞数据/raw_data/UCSF_Li_2026/STOMICS.h5ad", backed='r')

SAMPLES = ['011','012']
mask_full = adata.obs['sample_id'].isin(SAMPLES)
adata_full = adata[mask_full].to_memory()
mask_hb = (adata_full.obs['celltype'] == 'HB').values
adata_hb = adata_full[mask_hb].copy()
print(f"Full: {adata_full.shape}  HB: {adata_hb.shape}")

sc.pp.normalize_total(adata_full, target_sum=1e4)
sc.pp.log1p(adata_full)

pairs = [('SPP1','CD44'), ('FN1','CD47')]

def make_3panel(ax_row, s_full, g_de, g_act, sample):
    """Populate one row of 3 panels"""
    xy_full = s_full.obsm['spatial']
    if sample == '011':
        xy_full = np.column_stack([-xy_full[:,1], xy_full[:,0]])
    is_hb = (s_full.obs['celltype'] == 'HB').values
    xy_hb = xy_full[is_hb]
    expr_de  = s_full[is_hb, g_de].X.toarray().flatten()
    expr_act = s_full[is_hb, g_act].X.toarray().flatten()

    # Panel 1: Tissue + HB expression
    ax = ax_row[0]
    ax.scatter(xy_full[~is_hb,0], xy_full[~is_hb,1], c='lightgrey', s=0.1, rasterized=True)
    vmax = np.percentile(expr_de, 98)
    sc = ax.scatter(xy_hb[:,0], xy_hb[:,1], c=expr_de, s=2, cmap='Reds',
                    vmin=0, vmax=vmax, rasterized=True)
    cbar = plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02, shrink=0.6)
    ax.set_title(f"{sample}  {g_de} expr", fontsize=10, fontweight='bold')
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([]); ax.set_yticks([])

    # Panel 2: Colocalization
    ax = ax_row[1]
    de_n  = (expr_de  - expr_de.min())  / (expr_de.max()  - expr_de.min()  + 1e-10)
    act_n = (expr_act - expr_act.min()) / (expr_act.max() - expr_act.min() + 1e-10)
    colors = np.zeros((len(de_n), 3))
    colors[:,0] = de_n; colors[:,2] = act_n
    ax.scatter(xy_hb[:,0], xy_hb[:,1], c=colors, s=1.2, rasterized=True)
    ax.set_title(f"{sample}  {g_de}×{g_act}", fontsize=10, fontweight='bold')
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([]); ax.set_yticks([])
    ax.legend(handles=[Patch(color='red',label=g_de), Patch(color='blue',label=g_act),
                       Patch(color='purple',label='Both')],
              loc='lower right', fontsize=9, framealpha=0.8)

    # Panel 3: Density + top20%
    ax = ax_row[2]
    kde = gaussian_kde(xy_hb.T, bw_method=0.05)
    xmin, xmax = xy_full[:,0].min(), xy_full[:,0].max()
    ymin, ymax = xy_full[:,1].min(), xy_full[:,1].max()
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    density = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
    ax.imshow(density.T, extent=[xmin,xmax,ymin,ymax], origin='lower',
              cmap='Greys', aspect='auto', alpha=0.7)
    top_de  = expr_de  > np.percentile(expr_de, 80)
    top_act = expr_act > np.percentile(expr_act, 80)
    ax.scatter(xy_hb[top_de,0],  xy_hb[top_de,1],  c='#D73027', s=3, alpha=0.7, label=f'{g_de} top20%')
    ax.scatter(xy_hb[top_act,0], xy_hb[top_act,1], c='#4575B4', s=3, alpha=0.7, label=f'{g_act} top20%')
    ax.set_title(f"{sample}  density+top20%", fontsize=10, fontweight='bold')
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([]); ax.set_yticks([])
    ax.legend(fontsize=9, loc='lower right', framealpha=0.8)

# ---- Individual figures ----
for (g_de, g_act) in pairs:
    fig, axes = plt.subplots(2, 3, figsize=(18, 7))
    for si, sample in enumerate(SAMPLES):
        s_full = adata_full[adata_full.obs['sample_id'] == sample]
        make_3panel(axes[si], s_full, g_de, g_act, sample)

    fig.suptitle(f"STOMICS Spatial: {g_de}×{g_act} — De-repression vs Activation",
                 fontsize=13, fontweight='bold', y=0.99)
    fig.subplots_adjust(left=0.04, right=0.96, top=0.92, bottom=0.04, hspace=0.15, wspace=0.06)
    fig.savefig(OUT / f"FigS_stomics_{g_de}_{g_act}.png", dpi=300, facecolor='white')
    plt.close()
    print(f"Done: FigS_stomics_{g_de}_{g_act}.png")

# ---- Combined 12-panel ----
fig, axes = plt.subplots(4, 3, figsize=(16, 14))
row = 0
for si, sample in enumerate(SAMPLES):
    for (g_de, g_act) in pairs:
        s_full = adata_full[adata_full.obs['sample_id'] == sample]
        make_3panel(axes[row], s_full, g_de, g_act, sample)
        row += 1

fig.suptitle("STOMICS Spatial: De-repression vs Activation — All Pairs",
             fontsize=14, fontweight='bold', y=0.995)
fig.subplots_adjust(left=0.04, right=0.96, top=0.95, bottom=0.03, hspace=0.2, wspace=0.06)
fig.savefig(OUT / "FigS_stomics_spatial.png", dpi=300, facecolor='white')
plt.close()
print("Done: FigS_stomics_spatial.png")
