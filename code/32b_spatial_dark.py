#!/usr/bin/env python3
"""All 16 STOMICS slices — dark background version"""
import anndata, matplotlib.pyplot as plt

a = anndata.read_h5ad('/home/weijuan/文档/胎盘单细胞数据/raw_data/UCSF_Li_2026/STOMICS.h5ad')

ct_colors = {
    'HB': '#FF4444', 'SCT': '#6699CC', 'VCT': '#88BBEE',
    'fVEC': '#FF9955', 'FB': '#77CC66', 'DSC': '#EECC44',
    'Immune': '#CCEE88', 'PV': '#66AADD', 'EVT': '#99CC88',
    'iEVT': '#FFAA55', 'Epi': '#9977CC', 'mVEC': '#DD6688',
    'eEVT': '#FF7744', 'EVTpro': '#77CC99', 'DSC4': '#9999CC', 'DSC3': '#DD88BB'
}

OUTDIR = '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'
all_sids = sorted(a.obs['sample_id'].unique(), key=lambda x: int(x))
print(f"Plotting {len(all_sids)} dark slices...")

for sid in all_sids:
    sub = a[a.obs['sample_id']==sid].copy()
    hb_n = (sub.obs['celltype']=='HB').sum()
    hb_pct = hb_n / sub.shape[0] * 100
    spatial = sub.obsm['spatial']
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), facecolor='black')
    
    # Left: all cell types
    ax = axes[0]
    ax.set_facecolor('black')
    for ct in sorted(sub.obs['celltype'].unique()):
        mask = sub.obs['celltype'] == ct
        ax.scatter(spatial[mask,0], spatial[mask,1], s=0.3,
                   c=ct_colors.get(ct,'#666666'), label=ct, alpha=0.6, rasterized=True)
    ax.set_title(f'Sample {sid}  |  {sub.shape[0]:,} cells', fontsize=9, color='white')
    ax.legend(markerscale=6, fontsize=5, loc='center left', bbox_to_anchor=(1.01,0.5),
              frameon=False, title='Cell Type', title_fontsize=6,
              labelcolor='white')
    # Make legend title white
    ax.get_legend().get_title().set_color('white')
    ax.set_xticks([]); ax.set_yticks([])
    ax.tick_params(colors='white')
    
    # Right: HB only
    ax = axes[1]
    ax.set_facecolor('black')
    hb_mask = sub.obs['celltype'] == 'HB'
    ax.scatter(spatial[~hb_mask,0], spatial[~hb_mask,1], s=0.2, c='#1a1a1a', alpha=0.3, rasterized=True)
    ax.scatter(spatial[hb_mask,0], spatial[hb_mask,1], s=1.2, c='#FF4444', alpha=0.8, rasterized=True)
    ax.set_title(f'Hofbauer  |  n={hb_n} ({hb_pct:.1f}%)', fontsize=9, color='white')
    ax.set_xticks([]); ax.set_yticks([])
    
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/Fig3a_spatial_dark_{sid}.png', dpi=150, bbox_inches='tight', facecolor='black')
    plt.close()
    print(f'  {sid}: {sub.shape[0]:,d} cells, HB={hb_n} ({hb_pct:.1f}%)')

print(f"Done — {len(all_sids)} dark slices")
