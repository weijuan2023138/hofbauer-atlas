#!/usr/bin/env python3
"""16 STOMICS slices with gestational week labels"""
import anndata, matplotlib.pyplot as plt

a = anndata.read_h5ad('/home/weijuan/文档/胎盘单细胞数据/raw_data/UCSF_Li_2026/STOMICS.h5ad')

# Gestational weeks from Supplementary Table 4
gw_map = {
    '001': 'GW24+0', '002': 'GW22+3', '003': 'GW22+0', '004': 'GW22+0',
    '005': 'GW21+2', '006': 'GW20+0', '007': 'GW18+6', '008': 'GW21+5',
    '009': 'GW20+1', '010': 'GW20+3', '011': 'GW20+5', '012': 'GW24+0',
    '013': 'GW22+1', '014': 'GW24+1', '015': 'GW20+5', '016': 'GW21+0'
}

ct_colors = {
    'HB': '#D73027', 'SCT': '#4575B4', 'VCT': '#91BFDB',
    'fVEC': '#FC8D59', 'FB': '#99D594', 'DSC': '#FEE08B',
    'Immune': '#E6F598', 'PV': '#3288BD', 'EVT': '#ABDDA4',
    'iEVT': '#FDAE61', 'Epi': '#5E4FA2', 'mVEC': '#9E0142',
    'eEVT': '#F46D43', 'EVTpro': '#66C2A5', 'DSC4': '#8DA0CB', 'DSC3': '#E78AC3'
}

OUTDIR = '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures'
all_sids = sorted(a.obs['sample_id'].unique(), key=lambda x: int(x))

for sid in all_sids:
    sub = a[a.obs['sample_id']==sid].copy()
    hb_n = (sub.obs['celltype']=='HB').sum()
    gw = gw_map.get(sid, '?')
    spatial = sub.obsm['spatial']
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), facecolor='white')
    
    ax = axes[0]
    ax.set_facecolor('white')
    for ct in sorted(sub.obs['celltype'].unique()):
        mask = sub.obs['celltype'] == ct
        ax.scatter(spatial[mask,0], spatial[mask,1], s=0.3,
                   c=ct_colors.get(ct,'grey'), label=ct, alpha=0.5, rasterized=True)
    ax.set_title(f'Sample {sid} ({gw})  |  {sub.shape[0]:,} cells', fontsize=9)
    ax.legend(markerscale=6, fontsize=5, loc='center left', bbox_to_anchor=(1.01,0.5),
              frameon=False, title='Cell Type', title_fontsize=6)
    ax.set_xticks([]); ax.set_yticks([])
    
    ax = axes[1]
    ax.set_facecolor('white')
    hb_mask = sub.obs['celltype'] == 'HB'
    ax.scatter(spatial[~hb_mask,0], spatial[~hb_mask,1], s=0.2, c='lightgrey', alpha=0.2, rasterized=True)
    ax.scatter(spatial[hb_mask,0], spatial[hb_mask,1], s=1.0, c='#D73027', alpha=0.7, rasterized=True)
    ax.set_title(f'Hofbauer  |  n={hb_n} ({hb_n/sub.shape[0]*100:.1f}%)', fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/Fig3a_spatial_{sid}.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  {sid} ({gw}): {sub.shape[0]:,d} cells, HB={hb_n}')

print("Done")
