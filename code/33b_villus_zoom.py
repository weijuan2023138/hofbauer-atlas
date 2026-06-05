#!/usr/bin/env python3
"""Unified villus zoom plots for 010 and 011"""
import anndata, numpy as np, matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

a = anndata.read_h5ad('/home/weijuan/文档/胎盘单细胞数据/raw_data/UCSF_Li_2026/STOMICS.h5ad')

cc = {'HB':'#D73027','SCT':'#4575B4','VCT':'#91BFDB','fVEC':'#FC8D59','FB':'#99D594',
      'DSC':'#FEE08B','Immune':'#E6F598','PV':'#3288BD','iEVT':'#FDAE61',
      'EVT':'#ABDDA4','mVEC':'#9E0142'}

for sid, gw in [('010','GW20+3'), ('011','GW20+5')]:
    sub = a[a.obs['sample_id']==sid].copy()
    spatial = sub.obsm['spatial']
    fb = spatial[sub.obs['celltype']=='FB']
    cx, cy = fb.mean(axis=0)
    w = 3000
    inw = (spatial[:,0]>cx-w)&(spatial[:,0]<cx+w)&(spatial[:,1]>cy-w)&(spatial[:,1]<cy+w)
    sz, ct = spatial[inw], sub.obs['celltype'].values[inw]
    
    fig, ax = plt.subplots(figsize=(8,8), facecolor='white')
    ax.set_facecolor('white')
    for c in np.unique(ct):
        m = ct==c; s = 2 if c=='HB' else 0.8
        ax.scatter(sz[m,0],sz[m,1],s=s,c=cc.get(c,'grey'),alpha=0.7,rasterized=True)
    # SCT boundary
    sct_z = sz[ct=='SCT']
    if len(sct_z)>10:
        try:
            hull=ConvexHull(sct_z)
            for si in hull.simplices: ax.plot(sct_z[si,0],sct_z[si,1],'b--',lw=1.5,alpha=0.5)
        except: pass
    # Unified annotations
    ax.annotate('Villus Stroma\n(FB + HB + fVEC)', xy=(cx,cy), fontsize=11, ha='center', va='center',
                color='#333', bbox=dict(boxstyle='round',fc='white',alpha=0.8,ec='grey'))
    ax.annotate('SCT/VCT', xy=(cx+1700,cy+1400), fontsize=10, ha='center',
                color='#4575B4', bbox=dict(boxstyle='round',fc='white',alpha=0.8,ec='#4575B4'))
    ax.set_title(f'Sample {sid} ({gw})', fontsize=12, fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([])
    plt.tight_layout()
    plt.savefig(f'/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures/Fig3b_villus_zoom_{sid}.png',
                dpi=250, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'Sample {sid} ({gw}) saved')
