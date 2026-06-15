#!/usr/bin/env python3
"""Spatial dual-track: ECM vs Immune module scores on STOMICS HB spots"""
import scanpy as sc, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

OUT = Path("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/spatial")
OUT.mkdir(parents=True, exist_ok=True)

# Load STOMICS HB spots (backed mode)
adata = sc.read_h5ad("/home/weijuan/文档/胎盘单细胞数据/raw_data/UCSF_Li_2026/STOMICS.h5ad", backed='r')
SAMPLES = ['001','002','004','007','011','012','014','015']
mask = (adata.obs['celltype'] == 'HB') & (adata.obs['sample_id'].isin(SAMPLES))
adata_hb = adata[mask].to_memory()
sc.pp.normalize_total(adata_hb, target_sum=1e4)
sc.pp.log1p(adata_hb)

# ECM module
ecm_genes = ['FN1','SPP1','COL1A2','COL1A1','MMP9','TGFB1','FLT1','PAPPA','AEBP1','NOTUM']
ecm_avail = [g for g in ecm_genes if g in adata_hb.var_names]
ecm_expr = adata_hb[:, ecm_avail].X.toarray()
ecm_z = (ecm_expr - ecm_expr.mean(axis=0)) / (ecm_expr.std(axis=0) + 1e-8)
adata_hb.obs['ECM_score'] = ecm_z.mean(axis=1)

# Immune module
imm_genes = ['HLA-DRA','HLA-DRB1','CD74','FCGR3A','IL1B','TNF','CXCL8','CD44','CD47','IFITM1']
imm_avail = [g for g in imm_genes if g in adata_hb.var_names]
imm_expr = adata_hb[:, imm_avail].X.toarray()
imm_z = (imm_expr - imm_expr.mean(axis=0)) / (imm_expr.std(axis=0) + 1e-8)
adata_hb.obs['Immune_score'] = imm_z.mean(axis=1)

print(f"ECM: {len(ecm_avail)} genes, Immune: {len(imm_avail)} genes")

# Bivariate colormap (blue=ECM only, red=Immune only, purple=both)
def bivariate_cmap():
    c1 = plt.cm.Blues(np.linspace(0.2, 1, 128))[:, :4]
    c2 = plt.cm.Reds(np.linspace(0.2, 1, 128))[:, :4]
    colors = np.zeros((128, 128, 4))
    for i in range(128):
        for j in range(128):
            blend = (c1[i] * (128-j)/128 + c2[j] * i/128)
            colors[i, j] = np.clip(blend, 0, 1)
    return np.clip(colors, 0, 1)

spatial = adata_hb.obsm['spatial']

for sample in SAMPLES:
    s_mask = adata_hb.obs['sample_id'] == sample
    xy = spatial[s_mask.values]
    n_hb = s_mask.sum()
    
    ecm = adata_hb.obs['ECM_score'][s_mask].values
    imm = adata_hb.obs['Immune_score'][s_mask].values
    
    # Cap at 98th percentile
    ecm_cap = np.clip(ecm, 0, np.percentile(ecm, 98))
    imm_cap = np.clip(imm, 0, np.percentile(imm, 98))
    
    # Bin into 128×128 grid for bivariate colormap
    ecm_idx = (ecm_cap / max(ecm_cap.max(), 0.01) * 127).astype(int)
    imm_idx = (imm_cap / max(imm_cap.max(), 0.01) * 127).astype(int)
    cmap_colors = bivariate_cmap()
    point_colors = cmap_colors[imm_idx, ecm_idx]
    
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    plt.subplots_adjust(wspace=0.05)
    
    # Panel 1: ECM only
    ax = axes[0]
    scat = ax.scatter(xy[:,0], xy[:,1], c=ecm_cap, s=0.5, cmap='Blues', rasterized=True)
    ax.set_title('ECM module', fontsize=13, fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect('equal')
    plt.colorbar(scat, ax=ax, fraction=0.03, pad=0.02)
    
    # Panel 2: Immune only
    ax = axes[1]
    scat = ax.scatter(xy[:,0], xy[:,1], c=imm_cap, s=0.5, cmap='Reds', rasterized=True)
    ax.set_title('Immune module', fontsize=13, fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect('equal')
    plt.colorbar(scat, ax=ax, fraction=0.03, pad=0.02)
    
    # Panel 3: Bivariate (Blue=ECM, Red=Immune, Purple=dual)
    ax = axes[2]
    ax.scatter(xy[:,0], xy[:,1], c=point_colors, s=0.5, rasterized=True)
    ax.set_title('ECM (Blue) + Immune (Red)', fontsize=13, fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect('equal')
    
    fig.suptitle(f'STOMICS Dual-Track: Sample {sample}',
                 fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUT / f'stomics_dual_track_{sample}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Dual-track: {sample}")

# Also: correlation between ECM and Immune per sample
print("\nECM-Immune correlation per sample:")
for sample in SAMPLES:
    s_mask = adata_hb.obs['sample_id'] == sample
    r = np.corrcoef(adata_hb.obs['ECM_score'][s_mask], adata_hb.obs['Immune_score'][s_mask])[0,1]
    print(f"  {sample}: r={r:.3f}")

print("\nDual-track spatial analysis complete")
