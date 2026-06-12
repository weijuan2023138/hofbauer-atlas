#!/usr/bin/env python3
"""Spatial validation: remaining 8 samples — HB comm genes on STOMICS"""
import scanpy as sc, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures/FigS")
OUT.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad("/home/weijuan/文档/胎盘单细胞数据/raw_data/UCSF_Li_2026/STOMICS.h5ad", backed='r')

SAMPLES = ['003','005','006','008','009','010','013','016']
mask = (adata.obs['celltype'] == 'HB') & (adata.obs['sample_id'].isin(SAMPLES))
adata_hb = adata[mask].to_memory()
print(f"HB subset: {adata_hb.shape}")

sc.pp.normalize_total(adata_hb, target_sum=1e4)
sc.pp.log1p(adata_hb)

genes = ['SPP1','FN1','CD44','PTPRM','CD47']
spatial = adata_hb.obsm['spatial']

n_samples = len(SAMPLES)
n_genes = len(genes)

fig, axes = plt.subplots(n_samples, n_genes, figsize=(n_genes*4.2, n_samples*4))

for si, sample in enumerate(SAMPLES):
    s_mask = adata_hb.obs['sample_id'] == sample
    xy = spatial[s_mask.values]
    n_hb = s_mask.sum()
    for gi, gene in enumerate(genes):
        ax = axes[si, gi]
        expr = adata_hb[s_mask, gene].X.toarray().flatten()
        vmax = np.percentile(expr, 98)
        scat = ax.scatter(xy[:,0], xy[:,1], c=expr, s=0.3, cmap='RdYlBu_r',
                         vmin=0, vmax=vmax, rasterized=True)
        ax.set_title(f"{gene}  {sample}  (n={n_hb})", fontsize=11, fontweight='bold')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_aspect('equal')
        plt.colorbar(scat, ax=ax, fraction=0.03, pad=0.02)

fig.suptitle("STOMICS Spatial: HB Communication Genes — Samples 003-016",
             fontsize=16, fontweight='bold', y=0.995)
fig.tight_layout()
fig.savefig(OUT / "FigS_stomics_all2.png", dpi=300, bbox_inches='tight')
plt.close()
print(f"Done: {OUT / 'FigS_stomics_all2.png'}")
