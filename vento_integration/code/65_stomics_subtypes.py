#!/usr/bin/env python3
"""STOMICS spatial: Hofbauer subtype module scores on HB spots"""
import scanpy as sc, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/spatial")
OUT.mkdir(parents=True, exist_ok=True)

# Load STOMICS in backed mode (memory efficient)
adata = sc.read_h5ad("/home/weijuan/文档/胎盘单细胞数据/raw_data/UCSF_Li_2026/STOMICS.h5ad", backed='r')

# Filter to HB spots
SAMPLES = ['001','002','004','007','011','012','014','015']
mask = (adata.obs['celltype'] == 'HB') & (adata.obs['sample_id'].isin(SAMPLES))
adata_hb = adata[mask].to_memory()
print(f"HB spots: {adata_hb.shape}")

sc.pp.normalize_total(adata_hb, target_sum=1e4)
sc.pp.log1p(adata_hb)

# Subtype marker genes (top 20 per subtype from scRNA-seq)
subtype_markers = {
    'Pro-inflammatory': ['CCL3','CCL4','TNF','CXCL8','IL1B','EGR3','RELB','JUNB','NFKB1','IER3','NFKBIA','IER2','FOS','FOSB','JUN','DUSP1','DUSP2','ZFP36','NR4A1','AREG'],
    'MHCII+ Antigen-presenting': ['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1','CD74','FCGR3A','IFITM1','IFITM3','IFI27','IRF8','CIITA','CST3','LGMN','CTSD','CTSB','CTSS','LAMP1','FCGR2B'],
    'Homeostatic': ['FOLR2','CD163','DAB2','MAF','F13A1','MRC1','STAB1','SIGLEC1','VSIG4','CD5L','SEPP1','C1QA','C1QB','C1QC','APOE','APOC1','TREM2','GPNMB','LGALS3','CTSL'],
    'PRKN+ Autophagy': ['SQSTM1','BNIP3','PRKN','LC3B','GABARAP','ATG5','ATG7','BECN1','ULK1','OPTN','PINK1','TFEB','LAMP2','CTSB','CTSD','HSPA8','HSP90AA1','BAG3','VCP','UBB'],
    'Vascular remodeling': ['FN1','COL1A2','COL1A1','SPP1','MMP9','TGFB1','ITGB1','CD44','FLT1','PAPPA','AEBP1','NOTUM','SERPINE2','THBS1','VEGFA','BMP2','LOX','CTGF','FSTL1','IGFBP3'],
    'MKI67+ Proliferating': ['MKI67','TOP2A','CENPF','BIRC5','PCNA','MCM3','MCM5','UBE2C','CCNB1','CDK1','AURKA','AURKB','CKS1B','STMN1','TUBB','H2AFZ','HMGB2','NUSAP1','CKS2','TYMS']
}

# Compute module scores (mean z-score of available genes)
for st, genes in subtype_markers.items():
    avail = [g for g in genes if g in adata_hb.var_names]
    if len(avail) < 5:
        print(f"  {st}: only {len(avail)} genes, skipping")
        continue
    expr = adata_hb[:, avail].X.toarray()
    z = (expr - expr.mean(axis=0)) / (expr.std(axis=0) + 1e-8)
    adata_hb.obs[f'module_{st}'] = z.mean(axis=1)
    print(f"  {st}: {len(avail)} genes")

# Plot: one figure per sample, 3×2 subtype grid
spatial = adata_hb.obsm['spatial']

for sample in SAMPLES:
    s_mask = adata_hb.obs['sample_id'] == sample
    xy = spatial[s_mask.values]
    n_hb = s_mask.sum()
    
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    for i, st in enumerate(subtype_markers.keys()):
        ax = axes[i//3, i%3]
        col = f'module_{st}'
        if col not in adata_hb.obs.columns: continue
        vals = adata_hb.obs[col][s_mask].values
        vmax = np.percentile(vals, 98)
        scat = ax.scatter(xy[:,0], xy[:,1], c=vals, s=0.5, cmap='RdYlBu_r',
                         vmin=0, vmax=max(vmax, 0.1), rasterized=True)
        ax.set_title(f'{st}', fontsize=12, fontweight='bold')
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect('equal')
        plt.colorbar(scat, ax=ax, fraction=0.03, pad=0.02)
    
    fig.suptitle(f'STOMICS: Hofbauer Subtype Modules — Sample {sample} (n={n_hb} HB spots)',
                 fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUT / f'stomics_subtypes_{sample}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plotted: {sample}")

print("\nSTOMICS subtype mapping complete")
