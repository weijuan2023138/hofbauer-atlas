#!/usr/bin/env python3
"""pySCENIC-based regulon activity and TF-regulon coupling slopes"""
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Load merged h5ad
ad = sc.read_h5ad('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/final_10datasets_merged.h5ad')
print(f"Loaded: {ad.n_obs} cells, {ad.n_vars} genes")

# Disease labels
labels = pd.read_csv('/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/per_cell_disease_labels.csv', index_col=0)
detail = labels['disease_detail']

dc6 = np.full(ad.n_obs, '', dtype=object)
for i, ds in enumerate(ad.obs['dataset']):
    det = detail.iloc[i] if i < len(detail) else ''
    # Normal Early
    if ds in ('E-MTAB-12421','E-MTAB-6701') or (ds=='E-MTAB-12795' and det=='normal'):
        dc6[i] = 'Normal Early'
    elif ds == 'GSE214607':
        dc6[i] = 'Miscarriage'
    elif det in ('toxoplasmosis','listeriosis','Plasmodium malariae malaria'):
        dc6[i] = 'Infection'
    elif (ds=='GSE290578' and det=='Normal') or (ds=='GSE298602' and det=='Control'):
        dc6[i] = 'Normal Late'
    elif det in ('PE','PreE_SF','gHTN','GSE173193','GSE298119') or (ds=='GSE290578' and det=='PE'):
        dc6[i] = 'PE'
    elif det in ('PTL','PTNL'):
        dc6[i] = 'Preterm'

# PE 3-group
dc3 = np.full(ad.n_obs, '', dtype=object)
for i, ds in enumerate(ad.obs['dataset']):
    det = detail.iloc[i] if i < len(detail) else ''
    if (ds=='GSE290578' and det=='Normal') or (ds=='GSE298602' and det=='Control'):
        dc3[i] = 'Normal Late'
    elif ds=='GSE290578' and det=='PE':
        dc3[i] = 'Early PE'
    elif det in ('PreE_SF','gHTN','GSE173193','GSE298119'):
        dc3[i] = 'Late PE'

# Regulons
regulons = {
    'CEBPA': ['SPP1','FN1','COL1A2','PAPPA','FLT1','CD44','AEBP1','NOTUM','COL1A1','IGF1'],
    'IRF1': ['HLA-DRA','CD74','STAT1','IFITM1','CXCL10','GBP2','FCGR3A','IFI27','HLA-DQB1','HLA-DPA1','IRF8'],
    'KLF4': ['TGFB1','VEGFA','THBS1','ITGAV','PTPRM','MMP9','BMP2','COL4A1','CD47','ITGB1','ITGB5']
}
for tf in regulons:
    regulons[tf] = [g for g in regulons[tf] if g in ad.var_names]

# AUCell-like: rank-based regulon activity score
# For each cell, rank genes, compute fraction of regulon genes in top 5% of ranked genes
expr = ad.X.toarray() if hasattr(ad.X, 'toarray') else ad.X
from scipy.stats import rankdata
auc_scores = {}
for tf, genes in regulons.items():
    auc = np.zeros(ad.n_obs)
    for i in range(ad.n_obs):
        ranks = rankdata(-expr[i,:], method='average')  # descending
        top_k = int(ad.n_vars * 0.05)
        hits = sum(1 for g in genes if ad.var_names.get_loc(g) < ad.n_vars and ranks[ad.var_names.get_loc(g)] <= top_k)
        auc[i] = hits / len(genes)
    auc_scores[tf] = auc
    print(f"  {tf}: AUC mean={auc.mean():.3f}")

# Compute slopes
def compute_slopes(groups, target_groups):
    results = []
    for tf in regulons:
        for g in target_groups:
            mask = groups == g
            if mask.sum() < 50: continue
            expr_tf = expr[mask, ad.var_names.get_loc(tf)]
            r, p = spearmanr(expr_tf, auc_scores[tf][mask])
            results.append({'Disease': g, 'TF': tf, 'Slope': r})
    return pd.DataFrame(results)

groups6 = ['Normal Early','Miscarriage','Infection','Normal Late','PE','Preterm']
slope6 = compute_slopes(dc6, groups6)
print("\n=== 6-group slopes ===")
for _, row in slope6.iterrows():
    print(f"  {row['Disease']:12s} {row['TF']:6s} r={row['Slope']:.3f}")

groups3 = ['Normal Late','Early PE','Late PE']
slope3 = compute_slopes(dc3, groups3)
print("\n=== PE 3-group slopes ===")
for _, row in slope3.iterrows():
    print(f"  {row['Disease']:12s} {row['TF']:6s} r={row['Slope']:.3f}")

# Plot
FIGDIR = '/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/vento_integration/figures'
colors = {'CEBPA': '#4A7BB0', 'IRF1': '#D93829', 'KLF4': '#7B3294'}

def plot_slopes(df, groups, fname, title, w=10, h=5, rot=30):
    fig, ax = plt.subplots(figsize=(w/2.5, h/2.5))
    x = np.arange(len(groups))
    w_bar = 0.25
    for i, tf in enumerate(['CEBPA','IRF1','KLF4']):
        vals = [df[(df.Disease==g) & (df.TF==tf)].Slope.values[0] if len(df[(df.Disease==g) & (df.TF==tf)])>0 else 0 for g in groups]
        ax.bar(x + i*w_bar, vals, w_bar, color=colors[tf], label=tf, edgecolor='black', linewidth=0.5)
    ax.axhline(y=0, color='grey', linewidth=0.5)
    ax.set_xticks(x + w_bar)
    ax.set_xticklabels(groups, rotation=rot, ha='right' if rot>0 else 'center', fontsize=11, fontweight='bold')
    ax.set_ylabel('TF-regulon coupling (Spearman r)', fontsize=10)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(fontsize=11, loc='upper right')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/{fname}', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

plot_slopes(slope6, groups6, 'Fig7d_slope_6groups.png', 'Regulatory coupling strength')
plot_slopes(slope3, groups3, 'Fig7d_slope_PE3.png', 'Regulatory coupling — PE subtypes', w=6, rot=0)
print("\nBoth figures saved")
