#!/usr/bin/env python3
"""Re-classify Vento-Tormo 2018 from scratch using full classifier pipeline."""
import scanpy as sc
import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

BASE = '/home/weijuan/文档/胎盘单细胞数据'
CLASSIFIER_DIR = f'{BASE}/results/phase1_classifier'
OUT = f'{BASE}/ucsf_integration/results/classification'

# Load classifier genes
gene_df = pd.read_csv(f'{CLASSIFIER_DIR}/classifier_genes.csv')
hb_genes = gene_df[gene_df['direction']=='HB_up']['gene'].tolist()
mat_genes = gene_df[gene_df['direction']=='MAT_up']['gene'].tolist()
MAC_MARKERS = ['CD68','CD14','AIF1','CSF1R','CD163','ITGAM','FCGR3A','LYZ','C1QA','C1QB']

print("="*60)
print("Vento-Tormo 2018: Full Re-classification")
print("="*60)

# Load full dataset
ad = sc.read_h5ad(f'{BASE}/processed/vento_tormo_2018_processed.h5ad')
print(f"Loaded: {ad.n_obs} cells, {ad.n_vars} genes")

# Check if already log-normalized, normalize if not
if ad.X.max() > 100:
    print("Normalizing...")
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

# Macrophage pre-screen
mac_found = [g for g in MAC_MARKERS if g in ad.var_names]
sc.tl.score_genes(ad, gene_list=mac_found, score_name='Mac_score')
mac_mask = ad.obs['Mac_score'] > 0
print(f"Mac_score > 0: {mac_mask.sum()}/{ad.n_obs}")

# Classifier
hb_found = [g for g in hb_genes if g in ad.var_names]
mat_found = [g for g in mat_genes if g in ad.var_names]
print(f"Classifier genes: HB={len(hb_found)}, MAT={len(mat_found)}")

sc.tl.score_genes(ad, gene_list=hb_found, score_name='HBC_score')
sc.tl.score_genes(ad, gene_list=mat_found, score_name='MAT_score')
ad.obs['DIFF'] = ad.obs['HBC_score'] - ad.obs['MAT_score']

THRESH = 0.32
hbc_mask = mac_mask & (ad.obs['DIFF'] > THRESH)
mat_mask = mac_mask & (ad.obs['DIFF'] < -THRESH)

print(f"DIFF > 0.32: {(ad.obs['DIFF'] > 0.32).sum()}")
print(f"DIFF > 0.32 & Mac>0: {hbc_mask.sum()}")
print(f"DIFF < -0.32 & Mac>0: {mat_mask.sum()}")

# Extract Hofbauer
hb = ad[hbc_mask].copy()
print(f"\nVento-Tormo Hofbauer (classifier): {hb.n_obs} cells")

# Check marker expression
for g in ['FOLR2','CD163','DAB2','CGA','KRT8','CSH1','PECAM1','IGKC']:
    if g in hb.var_names:
        expr = hb[:, g].X.toarray().flatten() if hasattr(hb[:, g].X, 'toarray') else hb[:, g].X.flatten()
        print(f"  {g}: mean={expr.mean():.3f}, P95={np.percentile(expr,95):.3f}")

# Also try stricter threshold
for extra_thresh in [0.5, 0.8, 1.0]:
    n = (ad.obs['DIFF'] > extra_thresh).sum()
    pct = n/ad.n_obs*100
    print(f"  DIFF > {extra_thresh}: {n} ({pct:.1f}%)")

# Compare vs original cell_type_fine HB
orig_hb = ad[ad.obs['cell_type_fine'] == 'Hofbauer']
print(f"\nOriginal cell_type_fine 'Hofbauer': {orig_hb.n_obs}")
overlap = (ad.obs['cell_type_fine'] == 'Hofbauer') & hbc_mask
print(f"Overlap (both original HB AND classifier HB): {overlap.sum()}")
print(f"Original HB NOT passing classifier: {(ad.obs['cell_type_fine']=='Hofbauer').sum() - overlap.sum()}")
print(f"Classifier HB NOT in original: {hbc_mask.sum() - overlap.sum()}")

# Save
hb.write_h5ad(f'{OUT}/vento_tormo_reclassified_hofbauer.h5ad')
print(f"\nSaved: {OUT}/vento_tormo_reclassified_hofbauer.h5ad")
