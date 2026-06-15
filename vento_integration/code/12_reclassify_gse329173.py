#!/usr/bin/env python3
"""Re-classify GSE329173 from raw 10X data with careful ENSG→Symbol conversion."""
import scanpy as sc
import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

BASE = '/home/weijuan/文档/胎盘单细胞数据'
RAW = f'{BASE}/raw_data/gse329173'
CLASSIFIER_DIR = f'{BASE}/results/phase1_classifier'
OUT = f'{BASE}/ucsf_integration/results/classification'
os.makedirs(OUT, exist_ok=True)

gene_df = pd.read_csv(f'{CLASSIFIER_DIR}/classifier_genes.csv')
hb_genes = gene_df[gene_df['direction']=='HB_up']['gene'].tolist()
mat_genes = gene_df[gene_df['direction']=='MAT_up']['gene'].tolist()
MAC_MARKERS = ['CD68','CD14','AIF1','CSF1R','CD163','ITGAM','FCGR3A','LYZ','C1QA','C1QB']

# ENSG→Symbol from Arutyunyan
print("Loading ENSG→Symbol reference...")
arut = sc.read_h5ad(f'{BASE}/raw_data/Arutyunyan_2023/Arutyunyan_2023_spatial_multiomics_placenta.h5ad', backed='r')
ensg_map = {}
for i, g in enumerate(arut.var_names):
    s = arut.var['gene_symbols'].iloc[i]
    if isinstance(s, str) and s and s != 'nan':
        ensg_map[g] = s
print(f"  {len(ensg_map)} mappings")

samples = {'SPE1': 'GSM9698580', 'SPE2': 'GSM9698581', 'SPE3': 'GSM9698582'}
results = []

for name, gsm in samples.items():
    print(f"\n{'='*60}")
    print(f"Re-classifying {name}...")
    print(f"{'='*60}")
    
    # Load raw 10X
    ad = sc.read_10x_mtx(RAW, var_names='gene_ids', prefix=f'{gsm}_{name}_')
    print(f"  Raw: {ad.n_obs} cells, {ad.n_vars} genes")
    
    # ENSG→Symbol
    new_names = []
    seen = {}
    for g in ad.var_names:
        sym = ensg_map.get(g, g)
        if sym in seen:
            new_names.append(f'{sym}__dup{seen[sym]}')
            seen[sym] += 1
        else:
            new_names.append(sym)
            seen[sym] = 1
    ad.var_names = new_names
    print(f"  After ENSG→Symbol: {ad.n_vars} genes")
    
    # Check key markers exist
    for g in ['FOLR2','CD163','DAB2','CGA','KRT8','HLA-DRA']:
        if g in ad.var_names:
            print(f"    {g}: present")
        else:
            print(f"    {g}: MISSING!")
    
    # Normalize
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    
    # Mac pre-screen
    mac_found = [g for g in MAC_MARKERS if g in ad.var_names]
    print(f"  Mac markers found: {len(mac_found)}/{len(MAC_MARKERS)}")
    
    if len(mac_found) >= 5:
        sc.tl.score_genes(ad, gene_list=mac_found, score_name='Mac_score')
        mac_mask = ad.obs['Mac_score'] > 0
    else:
        mac_mask = np.ones(ad.n_obs, dtype=bool)
    print(f"  Mac_score > 0: {mac_mask.sum()}/{ad.n_obs}")
    
    # Classifier
    hb_found = [g for g in hb_genes if g in ad.var_names]
    mat_found = [g for g in mat_genes if g in ad.var_names]
    print(f"  Classifier: HB={len(hb_found)}/{len(hb_genes)}, MAT={len(mat_found)}/{len(mat_genes)}")
    
    sc.tl.score_genes(ad, gene_list=hb_found, score_name='HBC_score')
    sc.tl.score_genes(ad, gene_list=mat_found, score_name='MAT_score')
    ad.obs['DIFF'] = ad.obs['HBC_score'] - ad.obs['MAT_score']
    
    # Strict Hofbauer
    hbc_mask = mac_mask & (ad.obs['DIFF'] > 0.32)
    print(f"  Hofbauer (DIFF>0.32): {hbc_mask.sum()}")
    
    # Also check DIFF distribution
    for t in [0, 0.2, 0.32, 0.5, 0.8]:
        n = (mac_mask & (ad.obs['DIFF'] > t)).sum()
        print(f"    DIFF > {t}: {n}")
    
    # Check CGA in Hofbauer
    if hbc_mask.sum() > 0 and 'CGA' in ad.var_names:
        cga_vals = ad[hbc_mask, 'CGA'].X.toarray().flatten()
        print(f"  CGA in Hofbauer: mean={cga_vals.mean():.3f}, P95={np.percentile(cga_vals,95):.3f}")
    
    # Save reclassified
    ad.obs['dataset'] = 'gse329173'
    ad.obs['sample'] = name
    out_path = f'{OUT}/gse329173_{name}_reclassified_v2.h5ad'
    ad.write_h5ad(out_path)
    
    results.append({'sample': name, 'total': ad.n_obs, 'hofbauer': hbc_mask.sum()})

# Summary
print(f"\n{'='*60}")
print("GSE329173 Re-classification Summary")
print(f"{'='*60}")
df = pd.DataFrame(results)
print(df.to_string(index=False))
print(f"\nTotal Hofbauer: {df['hofbauer'].sum()}")
