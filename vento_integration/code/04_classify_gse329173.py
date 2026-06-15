#!/usr/bin/env python3
"""Classify Hofbauer cells in GSE329173 (3 severe PE samples)."""
import scanpy as sc
import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = '/home/weijuan/文档/胎盘单细胞数据'
CLASSIFIER_DIR = os.path.join(BASE_DIR, 'results/phase1_classifier')
RAW_DIR = os.path.join(BASE_DIR, 'raw_data/gse329173')
OUTPUT_DIR = os.path.join(BASE_DIR, 'ucsf_integration/results/classification')
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("GSE329173: Hofbauer Classification (3 severe PE)")
print("=" * 60)

# Load classifier genes
gene_df = pd.read_csv(os.path.join(CLASSIFIER_DIR, 'classifier_genes.csv'))
hb_genes = gene_df[gene_df['direction'] == 'HB_up']['gene'].tolist()
mat_genes = gene_df[gene_df['direction'] == 'MAT_up']['gene'].tolist()
print(f"Classifier: {len(hb_genes)} HB-up, {len(mat_genes)} MAT-up")

MAC_MARKERS = ['CD68', 'CD14', 'AIF1', 'CSF1R', 'CD163', 'ITGAM', 'FCGR3A', 'LYZ', 'C1QA', 'C1QB']

# Arutyunyan as ENSG→Symbol reference
print("\nLoading ENSG→Symbol reference from Arutyunyan...")
arut = sc.read_h5ad(os.path.join(BASE_DIR, 'raw_data/Arutyunyan_2023/Arutyunyan_2023_spatial_multiomics_placenta.h5ad'), backed='r')
ensg_map = {}
for i, g in enumerate(arut.var_names):
    s = arut.var['gene_symbols'].iloc[i]
    if isinstance(s, str) and s and s != 'nan':
        ensg_map[g] = s
print(f"  {len(ensg_map)} ENSG→Symbol mappings")

# Process each sample
samples = ['SPE1', 'SPE2', 'SPE3']
results = []

for sample in samples:
    print(f"\n{'='*60}")
    print(f"[Processing] {sample}")
    print(f"{'='*60}")
    
    try:
        # Load 10X data
        mtx_path = os.path.join(RAW_DIR, f'GSM969858{0 if sample=="SPE1" else 1 if sample=="SPE2" else 2}_{sample}_matrix.mtx.gz')
        bc_path = os.path.join(RAW_DIR, f'GSM969858{0 if sample=="SPE1" else 1 if sample=="SPE2" else 2}_{sample}_barcodes.tsv.gz')
        feat_path = os.path.join(RAW_DIR, f'GSM969858{0 if sample=="SPE1" else 1 if sample=="SPE2" else 2}_{sample}_features.tsv.gz')
        
        adata = sc.read_10x_mtx(RAW_DIR, 
                                var_names='gene_ids',
                                prefix=f'GSM969858{0 if sample=="SPE1" else 1 if sample=="SPE2" else 2}_{sample}_')
        print(f"  Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")
        
        # Convert ENSG → Symbol
        new_names = []
        seen = {}
        for g in adata.var_names:
            name = ensg_map.get(g, g)
            if name in seen:
                new_names.append(f'{name}__dup{seen[name]}')
                seen[name] += 1
            else:
                new_names.append(name)
                seen[name] = 1
        adata.var_names = new_names
        print(f"  ENSG→Symbol: {len(new_names)} genes")
        
        # Normalize
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Macrophage pre-screen
        mac_found = [g for g in MAC_MARKERS if g in adata.var_names]
        sc.tl.score_genes(adata, gene_list=mac_found, score_name='Mac_score')
        mac_mask = adata.obs['Mac_score'] > 0
        print(f"  Mac_score > 0: {mac_mask.sum()} / {adata.n_obs} cells")
        
        # Classifier scoring
        hb_found = [g for g in hb_genes if g in adata.var_names]
        mat_found = [g for g in mat_genes if g in adata.var_names]
        print(f"  Classifier genes found: HB={len(hb_found)}/{len(hb_genes)}, MAT={len(mat_found)}/{len(mat_genes)}")
        
        sc.tl.score_genes(adata, gene_list=hb_found, score_name='HBC_score')
        sc.tl.score_genes(adata, gene_list=mat_found, score_name='MAT_score')
        adata.obs['DIFF'] = adata.obs['HBC_score'] - adata.obs['MAT_score']
        
        OPT_THRESHOLD = 0.32
        hbc_mask = mac_mask & (adata.obs['DIFF'] > OPT_THRESHOLD)
        mat_mask = mac_mask & (adata.obs['DIFF'] < -OPT_THRESHOLD)
        
        adata.obs['cell_type_new'] = 'Other'
        adata.obs.loc[hbc_mask, 'cell_type_new'] = 'Hofbauer'
        adata.obs.loc[mat_mask, 'cell_type_new'] = 'Maternal_Macrophage'
        adata.obs.loc[mac_mask & ~hbc_mask & ~mat_mask, 'cell_type_new'] = 'Ambiguous'
        
        # Metadata
        adata.obs['dataset'] = f'gse329173_{sample}'
        adata.obs['disease'] = 'Severe_PE'
        adata.obs['disease_group'] = 'Severe Preeclampsia'
        
        n_hbc = hbc_mask.sum()
        n_mat = mat_mask.sum()
        n_amb = (mac_mask & ~hbc_mask & ~mat_mask).sum()
        
        print(f"  Classification:")
        print(f"    Hofbauer: {n_hbc}")
        print(f"    Maternal Mac: {n_mat}")
        print(f"    Ambiguous: {n_amb}")
        print(f"    Other: {(~mac_mask).sum()}")
        
        # Validation markers
        for g in ['FOLR2', 'CD163', 'DAB2', 'HLA-DRA']:
            if g in adata.var_names:
                hbc_expr = adata[hbc_mask, g].X.toarray().mean() if n_hbc > 0 else 0
                mat_expr = adata[mat_mask, g].X.toarray().mean() if n_mat > 0 else 0
                print(f"    {g}: HBC={hbc_expr:.3f}, MAT={mat_expr:.3f}")
        
        # Save
        output_path = os.path.join(OUTPUT_DIR, f'gse329173_{sample}_reclassified.h5ad')
        adata.write_h5ad(output_path)
        print(f"  Saved: {output_path}")
        
        # Extract Hofbauer only
        hb_cells = adata[hbc_mask].copy()
        hb_path = os.path.join(OUTPUT_DIR, f'gse329173_{sample}_hofbauer.h5ad')
        hb_cells.write_h5ad(hb_path)
        print(f"  Hofbauer saved: {hb_path}")
        
        results.append({
            'sample': sample,
            'total_cells': adata.n_obs,
            'hofbauer': n_hbc,
            'maternal_mac': n_mat,
            'ambiguous': n_amb,
            'other': (~mac_mask).sum(),
            'pct_hofbauer': n_hbc / adata.n_obs * 100 if adata.n_obs > 0 else 0
        })
        
    except Exception as e:
        print(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()
        continue

# Summary
print("\n" + "=" * 60)
print("GSE329173 Classification Summary")
print("=" * 60)
if results:
    df = pd.DataFrame(results)
    print(df.to_string(index=False))
    print(f"\nTotal Hofbauer: {df['hofbauer'].sum()}")
    
    # Merge all Hofbauer
    all_hb = []
    for sample in samples:
        hb_path = os.path.join(OUTPUT_DIR, f'gse329173_{sample}_hofbauer.h5ad')
        if os.path.exists(hb_path):
            all_hb.append(sc.read_h5ad(hb_path))
    
    if len(all_hb) > 1:
        common_genes = set(all_hb[0].var_names)
        for ad in all_hb[1:]:
            common_genes = common_genes.intersection(set(ad.var_names))
        all_hb_sub = [ad[:, list(common_genes)] for ad in all_hb]
        combined = sc.concat(all_hb_sub, join='inner')
        combined_path = os.path.join(OUTPUT_DIR, 'gse329173_all_hofbauer.h5ad')
        combined.write_h5ad(combined_path)
        print(f"Merged Hofbauer: {combined.shape[0]} cells → {combined_path}")

print("\nDone!")
