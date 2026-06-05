#!/usr/bin/env python3
"""
Complete Hofbauer Classifier Pipeline for ALL Datasets
Following the original Phase 4 batch classification approach

Steps:
1. Load data
2. Convert ENSG -> Symbol if needed
3. Macrophage pre-screen (Mac_score > 0)
4. Apply classifier (DIFF > 0.32)
5. Save results + batch effect assessment
"""

import scanpy as sc
import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

# Paths
BASE_DIR = '/home/weijuan/文档/胎盘单细胞数据'
CLASSIFIER_DIR = os.path.join(BASE_DIR, 'results/phase1_classifier')
OUTPUT_DIR = os.path.join(BASE_DIR, 'ucsf_integration/results/classification')
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("Hofbauer Classifier: Complete Pipeline for ALL Datasets")
print("=" * 60)

# 1. Load classifier genes
print("\n[1/10] Loading classifier genes...")
gene_df = pd.read_csv(os.path.join(CLASSIFIER_DIR, 'classifier_genes.csv'))
hb_genes = gene_df[gene_df['direction'] == 'HB_up']['gene'].tolist()
mat_genes = gene_df[gene_df['direction'] == 'MAT_up']['gene'].tolist()
print(f"  HB-up genes: {len(hb_genes)}")
print(f"  MAT-up genes: {len(mat_genes)}")

# Macrophage markers for pre-screen
MAC_MARKERS = ['CD68', 'CD14', 'AIF1', 'CSF1R', 'CD163', 'ITGAM', 'FCGR3A', 'LYZ', 'C1QA', 'C1QB']

# 2. Define datasets
datasets = {
    'Arutyunyan': {
        'path': os.path.join(BASE_DIR, 'raw_data/Arutyunyan_2023/Arutyunyan_2023_spatial_multiomics_placenta.h5ad'),
        'disease': 'Normal_1st',
        'disease_group': 'Normal 1st trimester',
        'gene_id_type': 'ENSG'
    },
    'hoo_2024': {
        'path': os.path.join(BASE_DIR, 'raw_data/hoo_2024/hoo_2024_placenta_pathogens.h5ad'),
        'disease': 'Normal/Listeria/Toxoplasma/Malaria',
        'disease_group': 'Infection',
        'gene_id_type': 'ENSG'
    },
    'UCSF_Li_2026': {
        'path': os.path.join(BASE_DIR, 'raw_data/UCSF_Li_2026/scPlacenta_host.h5ad'),
        'disease': 'Normal',
        'disease_group': 'Normal 1st/2nd/Term',
        'gene_id_type': 'Symbol'
    }
}

# 3. Process each dataset
results = []
all_hofbauer = []

for ds_name, ds_info in datasets.items():
    print(f"\n{'='*60}")
    print(f"[Processing] {ds_name}")
    print(f"{'='*60}")
    
    try:
        # Load data
        adata = sc.read_h5ad(ds_info['path'])
        print(f"  Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")
        
        # Delete .raw if exists
        if hasattr(adata, 'raw') and adata.raw is not None:
            del adata.raw
        
        # Convert ENSG -> Symbol if needed
        if ds_info['gene_id_type'] == 'ENSG':
            print(f"  Converting ENSG IDs to gene symbols...")
            
            # Check for gene_symbols column
            if 'gene_symbols' in adata.var.columns:
                # Direct mapping from var
                ensg_to_symbol = dict(zip(adata.var_names, adata.var['gene_symbols']))
                new_names = []
                seen = {}
                for g in adata.var_names:
                    name = ensg_to_symbol.get(g, str(g))
                    if name in seen:
                        new_names.append(f'{name}__dup{seen[name]}')
                        seen[name] += 1
                    else:
                        new_names.append(name)
                        seen[name] = 1
                adata.var_names = new_names
                print(f"  Converted {len(new_names)} genes")
            else:
                # Use Arutyunyan as reference for mapping
                print(f"  Using Arutyunyan as reference for ENSG->Symbol mapping...")
                arut_path = os.path.join(BASE_DIR, 'raw_data/Arutyunyan_2023/Arutyunyan_2023_spatial_multiomics_placenta.h5ad')
                arut_raw = sc.read_h5ad(arut_path, backed='r')
                
                ensg_map = {}
                for i, g in enumerate(arut_raw.var_names):
                    s = arut_raw.var['gene_symbols'].iloc[i]
                    if isinstance(s, str) and s and s != 'nan':
                        ensg_map[g] = s
                
                new_names = []
                seen = {}
                for g in adata.var_names:
                    name = ensg_map.get(g, str(g))
                    if name in seen:
                        new_names.append(f'{name}__dup{seen[name]}')
                        seen[name] += 1
                    else:
                        new_names.append(name)
                        seen[name] = 1
                adata.var_names = new_names
                print(f"  Converted {len(new_names)} genes")
        
        # Normalize if needed
        if adata.X.max() > 100:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        
        # Macrophage pre-screen
        mac_found = [g for g in MAC_MARKERS if g in adata.var_names]
        print(f"  Mac markers found: {len(mac_found)}/{len(MAC_MARKERS)}")
        
        if len(mac_found) >= 5:
            sc.tl.score_genes(adata, gene_list=mac_found, score_name='Mac_score')
            mac_mask = adata.obs['Mac_score'] > 0
            print(f"  Mac_score > 0: {mac_mask.sum()} cells")
        else:
            mac_mask = np.ones(adata.n_obs, dtype=bool)
            print(f"  WARNING: too few mac markers, using all cells")
        
        # Find classifier genes
        hb_found = [g for g in hb_genes if g in adata.var_names]
        mat_found = [g for g in mat_genes if g in adata.var_names]
        print(f"  Classifier genes: HB={len(hb_found)}/{len(hb_genes)}, MAT={len(mat_found)}/{len(mat_genes)}")
        
        if len(hb_found) < 50 or len(mat_found) < 50:
            print(f"  WARNING: too few classifier genes, skipping")
            continue
        
        # Score and classify
        sc.tl.score_genes(adata, gene_list=hb_found, score_name='HBC_score')
        sc.tl.score_genes(adata, gene_list=mat_found, score_name='MAT_score')
        adata.obs['DIFF'] = adata.obs['HBC_score'] - adata.obs['MAT_score']
        
        OPT_THRESHOLD = 0.32
        hbc_mask = mac_mask & (adata.obs['DIFF'] > OPT_THRESHOLD)
        mat_mask = mac_mask & (adata.obs['DIFF'] < -OPT_THRESHOLD)
        
        # Add classification
        adata.obs['cell_type_new'] = 'Other'
        adata.obs.loc[hbc_mask, 'cell_type_new'] = 'Hofbauer'
        adata.obs.loc[mat_mask, 'cell_type_new'] = 'Maternal_Macrophage'
        adata.obs.loc[mac_mask & ~hbc_mask & ~mat_mask, 'cell_type_new'] = 'Ambiguous'
        
        # Add metadata
        adata.obs['dataset'] = ds_name
        adata.obs['disease'] = ds_info['disease']
        adata.obs['disease_group'] = ds_info['disease_group']
        
        n_hbc = hbc_mask.sum()
        n_mat = mat_mask.sum()
        n_amb = (mac_mask & ~hbc_mask & ~mat_mask).sum()
        
        print(f"  Classification results:")
        print(f"    Hofbauer: {n_hbc}")
        print(f"    Maternal Mac: {n_mat}")
        print(f"    Ambiguous: {n_amb}")
        print(f"    Other: {(~mac_mask).sum()}")
        
        # Validation - check marker expression
        print(f"\n  Validation (marker expression in Hofbauer vs Maternal):")
        for g in ['FOLR2', 'CD163', 'DAB2', 'HLA-DRA']:
            if g in adata.var_names:
                hbc_expr = adata[hbc_mask, g].X.mean() if n_hbc > 0 else 0
                mat_expr = adata[mat_mask, g].X.mean() if n_mat > 0 else 0
                print(f"    {g}: HBC={hbc_expr:.3f}, MAT={mat_expr:.3f}")
        
        # Save results
        results.append({
            'dataset': ds_name,
            'total_cells': len(adata),
            'hofbauer': n_hbc,
            'maternal_mac': n_mat,
            'ambiguous': n_amb,
            'other': (~mac_mask).sum(),
            'pct_hofbauer': n_hbc / len(adata) * 100,
            'mean_diff': adata.obs['DIFF'].mean(),
            'mean_diff_hb': adata.obs.loc[hbc_mask, 'DIFF'].mean() if n_hbc > 0 else 0,
            'mean_diff_mat': adata.obs.loc[mat_mask, 'DIFF'].mean() if n_mat > 0 else 0
        })
        
        # Extract Hofbauer cells
        hb_cells = adata[hbc_mask].copy()
        
        # Save Hofbauer cells
        output_path = os.path.join(OUTPUT_DIR, f'{ds_name}_hofbauer.h5ad')
        hb_cells.write_h5ad(output_path)
        print(f"\n  Saved: {output_path}")
        
        all_hofbauer.append(hb_cells)
        
    except Exception as e:
        print(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()
        continue

# 4. Summary
print("\n" + "=" * 60)
print("Classification Summary")
print("=" * 60)

if results:
    df = pd.DataFrame(results)
    print("\nDataset | Total | Hofbauer | Maternal | Ambiguous | % HB")
    print("-" * 70)
    for _, row in df.iterrows():
        print(f"{row['dataset']:15} | {row['total_cells']:6} | {row['hofbauer']:8} | {row['maternal_mac']:8} | {row['ambiguous']:9} | {row['pct_hofbauer']:5.1f}%")
    
    print(f"\nTotal Hofbauer: {df['hofbauer'].sum()}")
    print(f"Total cells: {df['total_cells'].sum()}")
    
    # Save summary
    df.to_csv(os.path.join(OUTPUT_DIR, 'classification_summary.csv'), index=False)
    print(f"\nSaved: {OUTPUT_DIR}/classification_summary.csv")

# 5. Batch effect assessment (if we have multiple datasets)
if len(all_hofbauer) > 1:
    print("\n" + "=" * 60)
    print("Batch Effect Assessment")
    print("=" * 60)
    
    # Merge all Hofbauer cells
    print("\n[5/10] Merging all Hofbauer cells...")
    
    # Find common genes across all datasets
    common_genes = set(all_hofbauer[0].var_names)
    for ad in all_hofbauer[1:]:
        common_genes = common_genes.intersection(set(ad.var_names))
    
    print(f"  Common genes across all datasets: {len(common_genes)}")
    
    # Subset to common genes
    all_hofbauer_sub = [ad[:, list(common_genes)] for ad in all_hofbauer]
    
    # Merge
    combined = sc.concat(all_hofbauer_sub, join='inner')
    print(f"  Combined: {combined.shape[0]} cells")
    
    # Save combined
    combined_path = os.path.join(OUTPUT_DIR, 'all_hofbauer_combined.h5ad')
    combined.write_h5ad(combined_path)
    print(f"  Saved: {combined_path}")

print("\n" + "=" * 60)
print("Done!")
print("=" * 60)
