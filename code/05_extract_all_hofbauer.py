#!/usr/bin/env python3
"""
Extract Hofbauer cells from ALL datasets and merge
Input: All reclassified h5ad files
Output: all_hofbauer_complete.h5ad (all 9 datasets)
"""

import scanpy as sc
import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

# Paths
BASE_DIR = '/home/weijuan/文档/胎盘单细胞数据'
PROCESSED_DIR = os.path.join(BASE_DIR, 'processed')
OUTPUT_DIR = os.path.join(BASE_DIR, 'ucsf_integration/results/classification')
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("Extract Hofbauer Cells from ALL Datasets")
print("=" * 60)

# 1. Define all datasets with their Hofbauer extraction method
datasets = {
    'Arutyunyan': {
        'path': os.path.join(PROCESSED_DIR, 'arutyunyan_processed.h5ad'),
        'method': 'reclassify',  # Need to reclassify
        'disease': 'Normal_1st',
        'disease_group': 'Normal 1st trimester'
    },
    'GSE290578': {
        'path': os.path.join(BASE_DIR, 'results/phase3_gse290578/gse290578_reclassified.h5ad'),
        'method': 'cell_type_new',
        'disease': 'Normal/PE',
        'disease_group': 'Normal 3rd trimester / Preeclampsia'
    },
    'gse214607': {
        'path': os.path.join(PROCESSED_DIR, 'gse214607_reclassified.h5ad'),
        'method': 'cell_type_new',
        'disease': 'RM/NC',
        'disease_group': 'Miscarriage / Normal'
    },
    'hoo_2024': {
        'path': os.path.join(PROCESSED_DIR, 'hoo_2024_reclassified.h5ad'),
        'method': 'cell_type_new',
        'disease': 'Normal/Listeria/Toxoplasma/Malaria',
        'disease_group': 'Infection'
    },
    'gse173193': {
        'path': os.path.join(PROCESSED_DIR, 'gse173193_reclassified.h5ad'),
        'method': 'cell_type_new',
        'disease': 'PE',
        'disease_group': 'Preeclampsia'
    },
    'gse183338': {
        'path': os.path.join(PROCESSED_DIR, 'gse183338_reclassified.h5ad'),
        'method': 'cell_type_new',
        'disease': 'PE',
        'disease_group': 'Preeclampsia'
    },
    'gse298119': {
        'path': os.path.join(PROCESSED_DIR, 'gse298119_reclassified.h5ad'),
        'method': 'cell_type_new',
        'disease': 'PE',
        'disease_group': 'Preeclampsia'
    },
    'my_preterm_cohort': {
        'path': os.path.join(PROCESSED_DIR, 'my_cohort_processed.h5ad'),
        'method': 'cell_type_fine',
        'disease': 'PTL/TL',
        'disease_group': 'Preterm Labor / Term Labor'
    },
    'UCSF_Li_2026': {
        'path': os.path.join(OUTPUT_DIR, 'UCSF_Li_2026_hofbauer.h5ad'),
        'method': 'already_hofbauer',  # Already extracted
        'disease': 'Normal',
        'disease_group': 'Normal 1st/2nd/Term'
    }
}

# 2. Extract Hofbauer cells from each dataset
all_hofbauer = []
summary = []

for ds_name, ds_info in datasets.items():
    print(f"\n[Processing] {ds_name}...")
    
    try:
        # Load data
        adata = sc.read_h5ad(ds_info['path'])
        print(f"  Loaded: {adata.shape[0]} cells")
        
        # Extract Hofbauer cells based on method
        if ds_info['method'] == 'already_hofbauer':
            # Already extracted
            hb_cells = adata
            n_hb = adata.shape[0]
        elif ds_info['method'] == 'cell_type_new':
            # Use cell_type_new column
            if 'cell_type_new' in adata.obs.columns:
                hb_mask = adata.obs['cell_type_new'] == 'Hofbauer'
                hb_cells = adata[hb_mask].copy()
                n_hb = hb_mask.sum()
            else:
                print(f"  WARNING: cell_type_new not found, skipping")
                continue
        elif ds_info['method'] == 'cell_type_fine':
            # Use cell_type_fine column
            if 'cell_type_fine' in adata.obs.columns:
                hb_mask = adata.obs['cell_type_fine'] == 'Hofbauer'
                hb_cells = adata[hb_mask].copy()
                n_hb = hb_mask.sum()
            else:
                print(f"  WARNING: cell_type_fine not found, skipping")
                continue
        else:
            print(f"  WARNING: Unknown method {ds_info['method']}, skipping")
            continue
        
        print(f"  Hofbauer cells: {n_hb}")
        
        # Add metadata
        hb_cells.obs['dataset'] = ds_name
        hb_cells.obs['disease'] = ds_info['disease']
        hb_cells.obs['disease_group'] = ds_info['disease_group']
        
        # Save results
        summary.append({
            'dataset': ds_name,
            'hofbauer': n_hb,
            'disease': ds_info['disease'],
            'disease_group': ds_info['disease_group']
        })
        
        all_hofbauer.append(hb_cells)
        
    except Exception as e:
        print(f"  ERROR: {e}")
        continue

# 3. Summary
print("\n" + "=" * 60)
print("Hofbauer Extraction Summary")
print("=" * 60)

if summary:
    df = pd.DataFrame(summary)
    print("\nDataset | Hofbauer | Disease Group")
    print("-" * 50)
    for _, row in df.iterrows():
        print(f"{row['dataset']:20} | {row['hofbauer']:8} | {row['disease_group']}")
    
    print(f"\nTotal Hofbauer: {df['hofbauer'].sum()}")

# 4. Merge all Hofbauer cells
print("\n" + "=" * 60)
print("Merging All Hofbauer Cells")
print("=" * 60)

if len(all_hofbauer) > 1:
    # Find common genes
    common_genes = set(all_hofbauer[0].var_names)
    for ad in all_hofbauer[1:]:
        common_genes = common_genes.intersection(set(ad.var_names))
    
    print(f"Common genes across all datasets: {len(common_genes)}")
    
    # Subset to common genes
    all_hofbauer_sub = [ad[:, list(common_genes)] for ad in all_hofbauer]
    
    # Merge
    combined = sc.concat(all_hofbauer_sub, join='inner')
    print(f"Combined: {combined.shape[0]} cells")
    
    # Save combined
    combined_path = os.path.join(OUTPUT_DIR, 'all_hofbauer_complete.h5ad')
    combined.write_h5ad(combined_path)
    print(f"Saved: {combined_path}")

print("\n" + "=" * 60)
print("Done!")
print("=" * 60)
