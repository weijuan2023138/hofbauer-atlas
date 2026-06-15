#!/usr/bin/env python3
"""Extract Hofbauer cells from 9 datasets (+ Vento-Tormo 2018, - GSE183338, + GSE329173) and merge."""
import scanpy as sc
import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = '/home/weijuan/文档/胎盘单细胞数据'
PROCESSED_DIR = os.path.join(BASE_DIR, 'processed')
OUTPUT_DIR = os.path.join(BASE_DIR, 'ucsf_integration/vento_integration')
ORIG_CLASSIFICATION_DIR = os.path.join(BASE_DIR, 'ucsf_integration/results/classification')
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("Extract Hofbauer Cells from 9 Datasets")
print("(- GSE183338, + GSE329173 severe PE)")
print("=" * 60)

datasets = {
    'Vento_Tormo_2018': {
        'path': os.path.join(OUTPUT_DIR, 'vento_tormo_hofbauer.h5ad'),
        'method': 'already_hofbauer',
        'disease': 'Normal_1st',
        'disease_group': 'Normal 1st trimester'
    },
    'Arutyunyan': {
        'path': os.path.join(PROCESSED_DIR, 'arutyunyan_processed.h5ad'),
        'method': 'reclassify',
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
        'path': os.path.join(ORIG_CLASSIFICATION_DIR, 'UCSF_Li_2026_hofbauer.h5ad'),
        'method': 'already_hofbauer',
        'disease': 'Normal',
        'disease_group': 'Normal 1st/2nd/Term'
    },
    'gse329173': {
        'path': os.path.join(ORIG_CLASSIFICATION_DIR, 'gse329173_all_hofbauer.h5ad'),
        'method': 'already_hofbauer',
        'disease': 'Severe_PE',
        'disease_group': 'Severe Preeclampsia'
    }
}

all_hofbauer = []
summary = []

for ds_name, ds_info in datasets.items():
    print(f"\n[Processing] {ds_name}...")
    
    try:
        adata = sc.read_h5ad(ds_info['path'])
        print(f"  Loaded: {adata.shape[0]} cells")

        if ds_info['method'] == 'already_hofbauer':
            hb_cells = adata
            n_hb = adata.shape[0]
        elif ds_info['method'] == 'cell_type_new':
            if 'cell_type_new' in adata.obs.columns:
                hb_mask = adata.obs['cell_type_new'] == 'Hofbauer'
                hb_cells = adata[hb_mask].copy()
                n_hb = hb_mask.sum()
            else:
                print(f"  WARNING: cell_type_new not found, skipping")
                continue
        elif ds_info['method'] == 'cell_type_fine':
            if 'cell_type_fine' in adata.obs.columns:
                hb_mask = adata.obs['cell_type_fine'] == 'Hofbauer'
                hb_cells = adata[hb_mask].copy()
                n_hb = hb_mask.sum()
            else:
                print(f"  WARNING: cell_type_fine not found, skipping")
                continue
        elif ds_info['method'] == 'reclassify':
            hb_mask = adata.obs['cell_type_fine'] == 'Hofbauer'
            hb_cells = adata[hb_mask].copy()
            n_hb = hb_mask.sum()
        else:
            print(f"  WARNING: Unknown method {ds_info['method']}, skipping")
            continue

        print(f"  Hofbauer cells: {n_hb}")
        hb_cells.obs['dataset'] = ds_name
        hb_cells.obs['disease'] = ds_info['disease']
        hb_cells.obs['disease_group'] = ds_info['disease_group']
        
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

# Summary
print("\n" + "=" * 60)
print("Hofbauer Extraction Summary (9 Datasets)")
print("=" * 60)
if summary:
    df = pd.DataFrame(summary)
    print(f"\n{'Dataset':25} | {'Hofbauer':>8} | Disease Group")
    print("-" * 65)
    for _, row in df.iterrows():
        print(f"{row['dataset']:25} | {row['hofbauer']:8} | {row['disease_group']}")
    print(f"\nTotal Hofbauer: {df['hofbauer'].sum()}")

# Merge
print("\n" + "=" * 60)
print("Merging All Hofbauer Cells")
print("=" * 60)

if len(all_hofbauer) > 1:
    common_genes = set(all_hofbauer[0].var_names)
    for ad in all_hofbauer[1:]:
        common_genes = common_genes.intersection(set(ad.var_names))
    print(f"Common genes across all datasets: {len(common_genes)}")

    all_hofbauer_sub = [ad[:, list(common_genes)] for ad in all_hofbauer]
    combined = sc.concat(all_hofbauer_sub, join='inner')
    print(f"Combined: {combined.shape[0]} cells, {combined.shape[1]} genes")

    combined_path = os.path.join(OUTPUT_DIR, 'all_hofbauer_9datasets.h5ad')
    combined.write_h5ad(combined_path)
    print(f"Saved: {combined_path}")

print("\nDone!")
