#!/usr/bin/env python3
"""Extract Hofbauer cells from 10 datasets (-GSE183338, +GSE329173, +GSE298602) and merge."""
import scanpy as sc
import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = '/home/weijuan/文档/胎盘单细胞数据'
PROCESSED_DIR = os.path.join(BASE_DIR, 'processed')
OUTPUT_DIR = os.path.join(BASE_DIR, 'ucsf_integration/vento_integration')
CLASS_DIR = os.path.join(BASE_DIR, 'ucsf_integration/results/classification')
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("Extract Hofbauer Cells from 10 Datasets")
print("(-GSE183338, +GSE329173, +GSE298602)")
print("=" * 60)

datasets = {
    'Vento_Tormo_2018': {
        'path': os.path.join(OUTPUT_DIR, 'vento_tormo_hofbauer.h5ad'),
        'method': 'already_hofbauer',
        'disease': 'Normal_1st', 'disease_group': 'Normal 1st trimester'
    },
    'Arutyunyan': {
        'path': os.path.join(PROCESSED_DIR, 'arutyunyan_processed.h5ad'),
        'method': 'reclassify',
        'disease': 'Normal_1st', 'disease_group': 'Normal 1st trimester'
    },
    'GSE290578': {
        'path': os.path.join(BASE_DIR, 'results/phase3_gse290578/gse290578_reclassified.h5ad'),
        'method': 'cell_type_new',
        'disease': 'Normal/PE', 'disease_group': 'Normal 3rd trimester / Preeclampsia'
    },
    'gse214607': {
        'path': os.path.join(PROCESSED_DIR, 'gse214607_reclassified.h5ad'),
        'method': 'cell_type_new',
        'disease': 'RM/NC', 'disease_group': 'Miscarriage / Normal'
    },
    'hoo_2024': {
        'path': os.path.join(PROCESSED_DIR, 'hoo_2024_reclassified.h5ad'),
        'method': 'cell_type_new',
        'disease': 'Normal/Listeria/Toxoplasma/Malaria', 'disease_group': 'Infection'
    },
    'gse173193': {
        'path': os.path.join(PROCESSED_DIR, 'gse173193_reclassified.h5ad'),
        'method': 'cell_type_new',
        'disease': 'PE', 'disease_group': 'Preeclampsia'
    },
    'gse298119': {
        'path': os.path.join(PROCESSED_DIR, 'gse298119_reclassified.h5ad'),
        'method': 'cell_type_new',
        'disease': 'PE', 'disease_group': 'Preeclampsia'
    },
    'my_preterm_cohort': {
        'path': os.path.join(PROCESSED_DIR, 'my_cohort_processed.h5ad'),
        'method': 'cell_type_fine',
        'disease': 'PTL/TL', 'disease_group': 'Preterm Labor / Term Labor'
    },
    'UCSF_Li_2026': {
        'path': os.path.join(CLASS_DIR, 'UCSF_Li_2026_hofbauer.h5ad'),
        'method': 'already_hofbauer',
        'disease': 'Normal', 'disease_group': 'Normal 1st/2nd/Term'
    },
    'gse329173': {
        'path': os.path.join(CLASS_DIR, 'gse329173_all_hofbauer.h5ad'),
        'method': 'already_hofbauer',
        'disease': 'Severe_PE', 'disease_group': 'Severe Preeclampsia'
    },
    'gse298602': {
        'path': os.path.join(CLASS_DIR, 'gse298602_all_hofbauer.h5ad'),
        'method': 'already_hofbauer',
        'disease': 'PE/Control', 'disease_group': 'Preeclampsia / Control'
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
            hb_cells = adata; n_hb = adata.shape[0]
        elif ds_info['method'] == 'cell_type_new':
            hb_mask = adata.obs['cell_type_new'] == 'Hofbauer'
            hb_cells = adata[hb_mask].copy(); n_hb = hb_mask.sum()
        elif ds_info['method'] == 'cell_type_fine':
            hb_mask = adata.obs['cell_type_fine'] == 'Hofbauer'
            hb_cells = adata[hb_mask].copy(); n_hb = hb_mask.sum()
        elif ds_info['method'] == 'reclassify':
            hb_mask = adata.obs['cell_type_fine'] == 'Hofbauer'
            hb_cells = adata[hb_mask].copy(); n_hb = hb_mask.sum()
        else:
            continue

        print(f"  Hofbauer: {n_hb}")
        hb_cells.obs['dataset'] = ds_name
        hb_cells.obs['disease'] = ds_info['disease']
        hb_cells.obs['disease_group'] = ds_info['disease_group']
        summary.append({'dataset': ds_name, 'hofbauer': n_hb})
        all_hofbauer.append(hb_cells)
    except Exception as e:
        print(f"  ERROR: {e}")

# Summary
print("\n" + "=" * 60)
df = pd.DataFrame(summary)
for _, row in df.iterrows():
    print(f"  {row['dataset']:25s} | {row['hofbauer']:6d}")
print(f"  {'TOTAL':25s} | {df['hofbauer'].sum():6d}")

# Merge
print("\nMerging...")
common_genes = set(all_hofbauer[0].var_names)
for ad in all_hofbauer[1:]:
    common_genes &= set(ad.var_names)
print(f"Common genes: {len(common_genes)}")

# Make obs_names unique by prepending dataset name + dedup
for ad in all_hofbauer:
    ds = ad.obs['dataset'].iloc[0]
    ad.obs_names = [f"{ds}_{b}" for b in ad.obs_names]
    ad.obs_names_make_unique()

all_sub = [ad[:, list(common_genes)] for ad in all_hofbauer]
combined = sc.concat(all_sub, join='inner')
combined.obs_names_make_unique()
print(f"Combined: {combined.shape[0]} cells, {combined.shape[1]} genes")

out = os.path.join(OUTPUT_DIR, 'all_hofbauer_10datasets.h5ad')
combined.write_h5ad(out)
print(f"Saved: {out}")
print("Done!")
