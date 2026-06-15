#!/usr/bin/env python3
"""Extract per-cell disease labels from original classified h5ad files."""
import scanpy as sc
import pandas as pd
import numpy as np

BASE = '/home/weijuan/文档/胎盘单细胞数据'
PROCESSED = f'{BASE}/processed'
CLASS = f'{BASE}/ucsf_integration/results/classification'

# For each dataset with mixed disease, extract per-cell labels
datasets = {
    'GSE290578': {
        'path': f'{BASE}/results/phase3_gse290578/gse290578_reclassified.h5ad',
        'col': 'cell_type_new', 'label': 'Hofbauer',
        'disease_col': 'condition'
    },
    'GSE298602': {
        'path': f'{CLASS}/gse298602_all_hofbauer.h5ad',
        'col': None, 'label': None,
        'disease_col': 'disease'
    },
    'E-MTAB-12795': {
        'path': f'{PROCESSED}/hoo_2024_reclassified.h5ad',
        'col': 'cell_type_new', 'label': 'Hofbauer',
        'disease_col': 'condition'
    },
    'GSE333257': {
        'path': f'{PROCESSED}/my_cohort_processed.h5ad',
        'col': 'cell_type_fine', 'label': 'Hofbauer',
        'disease_col': 'condition'
    },
    'GSE214607': {
        'path': f'{PROCESSED}/gse214607_reclassified.h5ad',
        'col': 'cell_type_new', 'label': 'Hofbauer',
        'disease_col': 'condition'
    },
}

all_labels = {}

for ds_name, info in datasets.items():
    print(f"\n{ds_name}:")
    ad = sc.read_h5ad(info['path'])
    
    # Filter to Hofbauer if needed
    if info['col'] is not None:
        mask = ad.obs[info['col']] == info['label']
        ad = ad[mask]
    
    # Get disease labels
    if info['disease_col'] in ad.obs.columns:
        labels = ad.obs[info['disease_col']].astype(str).values
    else:
        # Try alternate columns
        for col in ['disease', 'disease_group', 'sample']:
            if col in ad.obs.columns:
                labels = ad.obs[col].astype(str).values
                break
    
    counts = pd.Series(labels).value_counts()
    print(f"  n={len(labels)}, groups: {counts.to_dict()}")
    all_labels[ds_name] = labels

# Save as CSV for R to load
import json
with open(f'{BASE}/ucsf_integration/vento_integration/per_cell_disease.json', 'w') as f:
    # Convert numpy arrays to lists
    out = {k: v.tolist() for k, v in all_labels.items()}
    json.dump(out, f)

print(f"\nSaved per-cell disease labels")
