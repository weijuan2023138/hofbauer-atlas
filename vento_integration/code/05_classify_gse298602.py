#!/usr/bin/env python3
"""Classify Hofbauer cells in GSE298602 (11 samples: 3 Control + 7 PE + 1 gHTN)."""
import scanpy as sc
import numpy as np
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = '/home/weijuan/文档/胎盘单细胞数据'
CLASSIFIER_DIR = os.path.join(BASE_DIR, 'results/phase1_classifier')
RAW_DIR = os.path.join(BASE_DIR, 'raw_data/GSE298602')
OUTPUT_DIR = os.path.join(BASE_DIR, 'ucsf_integration/results/classification')
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("GSE298602: Hofbauer Classification (11 samples)")
print("=" * 60)

# Load classifier genes
gene_df = pd.read_csv(os.path.join(CLASSIFIER_DIR, 'classifier_genes.csv'))
hb_genes = gene_df[gene_df['direction'] == 'HB_up']['gene'].tolist()
mat_genes = gene_df[gene_df['direction'] == 'MAT_up']['gene'].tolist()
print(f"Classifier: {len(hb_genes)} HB-up, {len(mat_genes)} MAT-up")

MAC_MARKERS = ['CD68', 'CD14', 'AIF1', 'CSF1R', 'CD163', 'ITGAM', 'FCGR3A', 'LYZ', 'C1QA', 'C1QB']

# ENSG→Symbol reference from Arutyunyan
print("\nLoading ENSG→Symbol reference...")
arut = sc.read_h5ad(os.path.join(BASE_DIR, 'raw_data/Arutyunyan_2023/Arutyunyan_2023_spatial_multiomics_placenta.h5ad'), backed='r')
ensg_map = {}
for i, g in enumerate(arut.var_names):
    s = arut.var['gene_symbols'].iloc[i]
    if isinstance(s, str) and s and s != 'nan':
        ensg_map[g] = s
print(f"  {len(ensg_map)} ENSG→Symbol mappings")

# Sample metadata from SOFT
sample_info = {
    '3658-OP-1': ('Control', 'Control'),
    '3716-OP-1': ('PreE_SF', 'Preeclampsia'),
    '3723-OP-1': ('gHTN', 'Gestational Hypertension'),
    '3804-OP-1': ('PreE_SF', 'Preeclampsia'),
    '5009-OP-1': ('PreE_SF', 'Preeclampsia'),
    '5065-OP-2': ('Control', 'Control'),
    '5065-OP-3': ('Control', 'Control'),
    '5109-OP-1': ('PreE_SF', 'Preeclampsia'),
    '5300-OP-1': ('PreE_SF', 'Preeclampsia'),
    '5723-OP-2': ('PreE_SF', 'Preeclampsia'),
    '6015-OP-1': ('PreE_SF', 'Preeclampsia'),
}

# Map GSM IDs to samples
gsm_map = {}
for f in os.listdir(RAW_DIR):
    if f.endswith('_raw_feature_bc_matrix.h5'):
        gsm = f.split('_')[0]
        sample_name = f.replace(f'_{gsm}_', '_').replace('_raw_feature_bc_matrix.h5', '')
        # sample_name looks like "3658-OP-1" etc
        # Actually the format is GSM9018282_3658-OP-1_raw_feature_bc_matrix.h5
        # Let's extract differently
        parts = f.split('_')
        sample_name = parts[1]  # e.g., "3658-OP-1"
        gsm_map[sample_name] = f

print(f"\nFound {len(gsm_map)} samples")

results = []

for sample_name in sorted(gsm_map.keys()):
    disease, disease_group = sample_info.get(sample_name, ('Unknown', 'Unknown'))
    h5_file = gsm_map[sample_name]
    h5_path = os.path.join(RAW_DIR, h5_file)
    
    print(f"\n{'='*60}")
    print(f"[Processing] {sample_name} ({disease})")
    print(f"{'='*60}")
    
    try:
        adata = sc.read_10x_h5(h5_path)
        # Make var_names unique
        adata.var_names_make_unique()
        print(f"  Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")
        
        # Check gene ID type
        first_gene = adata.var_names[0]
        if first_gene.startswith('ENSG'):
            print(f"  Gene IDs: ENSG → converting to Symbol")
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
        else:
            print(f"  Gene IDs: Symbol (first: {first_gene})")
        
        # Normalize
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Macrophage pre-screen
        mac_found = [g for g in MAC_MARKERS if g in adata.var_names]
        sc.tl.score_genes(adata, gene_list=mac_found, score_name='Mac_score')
        mac_mask = adata.obs['Mac_score'] > 0
        print(f"  Mac_score > 0: {mac_mask.sum()} / {adata.n_obs}")
        
        # Classifier
        hb_found = [g for g in hb_genes if g in adata.var_names]
        mat_found = [g for g in mat_genes if g in adata.var_names]
        print(f"  Classifier genes: HB={len(hb_found)}/{len(hb_genes)}, MAT={len(mat_found)}/{len(mat_genes)}")
        
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
        
        adata.obs['dataset'] = 'gse298602'
        adata.obs['sample'] = sample_name
        adata.obs['disease'] = disease
        adata.obs['disease_group'] = disease_group
        
        n_hbc = hbc_mask.sum()
        n_mat = mat_mask.sum()
        n_amb = (mac_mask & ~hbc_mask & ~mat_mask).sum()
        
        print(f"  Classification: HBC={n_hbc}, MAT={n_mat}, AMB={n_amb}, Other={(~mac_mask).sum()}")
        
        # Validation
        for g in ['FOLR2', 'CD163', 'DAB2', 'HLA-DRA']:
            if g in adata.var_names:
                hbc_expr = adata[hbc_mask, g].X.toarray().mean() if n_hbc > 0 else 0
                mat_expr = adata[mat_mask, g].X.toarray().mean() if n_mat > 0 else 0
                print(f"    {g}: HBC={hbc_expr:.3f}, MAT={mat_expr:.3f}")
        
        # Save
        output_path = os.path.join(OUTPUT_DIR, f'gse298602_{sample_name}_reclassified.h5ad')
        adata.write_h5ad(output_path)
        
        # Hofbauer only
        if n_hbc > 0:
            hb_cells = adata[hbc_mask].copy()
            hb_path = os.path.join(OUTPUT_DIR, f'gse298602_{sample_name}_hofbauer.h5ad')
            hb_cells.write_h5ad(hb_path)
        
        results.append({
            'sample': sample_name,
            'disease': disease,
            'total': adata.n_obs,
            'hofbauer': n_hbc,
            'maternal_mac': n_mat,
            'ambiguous': n_amb,
            'other': (~mac_mask).sum()
        })
        
    except Exception as e:
        print(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()
        continue

# Summary
print("\n" + "=" * 60)
print("GSE298602 Classification Summary")
print("=" * 60)
if results:
    df = pd.DataFrame(results)
    print(df.to_string(index=False))
    pe_hb = df[df['disease'].isin(['PreE_SF', 'gHTN'])]['hofbauer'].sum()
    ctrl_hb = df[df['disease'] == 'Control']['hofbauer'].sum()
    print(f"\nPE/gHTN Hofbauer: {pe_hb}")
    print(f"Control Hofbauer: {ctrl_hb}")
    print(f"Total: {df['hofbauer'].sum()}")
    
    # Merge all Hofbauer
    all_hb = []
    for _, row in df.iterrows():
        if row['hofbauer'] > 0:
            hb_path = os.path.join(OUTPUT_DIR, f'gse298602_{row["sample"]}_hofbauer.h5ad')
            if os.path.exists(hb_path):
                all_hb.append(sc.read_h5ad(hb_path))
    
    if len(all_hb) > 1:
        common_genes = set(all_hb[0].var_names)
        for ad in all_hb[1:]:
            common_genes = common_genes.intersection(set(ad.var_names))
        all_hb_sub = [ad[:, list(common_genes)] for ad in all_hb]
        combined = sc.concat(all_hb_sub, join='inner')
        combined_path = os.path.join(OUTPUT_DIR, 'gse298602_all_hofbauer.h5ad')
        combined.write_h5ad(combined_path)
        print(f"Merged: {combined.n_obs} Hofbauer → {combined_path}")

print("\nDone!")
