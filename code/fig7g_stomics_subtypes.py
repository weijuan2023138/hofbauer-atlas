#!/usr/bin/env python3
"""Fig7g: STOMICS spatial — HB subtype proximity analysis"""
import scanpy as sc, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from pathlib import Path

OUT = Path("/home/weijuan/文档/胎盘单细胞数据/ucsf_integration/figures/Fig7")
OUT.mkdir(parents=True, exist_ok=True)

# Load STOMICS HB cells
adata = sc.read_h5ad("/home/weijuan/文档/胎盘单细胞数据/raw_data/UCSF_Li_2026/STOMICS.h5ad", backed='r')
samples = ['007','012']
mask = (adata.obs['celltype'] == 'HB') & (adata.obs['sample_id'].isin(samples))
adata_hb = adata[mask].to_memory()
sc.pp.normalize_total(adata_hb, target_sum=1e4)
sc.pp.log1p(adata_hb)
print(f"HB: {adata_hb.shape}")

# Subtype signature genes (top markers from our atlas)
subtype_sigs = {
    "Vascular remodeling":     ["SPP1","FN1","COL1A1","COL1A2","MMP9"],
    "MHCII+ Antigen-presenting":["HLA-DRA","HLA-DRB1","CD74","HLA-DPA1","HLA-DPB1"],
    "Pro-inflammatory":        ["IL1B","TNF","CXCL8","CCL3","CCL4"],
    "Homeostatic":             ["TGFB1","IGF1","CSF1","CD163","MRC1"],
    "PRKN+ Autophagy":         ["PRKN","SQSTM1","MAP1LC3B","BNIP3","OPTN"],
    "MKI67+ Proliferating":    ["MKI67","TOP2A","PCNA","CDK1","CCNB1"]
}

# Score each subtype
subtype_order = ["Vascular remodeling","MHCII+ Antigen-presenting","Pro-inflammatory",
                 "Homeostatic","PRKN+ Autophagy","MKI67+ Proliferating"]
for st in subtype_order:
    genes = [g for g in subtype_sigs[st] if g in adata_hb.var_names]
    if len(genes) >= 3:
        sc.tl.score_genes(adata_hb, gene_list=genes, score_name=st)
        print(f"  {st}: {len(genes)} genes scored")

# Assign each spot to top subtype
scores = np.column_stack([adata_hb.obs[s] for s in subtype_order])
best_subtype = np.argmax(scores, axis=1)
adata_hb.obs['top_subtype'] = [subtype_order[i] for i in best_subtype]

# Determine Sender (Vascular remodeling) and Receiver (MHCII+) masks
sender_mask = adata_hb.obs['top_subtype'] == "Vascular remodeling"
receiver_mask = adata_hb.obs['top_subtype'] == "MHCII+ Antigen-presenting"

# Spatial proximity: for each Sender spot, find nearest Receiver spot
for si, sample in enumerate(samples):
    s_mask = adata_hb.obs['sample_id'] == sample
    xy = adata_hb[s_mask].obsm['spatial']
    
    sender_idx = np.where(s_mask.values & sender_mask.values)[0]
    receiver_idx = np.where(s_mask.values & receiver_mask.values)[0]
    
    if len(sender_idx) < 5 or len(receiver_idx) < 5:
        print(f"  {sample}: insufficient cells")
        continue
    
    # Build KDTree for Receiver spots (in sample subset coordinates)
    # Map global indices to sample-local
    s_positions = np.where(s_mask.values)[0]
    global_to_local = {g: l for l, g in enumerate(s_positions)}
    sender_local = np.array([global_to_local[g] for g in sender_idx])
    receiver_local = np.array([global_to_local[g] for g in receiver_idx])
    
    tree = cKDTree(xy[receiver_local])
    distances, _ = tree.query(xy[sender_local], k=1)
    
    # Null: shuffle labels and recompute
    n_perm = 100
    null_distances = []
    all_spots = np.arange(len(xy))
    for _ in range(n_perm):
        np.random.shuffle(all_spots)
        fake_receivers = all_spots[:len(receiver_local)]
        fake_senders = all_spots[len(receiver_local):len(receiver_local)+len(sender_local)]
        null_tree = cKDTree(xy[fake_receivers])
        null_d, _ = null_tree.query(xy[fake_senders], k=1)
        null_distances.append(np.mean(null_d))
    
    null_mean = np.mean(null_distances)
    null_std = np.std(null_distances)
    actual_mean = np.mean(distances)
    z_score = (actual_mean - null_mean) / null_std if null_std > 0 else 0
    
    print(f"\n{sample}: Sender→Receiver proximity")
    print(f"  Sender spots: {len(sender_local)}, Receiver spots: {len(receiver_local)}")
    print(f"  Actual mean distance: {actual_mean:.1f}")
    print(f"  Null mean distance: {null_mean:.1f} ± {null_std:.1f}")
    print(f"  Z-score: {z_score:.2f} ({'CLOSER than random' if z_score < -2 else 'random'})")
    
    # Plot: spatial scatter with Sender (red) and Receiver (blue)
    fig, ax = plt.subplots(figsize=(8, 8))
    # All HB in grey
    ax.scatter(xy[:,0], xy[:,1], c='lightgrey', s=0.5, rasterized=True)
    # Sender
    ax.scatter(xy[sender_local,0], xy[sender_local,1], c='#D73027', s=2, alpha=0.7, label='Sender (Vasc remodeling)')
    # Receiver
    ax.scatter(xy[receiver_local,0], xy[receiver_local,1], c='#4575B4', s=2, alpha=0.7, label='Receiver (MHCII+)')
    ax.set_title(f"{sample}: Sender vs Receiver proximity\nZ={z_score:.2f} | Actual={actual_mean:.0f} vs Null={null_mean:.0f}",
                 fontsize=10, fontweight='bold')
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    ax.legend(fontsize=8, loc='upper right')
    fig.savefig(OUT / f"Fig7g_proximity_{sample}.png", dpi=300, bbox_inches='tight')
    plt.close()

# Also plot all 6 subtypes spatial distribution
fig, axes = plt.subplots(2, 6, figsize=(24, 8))
for si, sample in enumerate(samples):
    s_mask = adata_hb.obs['sample_id'] == sample
    xy = adata_hb[s_mask].obsm['spatial']
    for sti, st in enumerate(subtype_order):
        ax = axes[si, sti]
        st_mask_local = adata_hb.obs[s_mask]['top_subtype'] == st
        # Grey background
        ax.scatter(xy[:,0], xy[:,1], c='lightgrey', s=0.3, rasterized=True)
        ax.scatter(xy[st_mask_local,0], xy[st_mask_local,1], s=2, alpha=0.7, rasterized=True)
        ax.set_title(f"{sample}  {st}", fontsize=8, fontweight='bold')
        ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("STOMICS Spatial: HB Subtype Distribution", fontsize=14, fontweight='bold')
fig.tight_layout()
fig.savefig(OUT / "Fig7g_subtype_spatial.png", dpi=300, bbox_inches='tight')
plt.close()
print("\nDone: Fig7g")
