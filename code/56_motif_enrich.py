#!/usr/bin/env python3
"""Motif enrichment: scan known TF motifs in ATAC peaks, Fisher test Term vs Mid"""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from pyjaspar import jaspardb
import matplotlib.pyplot as plt

# Load sequences
def load_fasta(path):
    seqs, current = [], []
    for line in open(path):
        if line.startswith('>'):
            if current: seqs.append(''.join(current))
            current = []
        else:
            current.append(line.strip().upper())
    if current: seqs.append(''.join(current))
    return seqs

print("Loading sequences...")
term_seqs = load_fasta("results/atac_term_peaks.fa")
mid_seqs  = load_fasta("results/atac_mid_peaks.fa")
print(f"Term: {len(term_seqs)} peaks, Mid: {len(mid_seqs)} peaks")

# Load JASPAR core vertebrate motifs
print("Loading JASPAR motifs...")
jdb = jaspardb()
motifs = jdb.fetch_motifs(collection="CORE", tax_group="vertebrates")
motifs = [m for m in motifs if 6 <= len(m) <= 20]
print(f"Total motifs: {len(motifs)}")

# Scan: simple PWM scoring
def scan_motif(motif, seqs, threshold=0.80):
    """Count sequences containing motif hit above relative score threshold"""
    pwm = {b: np.array([motif.pwm[b][i] for i in range(len(motif))]) for b in "ACGT"}
    max_score = sum(max(pwm[b][i] for b in "ACGT") for i in range(len(motif)))
    min_score = sum(min(pwm[b][i] for b in "ACGT") for i in range(len(motif)))
    count = 0
    for seq in seqs:
        best = -1e9
        for i in range(len(seq) - len(motif) + 1):
            sub = seq[i:i+len(motif)]
            if any(c not in "ACGT" for c in sub): continue
            score = sum(pwm[c][j] for j,c in enumerate(sub))
            best = max(best, score)
        rel_score = (best - min_score) / (max_score - min_score + 1e-9)
        if rel_score > threshold: count += 1
    return count

# Focus on immune/developmental TFs
tf_list = [
    # Early TFs
    "CEBPA","CEBPB","CEBPD","MAFB","MAF","MAFK",
    "FOS","FOSL1","FOSL2","JUN","JUNB","JUND",
    "ID2","SOX4","HMGA2","TCF3","TCF4","LEF1",
    # Late TFs (immune)
    "NFKB1","NFKB2","RELA","RELB","IRF1","IRF2","IRF3","IRF4",
    "IRF5","IRF7","IRF8","STAT1","STAT2","STAT3","STAT5A","STAT6",
    "SPI1","ETS1","BATF","BATF3","RUNX1","RUNX3",
    "NR4A1","NR4A2","NR4A3","KLF4","KLF6","EGR1","EGR2",
]

results = []
n_term = len(term_seqs)
n_mid  = len(mid_seqs)

for mi, m in enumerate(motifs):
    if m.name not in tf_list: continue
    if mi % 20 == 0: print(f"  {mi}/{len(motifs)}...")
    
    t_hit = scan_motif(m, term_seqs)
    m_hit = scan_motif(m, mid_seqs)
    
    # Fisher exact test
    tbl = [[t_hit, n_term - t_hit],
           [m_hit, n_mid  - m_hit]]
    try:
        _, pval = fisher_exact(tbl, alternative='two-sided')
    except:
        pval = 1.0
    
    t_pct = 100 * t_hit / n_term
    m_pct = 100 * m_hit / n_mid
    
    results.append({
        'name': m.name, 'id': m.matrix_id,
        'term_pct': t_pct, 'mid_pct': m_pct,
        'delta': t_pct - m_pct, 'pval': pval,
        'term_hit': t_hit, 'mid_hit': m_hit
    })

df = pd.DataFrame(results)
df = df.sort_values('pval')
df['neg_log10p'] = -np.log10(df['pval'].clip(1e-10))
df['sig'] = df['pval'] < 0.05

print(f"\nSignificant motifs (p<0.05): {df.sig.sum()}/{len(df)}")
print("\nTop 15:")
for _, r in df.head(15).iterrows():
    print(f"  {r['name']:15s}  Term={r['term_pct']:5.1f}%  Mid={r['mid_pct']:5.1f}%  "
          f"Δ={r['delta']:+5.1f}%  p={r['pval']:.1e}")

# Save
df.to_csv("results/atac_motif_enrichment.csv", index=False)

# Plot top 12
top = df.head(12).iloc[::-1]
colors = ['#D73027' if d>0 else '#4575B4' for d in top['delta']]

fig, ax = plt.subplots(figsize=(7, 5))
ax.barh(range(len(top)), top['delta'], color=colors, height=0.6, edgecolor='#4d4d4d', lw=0.3)
ax.set_yticks(range(len(top)))
ax.set_yticklabels(top['name'], fontsize=10, fontweight='bold')
ax.axvline(0, color='black', lw=0.5)
ax.set_xlabel('Δ Motif Enrichment (Term% - Mid%)', fontweight='bold')
ax.set_title('TF Motif Enrichment in ATAC Peaks', fontweight='bold', fontsize=13)
# Add p-value stars
for i, (_, r) in enumerate(top.iterrows()):
    stars = '***' if r['pval']<0.001 else '**' if r['pval']<0.01 else '*' if r['pval']<0.05 else ''
    x = r['delta'] + (1.5 if r['delta']>0 else -1.5)
    ax.text(x, i, stars, ha='center', va='center', fontsize=9, color='#333333', fontweight='bold')
plt.tight_layout()
plt.savefig("figures/Fig2/Fig2_motif_enrichment.png", dpi=300, bbox_inches='tight')
print("\nSaved figures/Fig2/Fig2_motif_enrichment.png")
