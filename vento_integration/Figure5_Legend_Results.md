# Figure 5: Disease-specific Hofbauer cell programs identify PE subtypes

## Figure Legend

**Figure 5. Disease-specific transcriptional programs identify two mechanistically distinct PE subtypes.**

(a) UMAP visualization of Hofbauer subtypes (Pro-inflammatory, MHCII+ Antigen-presenting, Homeostatic, PRKN+ Autophagy, Vascular remodeling, MKI67+ Proliferating) across corrected disease groups. PE is split into Early PE (GSE290578, GW29-34) and Late PE (GSE298602, GSE173193, GSE298119). Normal_Early: Arutyunyan 2023 (E-MTAB-12421) + Vento-Tormo 2018 (E-MTAB-6701) + Hoo 2024 normal controls. Normal_Late: GSE290578 normal + GSE298602 control. Disease labels corrected per per-cell metadata.

(b) Subtype proportion bar chart across disease groups. Colors match panel (a).

(c) GO Biological Process GSEA dotplot comparing each disease against its trimester-matched normal. Each point represents one GO:BP term; color indicates normalized enrichment score (NES), size indicates core gene count. Only pathways significant in at least two disease comparisons are shown.

(d) Key gene expression Z-score dotplot across disease groups. Genes grouped by functional categories.

(e) Violin plots of TF expression in PE subtypes (Normal_Late, Early_PE, Late_PE). CEBPA, IRF1, JUN, and KLF4 are shown; supplementary panels show STAT1, STAT3, IRF8, FOS, RELB, ID2, MAFB. Significance assessed by pairwise Wilcoxon test. ***P<0.001, **P<0.01, *P<0.05.

(f) UMAP visualization of FLT1 expression across PE groups.

Supplementary panels: UpSet plot of DEG overlaps; GSE298602 internal analysis (PreE_SF vs Control); PE subgroup bar charts; GSEA dotplot for PE subtypes; individual TF violin plots.

## Results

**Four diseases, four distinct molecular programs.** Disease-specific Hofbauer cell transcriptional responses were largely non-overlapping (Supplementary UpSet), with each condition activating a unique gene expression program.

**Miscarriage — acute translation burst.** Miscarriage Hofbauer cells (428 cells, vs Normal_Early) showed the strongest signal in cytoplasmic translation activation (NES=3.1) with mild mitochondrial suppression. Limited cell numbers restricted GSEA power but the dominant signal was ribosomal and translational.

**Infection — mitochondrial metabolic collapse.** Infection Hofbauer cells (923 cells, corrected to exclude within-study normal controls) showed systematic downregulation of oxidative phosphorylation, ATP synthesis, and mitochondrial translation (NES -2.0 to -2.5). Unlike PE, immune pathways were not significantly altered, supporting pathogen-driven metabolic paralysis rather than immune activation.

**PE — immune activation with subtype-specific divergence.** PE Hofbauer cells (6,659 cells) showed broad immune pathway upregulation in GSEA, but this was dominated by the early-onset subtype. Subtype proportion analysis revealed the mechanistic basis: Early PE cells were predominantly Pro-inflammatory (42%) and Vascular remodeling (24%), with MHCII+ cells suppressed to 5%. Late PE cells showed the opposite pattern — MHCII+ cells expanded to 42% with Pro-inflammatory cells reduced to 11%. This compositional switch explained the GSEA findings: ECM/anti-angiogenic pathways were early PE-specific, while antigen processing pathways were shared by both subtypes but more prominent in late PE.

**Preterm — immune hyperactivation mirroring PE.** Preterm Hofbauer cells (3,694 cells, PTL+PTNL, TL excluded) showed strong immune pathway upregulation (adaptive immunity, T cell activation, leukocyte adhesion), directly opposing PE. Mitochondrial pathways were also suppressed, shared with Infection.

**PE is two diseases.** Early PE (GW29-34) and late PE (GW37-40) showed nearly opposite Hofbauer subtype compositions. Direct GSEA comparison between the two confirmed: Early PE was biased toward anti-angiogenic programs (negative regulation of angiogenesis NES=2.16), ECM organization, and CD4+ T cell differentiation. Late PE was biased toward insulin-like growth factor signaling and epithelial morphogenesis (NES=-2.08). These differences were confirmed within the GSE298602 dataset alone (PreE_SF vs Control, within-batch), ruling out batch effects.

**TF drivers of PE subtypes.** KLF4 (logFC=0.68, P=1.1e-40) was the single most significant TF in PE, followed by IRF1 (logFC=0.82, P=7.5e-32), JUN (logFC=0.36, P=2.7e-28), and STAT1 (logFC=0.35, P=2.7e-08). CEBPA was modestly significant (logFC=0.20, P=0.02). NFKB1, IRF8, RELB, and MAFB were not significantly different in PE, contradicting previous models based on contaminated disease labels.

## Results Summary

The PE disease analysis corrected two major errors from prior analyses: (1) disease labels were corrected using per-cell metadata, removing Normal/PE contamination that had generated false positive immune activation signals; (2) PE was split into early- and late-onset subtypes, revealing two mechanistically distinct diseases. Early PE Hofbauer cells are ECM/anti-angiogenic effectors secreting sFlt-1, while late PE Hofbauer cells are antigen-presenting cells with local immune dysregulation.
