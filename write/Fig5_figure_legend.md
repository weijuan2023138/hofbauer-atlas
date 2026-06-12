# Figure 5. Transcriptional remodeling of Hofbauer cells in pregnancy diseases.

**(A)** UMAP visualization of Hofbauer cells colored by six disease conditions (Normal_Early, Miscarriage, Infection, Normal_Late, PE, Preterm), arranged in a 2×3 grid. Each panel shows one condition with Hofbauer subtypes distinguished by color. 

**(B)** Stacked bar plot of six Hofbauer subtype proportions across six disease groups. Miscarriage shows expansion of Vascular remodeling (36.4%) and depletion of PRKN+ Autophagy subtypes relative to Normal_Early. PE is characterized by elevated Pro-inflammatory (41.0%) and Vascular remodeling (23.8%) subtypes. Preterm exhibits the highest Pro-inflammatory proportion (58.0%). Infection is dominated by the PRKN+ Autophagy subtype (92.1%), likely reflecting pathogen-induced autophagic activation.

**(C)** GSEA Hallmark pathway enrichment dotplot comparing four diseases (Miscarriage, Infection, PE, Preterm) against their trimester-matched normal controls (Miscarriage/Infection vs Normal_Early; PE/Preterm vs Normal_Late_noTL). Dot size represents −log10(FDR); color represents normalized enrichment score (NES). TNFα/NF-κB signaling and inflammatory response are commonly upregulated in PE and Preterm, while Miscarriage shows distinct upregulation of apoptosis and p53 pathways. Infection exhibits the broadest functional activation including interferon-γ response and complement cascade.

**(D)** Gene expression Z-score dotplot showing 40 key genes across six disease groups. Genes span autophagy (SQSTM1, BNIP3, PRKN), immune signaling (NFKB1, RELB, TNF, IL1B, CXCL8), complement (C1QA-C), ECM/communication (SPP1, FN1, COL1A2, TGFB1), antigen presentation (HLA-DRA, HLA-DQA1, HLA-DQB1), and transcription factors (STAT1, STAT3, CEBPA, MAFB). Color represents Z-score (RdYlBu); dot size represents |Z-score|.

**(E)** Violin plots showing expression of four key transcription factors (STAT1, NFKB1, CEBPA, JUN) across six disease groups. Significance stars indicate FDR < 0.05 for disease vs trimester-matched normal comparisons. NFKB1 and STAT1 are globally activated in Infection and Preterm. CEBPA is downregulated in Infection. JUN shows Miscarriage-specific upregulation.

**Supplementary Figure 5A.** UpSet plot showing intersections of significant DEGs (|log2FC| > 0.5, FDR < 0.05) across Miscarriage (986), Infection (4,798), PE (2,142), and Preterm (2,985). 101 genes are differentially expressed in all four diseases, including immediate-early response (FOSB, JUNB, ZFP36) and inflammatory (CCL3, CXCL8, MIF) genes. PE and Preterm share 1,484 DEGs, reflecting shared late-gestation disease signatures.

**Supplementary Figure 5B.** CellChat LR pair communication probability comparison between control and infected conditions in the hoo_2024 dataset, showing the top 20 LR pairs by cumulative communication probability.

**Supplementary Figure 5C.** CellChat pathway-level communication strength change (log2[Infected/Control]) for the top 15 pathways, showing signaling pathway reprogramming in placental infection.

**Supplementary Figure 5D.** Violin plots of additional transcription factors (STAT3, RELB, MAFB, ID2, KLF4, FOS, IRF1, IRF8) across six disease groups.
