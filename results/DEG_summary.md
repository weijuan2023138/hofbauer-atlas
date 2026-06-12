# Hofbauer Atlas — 最终差异基因表（6亚型）

## 6个最终亚型

| 亚型 | 标志物 | 颜色 |
|------|--------|------|
| Pro-inflammatory | IL1B, CCL3, CXCL8 | #C62828 |
| MHCII+ Antigen-presenting | HLA-DRA, HLA-DPA1, CD74 | #E65100 |
| Homeostatic | FOLR2, CD163, LYVE1, DAB2 | #1565C0 |
| PRKN+ Autophagy | PRKN, SQSTM1, MAP1LC3B | #6A1B9A |
| Vascular remodeling | SPP1, FN1, VEGFA | #2E7D32 |
| MKI67+ Proliferating | MKI67, TOP2A, BIRC5 | #455A64 |

## 差异基因文件

- **all_markers_final.csv** — 全部显著差异基因 (padj < 0.05, |log2FC| > 0.5)
- **Hofbauer_Atlas_Final_top10_named.csv** — 每个亚型Top 10 markers

## 11→6 cluster映射

原始15个cluster经污染物剔除后合并为6个最终亚型。C4 (Trophoblast-associated), C7 (SPP1+ Remodeling), C8 (C1Q+ Complement) 等被合并/剔除，详见工作记录。
