# Hofbauer Cell Atlas — 全孕期发育图谱与疾病分析

整合8个scRNA-seq数据集（17,896 Hofbauer细胞）+ ATAC-seq + STOMICS空间转录组的胎盘驻留巨噬细胞综合分析管线。

## 目录结构

```
ucsf_integration/
├── code/            # 分析脚本 (R + Python)
├── results/         # 输出数据 (RDS, CSV, h5ad)
│   └── classification/  # 分类器输出
├── figures/         # 出版用图 (300dpi PNG/PDF)
│   ├── Fig1/        # 发育图谱 + 蛋白质验证
│   ├── Fig2/        # TF开关 + ATAC-seq
│   ├── Fig3/        # 空间转录组niche分析
│   ├── Fig4/        # CellChat通讯网络
│   ├── Fig5/        # 疾病扰动分析
│   ├── Fig6/        # 双轨TF调控模型
│   ├── Fig7/        # 疾病中双轨模型破坏
│   └── FigS/        # 补充图
├── write/           # 结果描述 (中文)
├── ref/             # 参考文件 (GMT基因集)
├── deprecated/      # 废弃/迭代中间文件
│   ├── code/
│   ├── figures/
│   └── results/
├── logs/            # 运行日志
└── 工作记录.md       # 分析决策记录
```

## 图 Panels 与脚本对应

| 图 | 脚本 | 说明 |
|----|------|------|
| Fig1A | `42_fig1a_model.R` | 数据概览模型图 |
| Fig1B | `17_final_comprehensive_figures.R` | 亚型UMAP |
| Fig1C | `fig5a_umap.R` | 孕期UMAP |
| Fig1D | `fig5b_subtype_proportions.R` | 亚型比例线图 |
| Fig1E | `fig5c_gsea.R` | GSEA功能转换 |
| Fig1F | `29_module_volcano.R` | 模块评分轨迹 |
| Fig1G | `fig5d_gene_dotplot.R` | 发育基因点图 |
| Fig2a | `30_fig2a_final.R` | TF点图 |
| Fig2b | `34_atac_analysis.R` | ATAC火山图 |
| Fig2c | `fig6c_motif.R` | Motif富集 |
| Fig2d | `37_atac_rna_joint.R` | ATAC-RNA散点图 |
| Fig2e | `41_fig2b_gsea_network.R` | GSEA网络 |
| Fig2f | `20_GSEA_trimester.R` | 热图 |
| Fig3 | `32_spatial_plot.py`, `33_neighborhood.py` | 空间转录组 |
| Fig4 | `45_fig4_cellchat.R`, `63_spatial_LR_validation.py` | CellChat + 空间验证 |
| Fig5 | `fig5a_umap.R` ~ `fig5e_tf_violins.R`, `40_disease_analysis.R` | 疾病分析 |
| Fig6 | `fig6a_tf_comm_corr.R` ~ `fig6e_model.R` | 双轨调控 |
| Fig7 | `fig7_disease_dual.R`, `fig7_subtype_comm.R`, `fig7f_cellchat.R`, `fig7g_stomics_subtypes.py` | 疾病破坏 |

## 关键数据文件

| 文件 | 内容 |
|------|------|
| `results/hofbauer_final_clean.rds` | 最终Seurat对象 (17,896 cells) |
| `results/Hofbauer_Atlas_Final.rds` | 带标注的Atlas对象 |
| `results/Hofbauer_ATAC_mid_term.rds` | ATAC Mid vs Term 对象 |
| `results/cellchat_ucsf_mid.rds` | 中孕期CellChat结果 |
| `results/cellchat_hoo2024_*.rds` | 感染vs对照CellChat |
| `results/classification/all_hofbauer_final_corrected.h5ad` | 分类器输出 |

## 运行环境

- R >= 4.3, Python >= 3.10
- R包: Seurat, Signac, harmony, CellChat, clusterProfiler, slingshot, monocle3
- Python包: scanpy, squidpy, pandas, numpy
- MACS3: `~/.local/bin/macs3`

## 已知局限

1. Hi-C/ABC模型数据不存在——手稿中相关句子已标注为推断性陈述
2. STOMICS数据GEO accession待补充
3. 仅中孕期CellChat数据（PE/流产CellChat未运行）
4. 感染组全部来自hoo_2024单一数据集 (n=1,893)
