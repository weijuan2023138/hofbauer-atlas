# Hofbauer Cell Atlas — Shiny Application

交互式数据探索网站，包含6个模块：

| 模块 | 功能 |
|------|------|
| 发育图谱 | UMAP按亚型/孕期/疾病/数据集着色 |
| 基因表达 | 搜索任意基因，查看UMAP + 小提琴图 |
| 疾病对比 | 四种疾病的火山图 + 亚型比例 |
| TF调控 | 9个关键TF的表达UMAP + 分组小提琴图 |
| ECM-免疫双维空间 | ECM vs Immune模块评分散点图 |
| 数据下载 | DEG表、元数据、模块评分CSV下载 |

## 运行方式

```bash
cd ucsf_integration
Rscript -e 'shiny::runApp("shiny", port=3838)'
```

浏览器打开 `http://localhost:3838`

## 数据文件

预计算数据位于 `shiny_data/`（总计约15MB），从 `results/Hofbauer_Atlas_Final.rds` 提取。
