# 图1. 人类妊娠全孕期Hofbauer细胞发育的单细胞转录组图谱。

（A）研究设计和数据整合策略的模式图。整合8个独立scRNA-seq数据集，覆盖孕早期（GW4.5–10）、孕中期（GW11–24）和孕晚期（GW32–38），经质控后保留17,896个Hofbauer细胞用于下游分析。

（B）17,896个Hofbauer细胞的UMAP可视化，按6个转录亚型着色：促炎性（Pro-inflammatory，标志基因IL1B、CCL3、CXCL8）、MHCII+抗原呈递（MHCII+ Antigen-presenting，HLA-DRA、HLA-DPA1、CD74）、稳态（Homeostatic，FOLR2、CD163、LYVE1、DAB2）、PRKN+自噬（PRKN+ Autophagy，PRKN、SQSTM1、MAP1LC3B）、血管重塑（Vascular remodeling，SPP1、FN1、VEGFA）和MKI67+增殖（MKI67+ Proliferating，MKI67、TOP2A、BIRC5）。亚型经Harmony批次校正后通过Louvain聚类（resolution = 0.15）鉴定。

（C）按孕期分组的UMAP图，展示Hofbauer细胞转录状态沿孕周的系统性漂移：孕早期细胞集中于左侧区域（促炎性/增殖性），孕中期细胞向中心扩展（稳态），孕晚期细胞右移（抗原呈递/血管重塑）。

（D）六种Hofbauer亚型在各孕期中的比例折线图。促炎性亚型从孕早期的42.1%骤降至孕晚期的7.3%，稳态亚型从8.4%扩张至35.7%，MHCII+抗原呈递亚型从3.2%增至28.9%。

（E）孕晚期与孕早期Hofbauer细胞的基因集富集分析（GSEA），采用Hallmark基因集。孕早期富集通路包括E2F靶基因、MYC靶基因和氧化磷酸化；孕晚期富集通路以NF-κB介导的TNFα信号、补体级联反应、干扰素γ应答和IL6-JAK-STAT3信号为特征。NES，标准化富集得分；FDR，错误发现率。

（F）五个功能模块（祖细胞、增殖、代谢、重塑、免疫）在各孕期中的模块评分。箱线图展示中位数和四分位距。祖细胞模块单调下降（中位0.42至−0.31），免疫模块呈最大幅度单调上升（−0.48至0.56），重塑模块在孕中期达到峰值。

（G）20个代表性基因的点图，按三个功能层级组织：发育调控层（TREM2、CEBPA、SOX4、BHLHE41）、组织重塑层（SPP1、FN1、TIMP1、MMP14）和免疫效应层（FCGR3A、C1QA、C1QB、C1QC、CXCL8）。点大小代表表达细胞比例，颜色深浅代表平均表达水平。C1q补体成分从孕早期到孕晚期上调超过100倍。

补充图1A. 所有胎盘细胞类型的UMAP可视化及分类器注释，确认Hofbauer细胞（FOLR2、CD163、CD68、CSF1R）与成纤维细胞、内皮细胞、滋养细胞及其他免疫群体形成界限清晰的独立群体。

补充图1B. 数据集组成概览，展示样本在疾病状态（正常、子痫前期、流产、感染）和孕周上的分布，形成多维度表型锚定矩阵。

补充图1C. 六种Hofbauer亚型中典型标志基因（FOLR2、CD163、CD68、HLA-DRA、SPP1、FN1、MKI67、PRKN等）的点图，确认各亚型的转录身份和特异性。

补充图1D. Hofbauer亚型扩展标志基因补充点图。

补充图1X. 8个Hofbauer特征蛋白（FOLR2、CD163、MRC1、CD68、FCGR3A、HLA-DRA、TREM2、C1QA）的免疫荧光和流式细胞术蛋白水平验证，确认与转录表达模式的一致性。
