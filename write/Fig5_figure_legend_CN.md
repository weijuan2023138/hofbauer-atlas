# 图5. 妊娠疾病中Hofbauer细胞的转录重塑。

（A）UMAP可视化，六种疾病条件（Normal_Early、Miscarriage、Infection、Normal_Late、PE、Preterm）2×3排列，每面板按Hofbauer亚型着色。

（B）六种Hofbauer亚型在六组疾病中的比例堆叠柱状图。流产组血管重塑亚型扩张（36.4%），PRKN+自噬亚型缺失；PE组促炎性亚型（41.0%）和血管重塑亚型（23.8%）升高；早产组促炎性比例最高（58.0%）；感染组以PRKN+自噬亚型为主（92.1%）。

（C）GSEA Hallmark通路富集点图，四种疾病vs孕周匹配正常对照。点大小为−log10(FDR)，颜色为NES。TNFα/NF-κB信号和炎症应答为PE及早产共同上调通路；流产以上调凋亡和p53通路为特征；感染呈现最广泛的功能激活（干扰素γ应答、补体级联）。

（D）40个关键基因在六组疾病中的Z-score表达点图。涵盖自噬（SQSTM1、BNIP3、PRKN）、免疫信号（NFKB1、RELB、TNF、IL1B、CXCL8）、补体（C1QA-C）、ECM/通讯（SPP1、FN1、COL1A2、TGFB1）、抗原呈递（HLA-DRA等）和转录因子（STAT1、STAT3、CEBPA、MAFB）。颜色为Z-score（RdYlBu），点大小为|Z-score|。

（E）四个关键转录因子（STAT1、NFKB1、CEBPA、JUN）在六组中的小提琴图。星号为疾病vs孕周匹配正常对照FDR<0.05。NFKB1和STAT1在感染及早产中全线激活；CEBPA在感染中下调；JUN在流产中特异性上调。

补充图5A. UpSet图展示四种疾病显著DEG（|log2FC|>0.5，FDR<0.05）交集：流产986、感染4798、PE 2142、早产2985。101个基因在四种疾病中共同差异表达，包括即早基因（FOSB、JUNB、ZFP36）和炎症基因（CCL3、CXCL8、MIF）。PE与早产共享1484个DEG。

补充图5B. hoo_2024数据集感染vs对照CellChat LR对通讯概率对比，展示top 20 LR对。

补充图5C. 感染vs对照通路层面通讯强度变化（log2[Infected/Control]），展示top 15通路。

补充图5D. 其余8个转录因子（STAT3、RELB、MAFB、ID2、KLF4、FOS、IRF1、IRF8）在六组中的小提琴图。
