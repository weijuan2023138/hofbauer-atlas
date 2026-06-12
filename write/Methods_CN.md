# 方法

## 数据收集与预处理

本研究整合了8个独立scRNA-seq数据集，覆盖孕早期至孕晚期的人胎盘组织：（1）Arutyunyan等人（2023）的早孕期正常胎盘多组学数据（n=5,002个Hofbauer细胞）；（2）GSE290578（n=6,704），包含晚孕期正常和子痫前期（PE）胎盘；（3）gse214607（n=874），早孕期流产（Miscarriage）胎盘；（4）hoo_2024（n=1,975），早孕期感染（Toxoplasma/Malaria/Listeria）及对照胎盘；（5-7）gse173193（n=145）、gse183338（n=239）和gse298119（n=361），均为PE胎盘；（8）my_preterm_cohort（n=6,331），包含早产（PTL/PTNL）和足月产（TL）胎盘。此外，UCSF_Li_2026数据集（n=625个Hofbauer细胞）的中孕期正常胎盘scRNA-seq与snATAC-seq数据用于空间转录组和染色质可及性分析。所有数据集经统一质控：过滤线粒体基因比例>20%、UMI<500或基因数<200的低质量细胞。

Hofbauer细胞的鉴定采用两阶段策略。首先，使用巨噬细胞标志物集（CD68、CD14、AIF1、CSF1R、CD163、ITGAM、FCGR3A、LYZ、C1QA、C1QB）计算每个细胞的巨噬细胞评分（Mac_score），保留Mac_score>0的细胞进入第二阶段。随后，应用基于Vento-Tormo等人（2018）胎盘细胞图谱训练的Hofbauer分类器（F1=0.9935），该分类器利用Hofbauer细胞上调基因（HB-up）与母体免疫细胞上调基因（MAT-up）的表达差异（DIFF score），保留DIFF>0.32的细胞作为高置信度Hofbauer细胞。最终共保留17,896个Hofbauer细胞用于下游分析。

## 数据整合与批次校正

各数据集的Hofbauer细胞在19,940个共同基因上合并为统一表达矩阵。使用Seurat v5进行标准化（LogNormalize，scale factor=10,000）、高变基因鉴定（n=2,000）和缩放。PCA降维后（npcs=50），使用Harmony进行批次校正，以'dataset'为分组变量，对前30个PCA维度进行校正。校正后的Harmony嵌入用于UMAP降维（dims=1:30）和Louvain聚类。

## 聚类与亚型鉴定

使用Louvain算法（resolution=0.15）对Harmony校正后的数据聚类，初始得到15个cluster。通过差异基因分析（FindAllMarkers，Wilcoxon秩和检验，min.pct=0.3，|log2FC|>0.5，调整后p<0.05）鉴定每个cluster的标志基因。根据标志基因表达模式、经典Hofbauer标志物（FOLR2、CD163、CD68、CSF1R）表达水平、以及滋养细胞/红细胞/成纤维细胞标志物（KRT7、PAGE4、HBB、COL1A1）的排除性分析，剔除7个污染物cluster（C5、C7-C14）。剩余6个cluster（C0、C1、C2、C3、C4、C6）被鉴定为真正的Hofbauer细胞亚型，分别命名为：Pro-inflammatory（C0，以IL1B、CCL3、CXCL8高表达为特征）、MHCII+ Antigen-presenting（C3，HLA-DRA、HLA-DPA1、CD74高表达）、Homeostatic（C1，FOLR2、CD163、LYVE1、DAB2高表达）、PRKN+ Autophagy（C2，PRKN、SQSTM1、MAP1LC3B高表达）、Vascular remodeling（C4，SPP1、FN1、VEGFA高表达）和MKI67+ Proliferating（C6，MKI67、TOP2A、BIRC5高表达）。

## 发育分析

根据孕周将正常妊娠样本分为三个孕期组：孕早期（Early，GW4.5–10）、孕中期（Mid，GW11–24）和孕晚期（Late，GW32–38）。使用Wilcoxon秩和检验进行孕期间差异表达分析（|log2FC|>0.5，FDR<0.01）。转录因子差异分析专注于TF列表（共鉴定47个显著差异TF）。

发育功能模块评分使用AddModuleScore计算：祖细胞模块（CEBPA、ID2、SOX4、TCF4等）、增殖模块（MKI67、TOP2A、PCNA等）、代谢模块（ENO1、LDHA、PKM等）、重塑模块（SPP1、FN1、TIMP1、MMP14等）和免疫模块（C1QA、C1QB、C1QC、FCGR3A、HLA-DRA等）。模块评分在三个孕期间使用Kruskal-Wallis检验进行比较。

基因集富集分析（GSEA）使用clusterProfiler和fgsea包，以Hallmark基因集（MSigDB v2023.2）为参考，基因按log2FC排序后进行预排序GSEA（minSize=10，maxSize=500），FDR<0.1视为显著富集。

## ATAC-seq分析

UCSF_Li_2026数据集中经FACS纯化的Hofbauer细胞进行snATAC-seq。每个样本独立使用MACS3进行peak calling（effective genome size=2.7×10⁹），各样本peak取并集后经reduce合并为consensus peak set。使用FeatureMatrix对consensus peaks重新计数，构建peak×cell矩阵。经TF-IDF归一化后，使用LSI降维（dims=2:30）和UMAP可视化。细胞按孕周分为中期（Mid，ZY012/014/020，n=1,230）和足月（Term，ZY011/019，n=925）两组。

差异可及性分析使用FindMarkers（test.use="LR"，latent.vars="nCount_ATAC"，min.pct=0.05），|log2FC|>0.25且FDR<0.05的peak定义为差异可及性peak。

De novo motif富集分析：对足月特异性开放peak，使用MEME Suite中的AME工具进行已知motif扫描（JASPAR 2022数据库），E-value<1×10⁻⁵视为显著富集。FIMO用于单个peak中motif出现频率的定量（p<1×10⁻⁴），统计每个motif在中期和足月peak中的覆盖率百分比，使用Fisher精确检验比较组间差异。

对于图6中通讯基因调控元件的motif分析，提取通讯基因±50kb范围内的差异peak，使用FIMO扫描STAT3（MA0144.3）、STAT1（MA0137.3）、RELB（MA1117.1）、CEBPA（MA0102.4）和NFKB1（MA0105.4）的结合基序。

ATAC-RNA整合分析：将ATAC差异peak对应的最近基因与RNA-seq差异表达基因取交集，比较染色质可及性变化（ATAC log2FC）和转录变化（RNA log2FC）的方向一致性，使用Cohen's κ评估一致性程度。

## 空间转录组分析

使用Li等人（2026）通过Stereo-seq技术生成的16例中孕期（GW18⁺⁶–24⁺¹）健康胎盘基底板空间转录组数据。以bin50（~50 μm直径）为分析单位，经无监督聚类和已知标志物注释鉴定10种胎盘细胞类型。

Hofbauer细胞的邻域分析：以每个Hofbauer细胞bin为中心，50 μm半径圆形区域为邻域范围。计算邻域内各细胞类型的观察比例，与10,000次空间随机置换的期望比例比较，计算富集倍数（observed/expected）和经验p值（Bonferroni校正）。

配体-受体空间共定位验证：对选定的L-R对，以高表达配体的bin为种子点、高表达受体的bin为靶点，计算50 μm半径内的观察共定位事件数，与100次空间置换比较计算富集倍数。

## 细胞间通讯分析

使用CellChat v2对中孕期正常胎盘scRNA-seq数据进行配体-受体互作分析。使用CellChatDB.human数据库，以默认参数运行：truncatedMean用于计算平均表达（trim=0.1），projectData用于投射到蛋白-蛋白互作网络，computeCommunProb使用triMean方法。对Hofbauer细胞与7种胎盘细胞类型（FB、vEC、fEC、SCT、VCT、dNK、CD14_M）之间的入方向和出方向通讯进行系统分析。通讯概率>0且p<0.05的L-R互作视为显著。感染条件下（hoo_2024数据集）的CellChat分析使用相同参数，分别对Control和Infected组独立运行后比较。

## 疾病分析

四种妊娠疾病分别与孕周匹配的正常对照比较：流产（n=874）和感染（n=1,975）与正常早孕期（Normal_Early）比较，PE（n=5,743）和早产（n=3,716）与正常晚孕期（Normal_Late_noTL，排除足月产）比较。差异表达分析使用FindMarkers（Wilcoxon检验，|log2FC|>0.5，FDR<0.05），GSEA使用fgsea（Hallmark基因集，minSize=10，maxSize=500，FDR<0.1）。

跨疾病保守差异基因分析使用UpSet图：四个疾病分别与各自对照比较，提取显著DEG（|log2FC|>0.5，FDR<0.05），取交集鉴定共同差异基因。

## 双轨TF调控模型分析

TF-通讯基因相关性：对CEBPA、STAT3、STAT1、NFKB1和RELB五个TF与30个核心通讯基因在13,582个Hofbauer细胞中计算Spearman秩相关系数。

TF分组验证：按TF表达水平将Hofbauer细胞分为高表达组（top 30%）和低表达组（bottom 30%），使用Wilcoxon检验比较两组间通讯基因表达差异，计算表达倍数变化。

TF活性分析：为每个关键TF构建靶基因调控子（regulon），使用已知靶基因集的平均表达作为TF活性读数。通过比较TF mRNA表达与调控子活性的线性回归斜率（R²）评估mRNA-活性耦合强度。

ECM-免疫双维空间分析：定义ECM模块评分（FN1、SPP1、COL1A1、MMP9的AddModuleScore）和免疫模块评分（IL1B、TNF、CXCL8、CD44、CD47），在各疾病组中计算两个维度的95%置信椭圆，比较组间偏移。

## 统计方法

所有统计分析在R 4.3环境下进行。两组比较使用双尾Wilcoxon秩和检验，多组比较使用Kruskal-Wallis检验，分类变量比较使用Fisher精确检验或卡方检验。批次效应校正使用Harmony（group.by.vars='dataset'）。p值经Benjamini-Hochberg方法进行多重检验校正。除非特别说明，显著性标准为*p<0.05、**p<0.01、***p<0.001。箱线图中，盒体示四分位距（IQR），横线示中位数，须线延伸至1.5×IQR范围内的最远端数据点。

## 数据可用性

本研究中使用的公开数据集可通过GEO数据库获取：GSE290578、gse214607、gse173193、gse183338、gse298119。Arutyunyan等人（2023）的多组学数据可通过其论文获取。hoo_2024数据待发表。STOMICS空间转录组数据来自Li等人（2026，GEO: XXXXXX）。所有分析代码可在GitHub获取：[待补充]。
