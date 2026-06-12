# 方法

## 数据收集与Hofbauer细胞鉴定

本研究整合了9个独立scRNA-seq数据集的人胎盘组织Hofbauer细胞（胎盘驻留巨噬细胞）。数据集包括：（1）Arutyunyan等人（2023）的早孕期正常胎盘多组学数据（n=5,002个Hofbauer细胞）；（2）GSE290578（n=6,704），包含晚孕期正常胎盘和子痫前期（preeclampsia，PE）胎盘；（3）GSE214607（n=874），早孕期流产（miscarriage）胎盘；（4）hoo_2024（n=1,975），包含早孕期感染胎盘（刚地弓形虫Toxoplasma gondii、恶性疟原虫Plasmodium falciparum、单核细胞增生李斯特菌Listeria monocytogenes）及正常对照；（5-7）GSE173193（n=145）、GSE183338（n=239）和GSE298119（n=361），均为PE胎盘；（8）my_preterm_cohort（n=6,331），包括自发性早产临产（preterm labor，PTL）、早产未临产（preterm no labor，PTNL）和足月临产（term labor，TL）三组；（9）UCSF_Li_2026（Li等人，2026，n=625），涵盖正常早、中、晚孕期胎盘，该数据集同时提供了snATAC-seq和Stereo-seq空间转录组数据。所有数据集经统一质控：剔除线粒体基因比例超过20%、UMI计数低于500或基因数低于200的低质量细胞。对于使用Ensembl基因ID的数据集，以Arutyunyan数据集的gene_symbols列为参考，将ENSG ID转换为基因符号后进行统一分析。

Hofbauer细胞的鉴定采用两阶段策略。首先使用巨噬细胞核心标志物集（CD68、CD14、AIF1、CSF1R、CD163、ITGAM、FCGR3A、LYZ、C1QA、C1QB）对每个细胞计算巨噬细胞评分（Mac_score），仅保留Mac_score大于0的细胞进入第二阶段。随后应用基于Vento-Tormo等人（2018）胎盘细胞图谱训练的Hofbauer分类器（F1=0.9935），该分类器利用Hofbauer细胞上调基因集（HB-up）与母体巨噬细胞上调基因集（MAT-up）的表达差异（DIFF score）区分胎儿Hofbauer细胞与母体免疫细胞，保留DIFF大于0.32的细胞作为高置信度Hofbauer细胞。my_preterm_cohort数据集最初使用cell_type_fine注释，经重新应用上述分类器后（Mac_score>0且DIFF>0.32）从10,036个候选细胞中鉴定出6,331个高置信度Hofbauer细胞。最终在所有数据集中共鉴定出22,256个Hofbauer细胞。

## 数据整合、批次校正与亚型鉴定

各数据集提取的Hofbauer细胞在19,940个共同基因上合并为统一表达矩阵，转换为Seurat v5对象后，使用LogNormalize方法（scale factor=10,000）进行标准化，鉴定前2,000个高变基因，并对表达矩阵进行缩放。经主成分分析降维（npcs=50）后，使用Harmony进行批次效应校正，以数据集来源（dataset）为分组变量，对前30个PCA维度进行校正。校正后的Harmony嵌入矩阵用于UMAP降维（dims=1:30）和Louvain聚类。在resolution=0.15的条件下，初步鉴定出15个cluster。使用Wilcoxon秩和检验（FindAllMarkers，min.pct=0.3，|log2FC|>0.5，Benjamini-Hochberg校正后p<0.05）鉴定每个cluster的标志基因。结合经典Hofbauer标志物（FOLR2、CD163、CD68、CSF1R）的表达水平、Hofbauer分类器评分（DIFF score）、以及滋养细胞标志物（KRT7、PAGE4、PEG10）、红细胞标志物（HBB、HBA1）、成纤维细胞标志物（COL1A1、DCN）的排除性分析，将9个cluster（C5、C7-C14，共4,360个细胞）判定为污染物并剔除。剩余6个cluster（C0、C1、C2、C3、C4、C6；共17,896个细胞）被保留为真正的Hofbauer细胞亚型。每个亚型的生物学身份根据其标志基因的表达模式和已知功能命名：C0为促炎性亚型（Pro-inflammatory，标志基因为IL1B、CCL3、CCL4、CXCL8、NLRP3），C3为MHCII+抗原呈递亚型（MHCII+ Antigen-presenting，标志基因为HLA-DRA、HLA-DPA1、HLA-DQB1、CD74），C1为稳态亚型（Homeostatic，标志基因为FOLR2、CD163、LYVE1、DAB2、MRC1），C2为PRKN+自噬亚型（PRKN+ Autophagy，标志基因为PRKN、SQSTM1、MAP1LC3B、BNIP3），C4为血管重塑亚型（Vascular remodeling，标志基因为SPP1、FN1、VEGFA、TIMP1），C6为MKI67+增殖亚型（MKI67+ Proliferating，标志基因为MKI67、TOP2A、BIRC5、PCNA）。

## 发育轨迹与功能模块分析

正常妊娠样本按孕周分为三个孕期组：孕早期（Early，GW4.5–10）、孕中期（Mid，GW11–24）和孕晚期（Late，GW32–38）。使用Wilcoxon秩和检验进行孕期间的差异表达分析（|log2FC|>0.5，FDR<0.01）。针对转录因子（TF）的专项分析使用AnimalTFDB和CIS-BP数据库的交集TF列表，以孕晚期与孕早期的Hofbauer细胞为比较对象，共鉴定出47个显著差异表达的TF。

为量化单细胞水平的功能程序，使用Seurat的AddModuleScore函数定义了五个功能模块：祖细胞模块（包含CEBPA、ID2、SOX4、TCF4等发育调控因子）、增殖模块（MKI67、TOP2A、PCNA等细胞周期基因）、代谢模块（ENO1、LDHA、PKM等糖酵解和氧化磷酸化基因）、重塑模块（SPP1、FN1、TIMP1、MMP14等ECM相关基因）和免疫模块（C1QA、C1QB、C1QC、FCGR3A、HLA-DRA等免疫效应基因）。模块评分在三个孕期间使用Kruskal-Wallis检验进行差异比较。

基因集富集分析（GSEA）使用fgsea包（1.28.0），以MSigDB Hallmark基因集（v2023.2）为参考数据库。将差异表达基因按log2FC降序排列后作为预排序基因列表，使用fgsea进行基因集排列检验（minSize=10，maxSize=500，nperm=10,000），Benjamini-Hochberg校正后FDR<0.1的基因集视为显著富集。

## ATAC-seq染色质可及性分析

利用UCSF_Li_2026数据集中经FACS纯化的Hofbauer细胞进行snATAC-seq分析。使用MACS3对每个样本的Tn5转座酶片段独立进行peak calling（effective genome size=2.7×10⁹），各样本的peak取并集后经GenomicRanges的reduce函数合并为consensus peak set。使用Signac的FeatureMatrix对consensus peaks在各样本中重新计数，构建peak×cell二值化矩阵。经TF-IDF归一化后，使用 latent semantic indexing（LSI，dims=2:30）进行降维和UMAP可视化。细胞按孕周分为中期组（Mid，样本ZY012/014/020，n=1,230）和足月组（Term，样本ZY011/019，n=925）。差异可及性分析使用Signac的FindMarkers函数（test.use="LR"，latent.vars="nCount_ATAC"，min.pct=0.05），以|log2FC|>0.25且FDR<0.05为显著差异可及性peak的阈值。

De novo motif富集分析针对足月特异性开放peak进行：使用MEME Suite中的AME工具（Analysis of Motif Enrichment）对JASPAR 2022脊椎动物核心motif数据库进行扫描，E-value<1×10⁻⁵视为显著富集。单个peak中motif出现频率的定量使用FIMO（p<1×10⁻⁴），统计每个motif在中期和足月peak中的覆盖率百分比，使用Fisher精确检验比较组间差异。对于图6中通讯基因调控元件的motif扫描，提取30个核心通讯基因转录起始位点±50kb范围内的差异peak序列，使用FIMO分别扫描STAT3（MA0144.3）、STAT1（MA0137.3）、RELB（MA1117.1）、CEBPA（MA0102.4）和NFKB1（MA0105.4）的结合基序。

ATAC-RNA整合分析：将ATAC差异peak按基因组坐标分配给最近的基因（GREAT方法），提取同时具有ATAC差异可及性和RNA差异表达数据的基因（n=2,141），比较染色质可及性变化方向（ATAC log2FC）与转录变化方向（RNA log2FC）的一致性，使用Cohen's κ系数评估一致性程度，Fisher精确检验评估显著性。

基因座水平的ATAC片段覆盖度可视化：对SPP1、CD47、PTPRM和TGFB1基因座，使用Signac的CoveragePlot函数在±50kb窗口内分别展示中期和足月组的Tn5插入片段堆叠密度，使用tabix索引的fragment文件进行快速读取。

## 空间转录组分析

使用Li等人（2026）通过Stereo-seq技术（华大基因，500 nm空间分辨率）生成的16例中孕期（GW18⁺⁶–24⁺¹）健康胎盘基底板空间转录组数据。以bin50（约50 μm直径，含10–30个细胞）为分析单位，经无监督聚类（Leiden算法）和已知标志物注释鉴定出10种胎盘细胞类型：Hofbauer细胞（FOLR2⁺CD163⁺CD68⁺）、成纤维细胞（FB，COL1A1⁺DCN⁺LUM⁺）、绒毛血管内皮细胞（vEC，PECAM1⁺CDH5⁺）、胎儿内皮细胞（fEC）、合体滋养细胞（SCT，CGA⁺CGB⁺）、绒毛细胞滋养细胞（VCT，PAGE4⁺PEG10⁺）、蜕膜NK细胞（dNK，KLRF1⁺PRF1⁺）、CD14⁺单核/巨噬细胞（CD14_M）等。

Hofbauer细胞的空间邻域分析以每个Hofbauer细胞bin为中心，50 μm半径圆形区域为邻域范围。计算邻域内各细胞类型的观察比例，与10,000次全组织范围内随机抽取等数量坐标的空间置换期望比较，计算富集倍数（observed/expected ratio）和双边经验p值（经Bonferroni校正，p<0.01视为显著富集）。

配体-受体空间共定位验证针对四对关键L-R对（SPP1→ITGAV、FN1→ITGB1、COL1A2→ITGA1、PTPRM↔PTPRM）进行。以高表达配体的bin为种子点、高表达受体的bin为靶点（表达量>第75百分位数），计算50 μm半径内的观察共定位事件数，与100次空间位置随机置换后的期望共定位比较，计算富集倍数和置换p值。

## 细胞间通讯分析

全细胞类型配体-受体互作分析使用CellChat v2（1.6.1），以中孕期正常胎盘scRNA-seq数据（UCSF_Li_2026数据集）为输入。使用CellChatDB.human配体-受体数据库，以默认参数运行：truncatedMean方法计算基因平均表达（trim=0.1），projectData投射到蛋白-蛋白互作网络（PPI.human），computeCommunProb使用triMean方法计算通讯概率，filterCommunication过滤通讯概率低于阈值的互作（min.cells=10）。对Hofbauer细胞与7种胎盘细胞类型（FB、vEC、fEC、SCT、VCT、dNK、CD14_M）之间的入方向和出方向通讯进行系统分析。通讯概率大于零且permutation test p<0.05的L-R互作视为具有统计学意义。分析结果以气泡图（netVisual_bubble）、弦图（netVisual_chord_gene，14×14英寸，300dpi）和通路级统计汇总展示。

感染条件下的CellChat对比分析使用hoo_2024数据集（n=158,978细胞，15种细胞类型），分别对Control组和Infected组独立运行CellChat（参数同上），比较两组间入/出方向L-R对数量、通讯概率排名前15的通路变化以及关键分子（SPP1、FN1、COLLAGEN、TGF-β）的通讯强度。使用Wilcoxon检验比较两组间通路通讯概率的差异。

Hofbauer细胞亚型间通讯分析使用CellChat对6个亚型间的配体-受体互作进行独立分析，输入数据为Hofbauer_Atlas_Final.rds中按亚型标签分组的表达矩阵。每个亚型随机下采样至最多200个细胞（set.seed=42），独立运行CellChat标准流程。

## 疾病差异分析

四种代表性妊娠疾病分别与孕周匹配的正常对照进行比较：流产（n=874）和感染（n=1,975）与早孕期正常对照（Normal_Early）比较；PE（n=5,743，来自GSE290578、GSE173193、GSE183338、GSE298119四个数据集）和早产（n=3,716，来自my_preterm_cohort的PTL和PTNL组）与晚孕期正常对照（Normal_Late，排除足月产TL组以避免临产效应的混杂）比较。差异表达分析使用Wilcoxon秩和检验（FindMarkers，|log2FC|>0.5，FDR<0.05）。GSEA使用fgsea（Hallmark基因集，参数同发育分析），以各疾病与其对照的差异表达基因log2FC排序值为输入。

跨疾病保守差异基因分析：将四种疾病各自鉴定的显著DEG（|log2FC|>0.5，FDR<0.05）取交集，使用UpSetR包可视化四种疾病间的DEG重叠模式，鉴定在两种及以上疾病中共同差异表达的基因。

## 双轨TF调控模型分析

TF-通讯基因相关性分析：选取CEBPA、STAT3、STAT1、NFKB1和RELB五个关键转录因子与30个核心通讯基因（包括ECM相关基因SPP1、FN1、COL1A1、COL1A2和免疫粘附相关基因PTPRM、CD44、CD47等），在13,582个Hofbauer细胞中计算Spearman秩相关系数。

TF分组验证分析：按各TF的mRNA表达水平将Hofbauer细胞分为高表达组（top 30%分位）和低表达组（bottom 30%分位），使用Wilcoxon秩和检验比较两组间通讯基因的表达差异，以表达倍数变化（fold change）量化调控效应的大小。同时比较两组间基因表达零值比例的差异以评估TF对表达概率的调控。

TF活性-表达耦合分析：为每个关键TF构建靶基因调控子（regulon），利用AnimalTFDB和TRRUST数据库获取已知靶基因集，以靶基因集的平均表达作为TF活性的功能性读数。通过线性回归比较各TF的mRNA表达水平与调控子活性之间的斜率（R²），评估mRNA表达向靶基因程序的传导效率。mRNA-活性耦合在四种疾病和正常对照中分别计算。

ECM-免疫双维空间分析：以ECM模块评分（FN1、SPP1、COL1A1、MMP9的AddModuleScore均值）和免疫模块评分（IL1B、TNF、CXCL8、CD44、CD47的AddModuleScore均值）为两个维度，在各疾病组中计算二维核密度估计和95%置信椭圆（使用ggplot2的stat_ellipse），比较不同疾病组在ECM-免疫双维空间中的整体偏移方向和幅度。

## 统计方法与可视化

所有统计分析在R 4.3.0环境下进行。两组间比较使用双尾Wilcoxon秩和检验（wilcox.test），多组比较使用Kruskal-Wallis检验，分类变量关联性使用Fisher精确检验。相关性分析使用Spearman秩相关系数。多重检验校正统一使用Benjamini-Hochberg方法。除非特别说明，显著性水平标注为：*p<0.05、**p<0.01、***p<0.001，"ns"表示不显著（p≥0.05）。箱线图中，盒体表示四分位距（IQR），盒内横线表示中位数，须线延伸至1.5倍IQR范围内的最远端数据点，超出须线范围的数据点以散点单独显示。

所有图形使用ggplot2（3.5.0）生成，配色方案统一：亚型配色为促炎性#C62828、MHCII+#E65100、稳态#1565C0、PRKN+#6A1B9A、血管重塑#2E7D32、增殖#455A64；孕期配色为早期#4575B4、晚期#D73027；火山图中疾病上调为#D73027、正常上调为#4575B4、不显著为灰色；GSEA热图使用蓝-白-红渐变。除特殊说明外，所有主图以300dpi分辨率输出为PNG和PDF双格式，白色背景，单行粗体标题居中，图例置于右侧。

## 数据可用性

本研究中使用的公开数据集可通过GEO数据库获取：GSE290578、GSE214607、GSE173193、GSE183338、GSE298119。Arutyunyan等人（2023）的多组学数据可通过其论文获取。hoo_2024数据待发表。STOMICS空间转录组数据来自Li等人（2026，GEO accession待补充）。my_preterm_cohort为自有数据。所有分析代码可在GitHub获取（https://github.com/待补充）。论文中未包含的补充数据和中间分析结果可向通讯作者合理请求获取。
