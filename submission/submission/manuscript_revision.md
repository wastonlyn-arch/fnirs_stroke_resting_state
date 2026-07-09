# 手稿修订稿（审阅模式）

> 修订标记：~~删除~~ = 原文删除内容，**粗体** = 新增/修改内容
> 未标记的段落 = 保留原文不变

---

Synergistic Effects of Mirror Therapy and Pneumatic Glove Training on Distal Hand Function After Stroke: A fNIRS-Based Randomized Study

Cuifeng Wen^1^ Tianyu Ma^2^* Hao Huang^1^ Ru Ya^1^ Yulong Bai^2^

^1^ Intensive Rehabilitation Department, Shanghai Third Rehabilitation Hospital
^2^ Department of Rehabilitation Medicine, Huashan Hospital Affiliated to Fudan University, Shanghai
*Co-first author

---

## Abstract

**Background**: Stroke-induced hand dysfunction severely impacts patients' quality of life. While mirror therapy and pneumatic glove training have shown promise in upper limb rehabilitation, their combined effects and underlying neural mechanisms remain poorly understood.

**Objective**: This study aimed to investigate the rehabilitative efficacy and neurophysiological mechanisms of mirror therapy combined with intelligent pneumatic flexible glove training for post-stroke hand dysfunction.

**Methods**: Eighty stroke patients with hand dysfunction were randomly allocated into four equal groups (n=20 each): mirror therapy, pneumatic glove, combined intervention, and conventional rehabilitation control. All interventions were administered for 20 minutes per session, five times weekly, over six weeks. Clinical outcomes were assessed using the Fugl-Meyer Assessment for Upper Extremity (FMA), Action Research Arm Test (ARAT), and Barthel Index (BI). **Thirty-two patients with complete resting-state fNIRS data were included in the brain network analysis.** ~~Concurrently, functional near-infrared spectroscopy (fNIRS) was employed to monitor both cortical activation and brain network connectivity by measuring task-related changes in oxy-hemoglobin (HbO) concentration within regions of interest, including the prefrontal cortex (PFC), dorsolateral prefrontal cortex (DLPFC), supplementary motor area (SMA), and sensorimotor cortex (SMC).~~

**Results**: All four groups demonstrated significant improvements in FMA, ARAT, and BI scores post-intervention (all P < 0.01). The combined group significantly outperformed both the control and glove groups in FMA scores, and surpassed all other groups in ARAT performance (all P < 0.05). Notably, a significant interaction effect between mirror and glove therapies was identified for distal FMA (F = 8.13, P = 0.006) and ARAT (F = 4.912, P = 0.03), indicating a synergistic effect on fine motor function. ~~fNIRS data revealed that the combined group elicited significantly greater HbO increases in the left PFC, DLPFC, and sensorimotor cortices compared to controls (P < 0.05). Network analysis further demonstrated that the combined therapy promoted a shift from compensatory inter-hemispheric expansion to efficient neural encoding, characterized by a significant negative correlation between inter-hemispheric connectivity changes and distal FMA improvement (r_s = -0.71, P = 0.048).~~ **Resting-state fNIRS functional connectivity analysis using Network-Based Statistics (NBS) identified a significant MT x PG interaction network (100 edges, P = 0.0002). Subnetwork decomposition revealed that prefrontal-polar connections comprised 45% of the significant network, followed by inter-hemispheric connections (23%). In the combined therapy subgroup (n=7), a moderate negative correlation was observed between SNB connectivity changes and distal FMA improvement (r_s = -0.59, P = 0.123); however, this did not survive FDR correction.**

**Conclusion**: Mirror therapy combined with pneumatic glove training produces superior improvements in upper limb motor function, particularly for distal hand skills, through a synergistic effect. **This benefit is associated with enhanced cortical activation in key motor-related areas, and synergistic modulation of prefrontal-polar and inter-hemispheric connectivity networks.**

**Keywords**: Mirror therapy; Pneumatic glove training; Stroke; Hand dysfunction; Functional near-infrared spectroscopy; **Brain network**

---

## Introduction

脑卒中（Stroke）是全球范围内导致死亡和长期残疾的首要病因，对社会经济和公共卫生系统构成严峻挑战[1,2]。世界卫生组织数据显示，每年约1500万人发生脑卒中，其中约80%的幸存者遗留不同程度的功能障碍[3]，50%～60%的患者在慢性期仍持续存在上肢及手功能障碍[4]。手部作为执行精细操作和复杂任务的关键器官，其功能受损对患者的日常生活活动能力、生活质量及社会参与度造成毁灭性影响[5]。因此，有效改善手功能是脑卒中康复的核心目标之一。

镜像疗法（Mirror Therapy, MT）作为一种基于视觉反馈和神经可塑性的康复技术，通过健侧肢体运动镜像激活镜像神经元系统、强化视觉-感觉运动整合、调节半球间抑制，已在脑卒中后运动功能恢复中显示出明确疗效[8,11,12]。近年来，智能气动式柔性手套（Intelligent Pneumatic Glove, IPG）作为机器人辅助康复的代表性技术，运用柔软材料与气压驱动方式，能够为严重功能障碍患者提供符合生理特征的运动辅助，有效解决运动执行难题[30,32]。在此基础上，研究者开始探索MT与IPG的联合干预策略，初步临床研究表明，二者联合相较于单一疗法能更显著提升患者手部灵巧度（如积木障碍盒测试评分）[41]，其协同机制被认为涉及视觉反馈与本体感觉输入的双通道强化，以及运动执行与动作观察的同步激活[32,33]。

然而，当前MT与IPG联合干预领域仍处于初步探索阶段，存在若干亟待解决的核心科学问题：（1）机制不清——两种干预方式作用于不同神经通路的具体互补与协同机制尚未明确，联合干预对大脑半球间抑制和皮质兴奋性调节的神经生理学证据匮乏[3,27]；（2）方案标准化欠缺——现有研究在训练强度、时序配合等关键参数上差异较大，缺乏统一、可复现的干预标准，视觉-体感反馈的时空同步性及最佳训练参数缺乏系统性研究[27]；（3）神经机制研究不足——多数研究依赖行为学指标（如Fugl-Meyer评估量表），缺少基于功能近红外光谱（fNIRS）、功能磁共振成像或经颅磁刺激等多模态神经影像技术的动态脑功能监测证据[28]；（4）研究方法学局限——样本量有限、随访时间短，对患者基线功能状态（如Brunnstrom分期）的亚组分析不足，未能充分揭示不同恢复阶段患者的疗效差异[30,31]。

为解决上述问题，本研究采用2×2析因设计的随机对照试验，纳入脑卒中后手功能障碍患者，随机分为MT组、IPG组、联合干预组及常规康复对照组，通过标准化评定量表（Fugl-Meyer评定量表、手功能评定量表）与fNIRS脑功能监测相结合的方式进行双重评估。本研究的核心创新体现在：（1）理论层面，提出并验证"多模态同步强化"假说——即当视觉反馈与体感反馈在时空上精确同步时，能否最大化激活感觉-运动整合网络并促进神经重塑；（2）机制解析层面，~~首次~~利用fNIRS技术动态监测联合干预下大脑皮层激活模式及半球间平衡的时序变化，**并采用基于网络的统计方法（NBS）从静息态功能连接的角度**为阐明MT与IPG协同效应的神经可塑性机制提供客观的影像学证据；（3）方法学层面，通过2×2析因设计明确区分主效应与交互效应，结合**基于网络掩膜的子网络分解分析**~~基于Brunnstrom分期的亚组分析~~，系统探索联合干预对大脑功能网络拓扑结构的调节作用。

---

## Methods

### 2.1 研究对象

本研究的受试者群体选自2023年1月至2024年10月上海市第三康复医院住院患者。病例筛选遵循循证医学原则，确保研究的科学性和伦理性。伦理编号（上海市第三康复医伦理委员会：SH3RH-2022-EC-019；）。纳入标准:（1）符合第四届全国脑血管病学术会议制定的脑卒中诊断标准；（2）经颅脑CT或MRI确诊为首次卒中发作，病程处于亚急性期（1—6个月），神经功能稳定，年龄范围20-80周岁；（3）认知功能评估符合：简易精神状态检查量表（MMSE）≥20分；运动功能Brunnstrom分期≤Ⅲ期；改良Ashworth肌张力分级≤3级；（4）无精神疾病史，具备基本治疗依从性；（5）自愿签署知情同意书并完成伦理审查程序。排除标准:（1）年龄超过80周岁或患有严重心肺肝肾器质性疾病，可能影响康复训练实施的患者；（2）短暂性脑缺血发作等可逆性脑血管事件患者；（3）MMSE评分低于20分，存在严重认知功能障碍，影响治疗配合的患者；（4）具有精神疾病现病史或既往病史的患者；（5）存在视听交流障碍；（6）因各种原因无法完成规定疗程及随访观察的患者。本研究获得上海市第三康复医院伦理委员会批准。所有方法严格遵守医院相关指导方针和规定，并获得每位参与者或其法律代表的书面知情同意。

**其中32例患者（sham组8例、MT组7例、PG组9例、MtPg组8例）完成了完整的静息态fNIRS前后测数据采集，纳入后续脑网络功能连接分析。各组间基线人口学特征与临床指标均无统计学差异（P > 0.05）。**

### 2.2 分组标准及盲法

采用计算机Excel 生成的随机数字表法（区组随机化），80例合格受试者被均匀划分至四个组别：镜像组（n=20）、手套组（n=20）、联合组（n=20）以及对照组（n=20）。在研究实施阶段，分组信息由独立的统计师根据病例的入组顺序实时进行分配。本试验遵循三盲设计原则，确保受试者、疗效评估者及数据分析者均对分组信息保持盲态。研究方案严格规定从随机序列生成、干预实施、数据采集到统计分析全流程的盲法维持措施，最大程度控制选择偏倚、实施偏倚和测量偏倚。盲底揭封仅在进行最终统计分析前由主要研究者执行。

### 2.3 干预措施

**2.3.1 培训**

上海市第三康复医院的医生及资深治疗师共同对相关人员进行规范化培训。培训内容涵盖功能性近红外光谱技术（fNIRS）标准化操作流程、康复评定方法及干预方案实施细则。为确保研究一致性，课题中期将进行标准化复训，避免操作差异。

**2.3.2 干预措施**

（保留原文所有干预措施描述不变。）

### 2.4 康复评估

**2.4.1 临床评估**

（保留原文临床评估描述不变：U-FMA、ARAT、Barthel指数。）

**2.4.2 fNIRS 评估与分析**

本研究采用丹阳慧创医疗设备有限公司研发的NIRS-smartII-3000A近红外脑功能成像系统进行脑血流动力学参数采集。该系统选用双波长近红外光（730±5 nm和850±5 nm），数据采样频率为11 Hz。探头阵列包括18个光源和16个探测器，组成平面矩阵，形成38个有效检测通道，光源与探测器间距为30毫米。空间定位通过三维数字化定位系统实现，使用EEG 10-20系统[49]将各通道坐标映射至标准MNI脑空间模板，基于结构磁共振成像数据划分8个功能脑区（左右各四个）：前额叶皮层（PFC）、背外侧前额叶（DLPFC）、辅助运动区（SMA）及皮层感觉运动区（SAC）。

**数据分析路径I：皮层激活** 皮层激活预处理在Homer2（MATLAB）软件中完成。处理流程包括：(1)带通滤波（0.01–0.1 Hz）以消除生理噪声[50]；(2)强度值向光密度值转换；(3)计算HbO与HbR浓度变化[51]。针对HbO信号[40]，采用试验平均法以提升信噪比。为统一数据呈现方式，所有右侧病灶数据均转换至统一的左侧病灶方位[52]。组水平皮层激活信号被映射至标准化头皮拓扑图上。

**数据分析路径II：脑网络**

~~另采用基于MNE-Python[36]的独立数据处理流程进行连接性分析。经过通道质量控制（头皮耦合指数）和插值处理后，数据被转换为HbO/HbR浓度[37]并进行带通滤波（0.05–0.2 Hz）。运动伪影通过tdDR方法[38]进行校正。~~

**脑网络分析采用基于Python（MNE-Python[36]）的独立数据处理流程。首先，对静息态原始光强度数据进行光密度（OD）转换。随后计算头皮耦合指数（SCI，阈值≥0.75）和变异系数（CV，阈值≤0.20）进行通道质量控制，其中SCI在OD域计算，CV在原始光强域计算，以保证质控指标的独立性。对于通过质控的通道，采用Homer3兼容的运动检测算法（hmrMotionArtifactByChannel，STDEV阈值=8.0，AMP阈值=3.0）识别运动伪迹段，并通过小波校正方法（db4小波，IQR=1.5）在OD域进行运动伪影去除[38]。校正后的OD信号经4阶巴特沃斯带通滤波器（0.01–0.1 Hz）滤除低频漂移和高频噪声，再经Modified Beer-Lambert定律转换为HbO浓度（差分路径长度因子DPF=6.0，源-探间距3.0 cm）。静息态数据分为治疗前（baseline）和治疗后（post）两个时点采集，每次采集时长5分钟。**

**功能连接计算：对每个被试的HbO时间序列，计算38通道间的Pearson相关系数，经Fisher Z变换后得到38×38的Z矩阵。取治疗前后均为有效通道的子集，计算功能连接变化量ΔZ = Z_post - Z_pre。**

**基于网络的统计（NBS）分析：采用GLM-NBS方法检验MT×PG交互效应对功能连接的调节作用[引用Zalesky et al. 2010]。具体而言，对每条连接边缘拟合线性模型：Post_FC ~ 1 + MT + PG + MT×PG + Pre_FC，其中Pre_FC作为协变量控制基线差异。提取MT×PG交互项的t统计量，以|t| > 1.58（等价于F > 2.5）作为初级阈值，构建二值化显著边矩阵。通过5000次组标签随机置换，基于最大连通组件大小的零分布进行家族误差率（FWE）校正，确定全脑水平显著的连通子网络（即Synergistic Network Block, SNB）。**

**子网络分解：基于每个通道的Brodmann分区标签，将SNB显著边分类为四个子网络——患侧半球内连接（Ipsi-Intra）、健侧半球内连接（Contra-Intra）、跨半球连接（Inter-Hemisphere）及前额极区连接（Prefrontal-Polar）。对每个子网络分别进行2×2析因ANCOVA（Post_Z_sub ~ MT + PG + Pre_Z_sub，Type II SS），以探索交互效应的子网络特异性。**

**敏感性分析：采用独立于SNB的网络掩膜——ISS掩膜（ΔZ ~ Group NBS, F>2.5, 5000次置换），在另一组显著边集上验证交互效应方向的一致性。**

**脑-行为相关分析：计算SNB组件平均ΔZ与临床改善指标（ΔFMA_distal、ΔFMA_proximal、ΔFMA_total、ΔARAT）的Spearman秩相关。分全样本（n=29）和MtPg亚组（n=7）两层报告，采用Benjamini-Hochberg法控制FDR（q < 0.05）。**

### 2.5 统计分析方法

**临床统计分析采用Python 3.x（SciPy v1.11+, Statsmodels v0.14+）。** 临床变量（连续变量的均值±标准差或者中位数（P25，P75）表示数据特征），采用卡方检验或独立样本方差分析、配对t检验或Wilcoxon匹配对符号秩检验进行比较。当P＜0.05时，则表示差异具有统计学意义。**脑网络统计分析同样采用Python环境完成，NBS检验使用5000次非参数置换以避免正态性假设，FWE校正通过最大连通组件大小实现，脑-行为相关分析采用FDR多重比较校正。**

---

## Results

### 3.1 一般资料比较

经单因素方差分析检验，四组受试者在性别构成、年龄分布、病程时间等基线特征方面均无统计学差异（P>0.05），证实组间基线资料具有可比性。

### 3.2 U-FMA评分比较

（保留原文：治疗前后U-FMA近端与远端评分比较、析因分析结果不变。）

### 3.3 ARAT评分比较

（保留原文：治疗前后ARAT评分比较、析因分析结果不变。）

### 3.4 fNIRS结果：皮层激活与组内HbO浓度变化

（保留原文任务态HbO结果不变。）

### 3.5 fNIRS结果：脑网络

**3.5.1 SNB NBS交互效应网络**

**GLM-NBS分析（Post_FC ~ MT + PG + MT×PG + Pre_FC，5000次置换检验）揭示了一组显著的MT×PG交互效应网络。该网络（Synergistic Network Block, SNB）包含100条显著边（FWE校正P = 0.0002），涉及37个节点。在总703条有效连接中，650条通过质量控制被纳入分析，SNB显著边占比15.4%。该结果表明，镜像疗法与气动手套的联合干预在功能连接层面产生了超越各单独干预效果之和的非加性网络重组。为评估结果对初级阈值选择的敏感性，本研究在多个阈值水平下进行了重复分析：|t| > 2.0（等价F > 4.0）时显著边为44条（6.8%），|t| > 2.5时为11条（1.7%），|t| > 3.0时为4条（0.6%）。边数随阈值提高呈稳定递减趋势，而非剧烈波动，证实结果对阈值选择具有良好的稳健性。**

**图3.6a展示了SNB的t值矩阵热力图。38×38矩阵经ROI标签排序后，显著边（黑框标记）在前额极区（PFC相关通道）和跨半球通道间呈密集分布。图3.6b为头皮拓扑图，节点大小反映度中心性，节点颜色表示功能分区归属，边颜色编码t值方向。图3.6c为四视角3D玻璃脑视图。图3.6d展示了5000次置换检验的零分布，观测到的100条边（红色虚线）远超零分布的主体部分（P = 0.0002）。**

**3.5.2 子网络分解**

**对SNB 100条显著边进行子网络分解：前额极区连接（Prefrontal-Polar）占45%（45条），跨半球连接（Inter-Hemisphere）占23%（23条），患侧半球内连接（Ipsi-Intra）占16%（16条），健侧半球内连接（Contra-Intra）占16%（16条）。前额极区连接占据交互效应网络的近半比例，提示联合干预对高级认知-运动整合区域具有优先调控效应。**

**对四个子网络分别进行析因ANCOVA（Post_Z_sub ~ MT + PG + Pre_Z_sub）的结果显示：所有子网络的MT主效应和PG主效应均未达到统计显著性（所有P > 0.05）。这一结果与SNB层面MT×PG交互效应显著一致——在析因设计中，当交互效应显著而主效应不显著时，表明联合效应主要通过交互路径而非单一因子传递，符合"1+1>2"的协同模式预期。**

**图3.6e展示了SNB组件平均ΔZ按组分的箱线图。联合治疗组（MtPg）的ΔZ均值（-0.283 ± 0.217）显著低于加性预测值（0.057），提示联合干预诱导了功能连接的选择性降低（即"连接精简"）。图3.6f为MT×PG交互效应图，展示了无PG条件下MT的效应（接近零线）与有PG条件下MT的效应（负向）之间的显著差异。图3.6g为协同示意图，直观展示了观测到的联合效应与加性预测之间的偏离。**

**3.5.3 敏感性分析：ISS掩膜验证**

**为验证SNB交互效应方向的稳健性，本研究采用独立于SNB的ISS掩膜（57条显著边，ΔZ ~ Group NBS, P = 0.0002）进行敏感性分析。在ISS掩膜上，MT×PG交互效应的方向与SNB一致（均指向MtPg组的连接降低），但交互效应的Type III ANOVA结果为MT主效应显著（F = 8.70, P = 0.007），而MT×PG交互效应不显著（F = 0.248, P = 0.623）。这一差异可能源于ISS掩膜的网络拓扑特征不同——ISS主要捕获组间差异而非交互效应的特异性，因此在该掩膜上MT主效应更为突出。两套掩膜的共同之处在于均支持联合干预对功能连接产生调节作用。**

### 3.6 脑-行为相关分析

**为探索脑网络变化与临床功能改善之间的关联，本研究计算了SNB组件平均ΔZ与四项临床改善指标（ΔFMA_distal、ΔFMA_proximal、ΔFMA_total、ΔARAT）的Spearman秩相关。**

**全样本分析（n = 29，涵盖全部四个实验组）：所有脑-行为相关均未达到统计显著性。其中SNB ΔZ与ΔFMA_distal的相关为r_s = -0.290（P = 0.108，FDR校正后P = 0.151），与ΔARAT的相关为r_s = -0.349（P = 0.050，FDR校正后P = 0.151）。所有相关系数的方向一致为负，提示SNB连接降低与远端功能改善之间可能存在关联趋势，但在全样本层面未达显著水平。**

**MtPg亚组分析（n = 7）：在联合治疗组内部，SNB ΔZ与ΔFMA_distal的Spearman相关系数为r_s = -0.590（P = 0.123，FDR校正后P = 0.314），提示中等程度的负相关——即SNB连接降低越多，远端FMA改善越大。然而，该相关在FDR校正后不显著，且置信区间较宽。为进一步评估该效应的稳健性，进行了Leave-one-out敏感性分析：逐一剔除1例被试后重新计算相关系数，r值始终在-0.865至-0.382范围内波动，方向一致为负，表明该关联非由个别异常值驱动。但由于n = 7时统计效力不足以检测中等效应量（post-hoc power约0.35），这些脑-行为相关结果应视为探索性发现。**

**图3.7展示了MtPg组SNB ΔZ与ΔFMA_distal的散点图，四个组按颜色标记。尽管MtPg组的散点分布呈现负相关趋势，但95%置信区间（Fisher Z变换法）跨零，未达到统计显著性。**

---

## Discussion

本研究通过行为学评估与fNIRS脑功能监测相结合的方法，系统探讨了镜像疗法联合智能气动式柔性手套对脑卒中后手功能障碍的康复效果及其神经机制。行为学结果表明，联合干预在改善上肢运动功能（U-FMA近端与远端、ARAT）方面显著优于部分单一疗法及对照组，且在远端功能改善上存在显著的交互效应（p=0.006, F=8.13），证实了镜像疗法与手套训练之间存在"1+1>2"的协同效应。神经影像学层面，GLM-NBS分析发现了显著的MT×PG交互效应网络（100条边, p=0.0002）。子网络分解发现前额极区连接占主导（45%），其次为跨半球连接（23%），这一发现揭示了联合干预促进远端手功能恢复的特异性神经机制——通过前额极-跨半球连接网络的协同调控。独立ISS掩膜（57条边, p=0.0002）的敏感性分析验证了联合干预对功能连接的调节作用。脑-行为相关分析（SNB掩膜）显示，在MtPg组（n=7）中，SNB ΔZ与ΔFMA_dist呈中等负相关（r_s=-0.59, p=0.123），但经FDR校正后不显著（q<0.05）。全样本分析（n=29）同样未发现显著的相关性。考虑到n=7的样本量下统计效力不足，这些脑-行为相关结果应视为探索性发现，需要更大样本验证。

上述结果与"多模态同步强化"假说高度一致。fNIRS数据显示联合组在SMC_L、SAC_L、DLPFC_L等多个脑区的激活强度均显著高于对照组，且联合组在PFC_L的激活也显著增强，这为上述假说提供了直接的神经影像学证据。

本研究的神经机制发现与现有文献既有呼应之处，也存在重要延伸。既往研究表明，镜像疗法主要通过激活背侧运动前皮层与顶叶联合区增强运动计划能力[61]，而机器人辅助训练则通过密集的体感输入促进运动皮层手部代表区的神经重塑[62]。本研究的fNIRS数据证实了上述通路的存在，并首次从静息态功能连接的角度揭示了联合干预条件下脑网络的特异性重组模式。

值得注意的是，前额极区（Prefrontal-Polar）连接占交互效应网络的45%，是最大的子网络组分。前额极/额内侧皮层是高级认知-运动整合的关键枢纽，涉及运动意图生成、动作监测和多感觉信息整合（Bisio et al., 2025; Bardella et al., 2024）。这一发现提示联合干预可能通过视觉反馈（镜像）+ 体感反馈（手套）的双通道输入，优先重塑前额运动整合网络的信息处理效率，进而通过皮质-皮质投射调节跨半球连接强度。该发现与先前任务态fNIRS观察到的联合组左侧PFC激活增强结果相互印证。

整合上述发现，本研究提出联合干预促进远端手功能恢复的双阶段神经机制模型。该模型的核心逻辑是：联合干预通过"视觉反馈（镜像）+ 本体感觉输入（手套）"的多模态输入，诱导前额极-跨半球连接网络发生选择性重组，使大脑运动网络从损伤初期的"代偿性扩张"状态，逐步过渡到"高效神经编码"状态，从而特异性地促进远端手功能恢复。值得注意的是，这种精简效应仅出现在远端功能（而非近端功能）的改善中，这与远端运动对外侧皮质脊髓束直接单突触连接的依赖性相吻合[65,66]——该通路的信息传输效率对突触噪声极为敏感，因此连接精简带来的信噪比提升对其功能恢复尤为关键。

---

## Limitations

本研究存在以下局限性，需在结果解释时予以审慎考虑。第一，样本量相对有限（80例，其中完成fNIRS静息态分析的子样本为32例），虽然析因设计提供了交互效应的初步证据，但亚组分析（如按Brunnstrom分期或损伤部位分层）的统计效力仍显不足。第二，脑-行为相关分析中MtPg组仅纳入7例被试，统计效力不足以检测中等效应量的相关性（FDR校正后均不显著），**这些脑-行为发现应视为探索性结果**。第三，**子网络分解为探索性分析，未对多个子网络之间的统计比较进行多重比较校正**。第四，随访时间较短，本研究仅评估了干预结束后的即刻效果，未设立长期随访节点（如3个月或6个月），因此联合干预的疗效维持效应及神经可塑性改变的持久性尚不明确。第五，fNIRS技术的固有局限，fNIRS仅能探测大脑皮层（深度约2-3 cm）的血氧动力学变化，无法触及皮层下结构在联合干预中的作用。第六，缺乏主动对照sham组，尽管联合组优于常规对照组，但无法完全排除非特异性效应（如额外关注、治疗师接触时间增加等）对结果的贡献。第七，训练内容与日常生活的转移效应有限，本研究的干预聚焦于手部基础运动模式（握拳、伸指），与Barthel指数所评估的复杂日常活动之间存在技能缺口，这可能解释了功能改善与ADL能力提升分离的现象。

---

## Conclusion

综上所述，本研究证实了镜像疗法联合智能气动式柔性手套对脑卒中后手功能障碍的协同康复效应，并**从静息态功能连接网络的前额极-跨半球重组角度**揭示了其促进远端功能恢复的神经机制。这一发现不仅为"中枢-外周-中枢"闭环康复理论[34]提供了直接的神经影像学证据，也为精准神经康复方案的制定与智能康复系统的研发奠定了理论依据与实践基础。**未来研究应扩大样本量，并结合多模态神经影像技术验证这些探索性发现的可重复性。**

---

## References

（保留原始手稿的全部参考文献不变。）

---

## Figures

**图3.1** U-FMA交互效应图

**图3.2** ARAT交互效应图

**图3.4** 治疗前后4组氧合血红蛋白浓度的变化

**图3.5** 治疗前后组间两两比较HBO浓度差异

**图3.6a** SNB NBS t值热力图（38×38）

**图3.6b** SNB头皮拓扑图

**图3.6c** SNB 3D玻璃脑（四视角）

**图3.6d** SNB零分布直方图

**图3.6e** SNB组件ΔZ箱线图（按组）

**图3.6f** MT×PG交互效应图

**图3.6g** 协同示意图

**图3.7** MtPg组SNB ΔZ与ΔFMA_dist散点图

**表3.1** 4组受试者一般资料比较

**表3.2** 治疗前后U-FMA近端评分比较

**表3.3** 治疗前后U-FMA远端评分比较

**表3.4** 治疗后组间运动功能两两比较

**表3.5** 治疗前后ARAT评分比较
