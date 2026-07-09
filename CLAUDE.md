# fNIRS Mirror Therapy Project

## 项目概述

镜像疗法（Mirror Therapy, MT）联合智能气动式柔性手套（Pneumatic Glove, PG）训练对脑卒中后手部功能障碍的康复疗效及神经机制研究。基于 fNIRS 功能神经影像，结合临床量表（FMA、ARAT、BI），采用随机对照设计（4组×20人）。

**目标期刊**: JCR Q2/Q3（如 Neurorehabilitation and Neural Repair, Frontiers in Neurology, Brain Topography 等）

**当前阶段**: 修订重投 — 已有静息态分析结果，需补全任务态分析管线，完善文献调研，撰写/修改手稿。

## 项目结构 (BIDS-compliant)

```
fNIRS-mirrorTherapy/
├── CLAUDE.md                         # 本文件
├── bids/                             # BIDS 结构（281 symlinks, 41 subjects）
│   ├── sub-XX/ses-{pre,post}/func/   # BIDS 命名规范
│   └── dataset_description.json
├── sourcedata/                       # → rawData/ + FX/ 的 symlink
├── pheno/
│   ├── participants.tsv              # 被试人口学+分组
│   └── participants.json             # 字段元数据
├── derivatives/
│   └── pipeline-v1/                  # 预处理后数据存放
├── code/
│   ├── 0_data_prep/                  # bids_setup, hemisphere_flip, inventory
│   ├── 1_clinical/                   # factorial_anova
│   ├── 2_preprocessing/              # preprocess_resting, qc_metrics
│   ├── 3_resting_state/              # NBS 分析（从旧 PC 重建）
│   ├── 4_task_state/                 # task_glm（新增）
│   ├── 5_visualization/              # 图表生成
│   ├── utils/                        # fnirs_loader, etc.
│   └── master_analysis.py            # 主编排脚本
├── output/                           # 分析产出
│   ├── figures/ tables/ reports/
│   └── processed/                    # .npy, .csv 中间数据
├── docs/                             # 文献、手稿、笔记
├── rawData/ FX/                      # 原始数据（保留）
└── submission/                       # 上一版投稿材料
```

## 数据清单

### 1. 静息态 fNIRS（rawData/）
- **文件数**: ~1200 (.nirs + .hcx)
- **被试**: ~34 名卒中患者
- **范式**: 10分钟静息态
- **命名**: `1静息态10min_日期_时间_编号_姓名_性别_出生日期_序号.nirs`
- **设备**: NIRS-smartII-3000A, 双波长 730/850 nm, 11 Hz, 38 通道
- **已分析**: GLM-NBS 脑网络（SNB/ISS 掩膜）, 脑-行为相关

### 2. 任务态 fNIRS（FX/TaskState/）
- **治疗前**: 29 subjects, 3 tasks (左手握拳, 右手握拳, 右手被动活动)
- **治疗后**: 29 subjects, 1 task (左手握拳)
- **未分析**: 需要建立完整的任务态分析管线

### 3. 临床数据
- `脑卒中患者康复评估数据-分组对比.xlsx`: 分组对比数据
- `量表最终原始数据2.xlsx`: 原始量表数据
- **指标**: FMA (总分/近端/远端), ARAT, Barthel Index (BI)

### 4. 研究设计
- **4组**: 镜像疗法(MT), 气动手套(PG), 联合(MtPg), 常规对照(Sham)
- **各 n=20** (临床), **n=32 完成静息态 fNIRS** (sham=8, MT=7, PG=9, MtPg=8)
- **干预**: 20min/次, 5次/周, 6周
- **时间点**: 治疗前(T0) / 治疗后(T1)

## 分析管线

### 旧管线（旧PC，需重建）
- 静息态预处理 → tdDR/wavelet → BPF → FC矩阵 → NBS → 脑-行为相关
- 参考: `submission/submission/manuscript_revised.md`

### 新管线（需构建）
1. **0_data_prep**: 读取 .nirs/.hcx 文件, 整合临床数据, 被试匹配
2. **1_clinical**: 临床结果复现 — ANOVA/LMM, 交互效应, 效应量
3. **2_preprocessing**: 统一预处理管线（wavelet motion correction, BPF 0.01-0.1 Hz, HbO/HbR 提取）
4. **3_resting_state**: 复现 NBS 结果（GLM-NBS, SNB/ISS, 5000 permutations）
5. **4_task_state**: **新建** — Block-design GLM 分析, ROI 激活, 组间比较, 偏侧化指数
6. **5_visualization**: 发表级图表（脑地形图、NBS 网络图、脑-行为散点图）

## 关键技术栈

- **Python 3.13** (miniconda)
- **核心**: numpy, scipy, pandas, statsmodels, pingouin
- **fNIRS**: mne, mne-nirs, nilearn, nistats
- **可视化**: matplotlib, seaborn
- **统计**: scipy.stats, statsmodels, pingouin (ANOVA/ANCOVA/LMM)

## 当前优先级

1. **任务态分析管线** — 从零构建，这是新增贡献的核心
2. **文献调研** — 任务态 fNIRS + 镜像疗法 + 卒中康复，找 gap
3. **手稿修改** — 在静息态基础上补充任务态发现

## 工作约定

- 分析代码用 Python scripts + Jupyter notebooks（探索性分析用 notebook，最终管线用 .py）
- 所有输出图表放在 `output/figures/`，用描述性文件名
- 预处理后的数据存到 `output/processed/`
- 文献调研用 Zotero 管理，通过 MCP 工具检索
- 手稿用 Markdown 撰写，图表引用 `output/figures/` 路径
