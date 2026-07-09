# Submission Package — 投稿文件包

> 基于最新 SNB 分析结果 (100条显著边, p=0.0002)
> 更新日期: 2026-06-22

---

## 文件说明

### 核心文件

| 文件 | 说明 |
|------|------|
| `manuscript_revision.md` | **审阅修订稿**（学术段落体，可直接投稿参考） |
| `manuscript_revised.md` | 技术修订稿（含修订标记和数据对照，供内部审查） |
| `manuscript.docx` | 原始 Word 手稿（保留作为参考） |

---

## 文件说明

### 核心文件

| 文件 | 说明 |
|------|------|
| `manuscript_revised.md` | **修订稿 Markdown**，含 `~~删除线~~` 和 `**粗体**` 修订标记，可直接查看所有修改 |
| `manuscript.docx` | 原始 Word 手稿（保留作为参考） |

### 新管线图片 (figures_new/)

| 类别 | 图片 | 手稿引用位置 |
|------|------|------------|
| **SNB NBS 网络** | `MirrorTherapy_SNB_NBS_Heatmap_*.png` | §3.6.1 图 3.6a |
| | `MirrorTherapy_SNB_NBS_Topoplot_*.png` | §3.6.1 图 3.6b |
| | `MirrorTherapy_SNB_NBS_3DGlass_*.png` | §3.6.1 图 3.6c |
| | `MirrorTherapy_SNB_NBS_NullDist_*.png` | §3.6.1 图 3.6d |
| | `MirrorTherapy_SNB_NBS_3D_interactive.html` | 交互式 3D 连接组 |
| **SNB ΔZ 效应** | `MirrorTherapy_SNB_deltaZ_Boxplot_*.png` | §3.6.2 图 3.6e |
| | `MirrorTherapy_SNB_deltaZ_Interaction_*.png` | §3.6.2 图 3.6f |
| | `MirrorTherapy_SNB_deltaZ_Synergy_*.png` | §3.6.2 图 3.6g |
| **ISS 敏感性** | `MirrorTherapy_ISS_deltaZ_Boxplot_*.png` | §3.6.4 图 3.6h |
| | `MirrorTherapy_ISS_deltaZ_Interaction_*.png` | §3.6.4 图 3.6i |
| **脑-行为** | `MirrorTherapy_BrainBehavior_*_vs_*.png` (24张) | §3.7 补充材料 |
| | `mtpg_scatter_snb_dz_vs_dFMA_dist.png` | §3.7.2 图 3.7 |
| **子网络** | `nbs_subnetwork_bars.png` | §3.7 图 3.8 |

### PDF 报告

| 文件 | 说明 |
|------|------|
| `MirrorTherapy_FinalReport.pdf` | 完整主分析报告 (SNB+NBS+ΔZ ANOVA) |
| `MirrorTherapy_SupplementaryReport.pdf` | 补充材料 (ISS/CATE/脑-行为散点图/子网络分解) |
| `MirrorTherapy_StatisticalProcess.pdf` | 统计分析过程说明 (中英双语) |
| `MirrorTherapy_DataManifest.pdf` | 输出文件清单 |

---

## 修订摘要

| # | 修改内容 | 位置 |
|---|---------|------|
| 1 | 运动校正: ~~tdDR~~ → **wavelet (db4, IQR=1.5)** | §2.5.2 |
| 2 | 新增静息态预处理管线描述 | §2.5.2 |
| 3 | 新增 Brain Network 方法学 (NBS/ANCOVA/脑-行为) | §2.6 |
| 4 | 统计软件: ~~SPSS 26.0~~ → **Python 3.x (SciPy/Statsmodels)** | §2.7 |
| 5 | Brain Network 结果: ~~220条边~~ → **100条边, p=0.0002** | §3.6 |
| 6 | 脑-行为: ~~r=-0.71 (ISS掩膜)~~ → **r=-0.59 (SNB掩膜), FDR ns** | §3.7 |
| 7 | Discussion 修正: 三路径→两条验证路径 | §4 |
| 8 | 新增前额极区发现的讨论 | §4 |
| 9 | 格式: `~~删除线~~`=删除, `**粗体**`=新增 | 全文 |

---

## 代码引用

- 预处理管线: [`fNIRS_mirror_therapy/code/3_fnirs_analysis/pipeline/`](../code/3_fnirs_analysis/pipeline/)
- NBS 分析: [`fNIRS_mirror_therapy/code/final_analysis/`](../code/final_analysis/)
- 配置: [`fNIRS_mirror_therapy/config/`](../config/)
