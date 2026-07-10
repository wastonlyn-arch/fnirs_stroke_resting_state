# MNE-NIRS 0.7.3 源码审计

> 日期: 2026-07-09 | 状态: 进行中

## 架构概览

```
mne_nirs/
├── io/                    # I/O: 只支持 SNIRF
│   ├── snirf/             #   read/write SNIRF (.snirf)
│   └── fold/              #   折叠特异性分析
├── preprocessing/         # 预处理
│   ├── _scalp_coupling_segmented.py  # SCI 头皮耦合指数
│   ├── _peak_power.py                # 峰值功率 QC
│   └── _mayer.py                     # Mayer 波检测
├── signal_enhancement/    # 信号增强（SPN 去除）
│   ├── _short_channel_correction.py  # 短距通道回归（核心）
│   └── _negative_correlation.py      # 负相关增强
├── statistics/            # GLM 统计
│   ├── _glm_level_first.py           # 1st-level GLM
│   ├── _statsmodels.py               # 组分析（混合效应）
│   └── _roi.py                       # ROI 级 GLM
├── channels/              # 通道/ROI 管理
│   ├── _channels.py                  # 通道距离计算
│   ├── _short.py                     # 短距通道检测
│   └── _roi.py                       # ROI 汇总
├── experimental_design/   # 实验设计工具
├── visualisation/         # 可视化（地形图、3D、GLM surface）
├── datasets/              # 内置示例数据集（6个）
├── simulation/            # 仿真数据生成
└── utils/                 # 共享工具
```

## 关键发现（总览）

### 与我们手写管线的差异

| 维度 | MNE-NIRS | 我们的手写管线 | 差异 |
|------|------|------|------|
| I/O | 只支持 SNIRF | 只支持 .nirs (Homer2) | **格式互斥** — 需转换才能用 |
| MBLL | 委托给 MNE `raw._data`（SNIRF 已含 HbO/HbR） | 自己手写 MBLL | MNE-NIRS 不自己做 MBLL，假设数据已转换 |
| 运动校正 | **无！** 依赖外部预处理 | 小波 IQR=1.5 | MNE-NIRS 没有 motion correction |
| 短距通道 | ✅ `short_channel_correction`（PCA/回归） | ❌ | MNE-NIRS 有这个但我们用不了（无硬件） |
| SCI QC | ✅ `scalp_coupling_index`（分段法） | ❌ 未实现 | 可引入做 QC |
| 1st-level GLM | OLS + AR-IRLS（自相关校正）| OLS | MNE-NIRS 额外处理残差自相关 |
| 组分析 | `statsmodels` 混合效应 | 简单减法 | 方向一致但 MNE-NIRS 更完整 |
| 可视化 | ✅ 3D montage, GLM topo, surface projection | ❌ 未实现 | 可借鉴 |
| NBS | ❌ | ✅ 从旧管线 | 互补 — MNE-NIRS 没有网络分析 |

### 核心设计理念

MNE-NIRS 假设的管线：
```
SNIRF (已有 HbO/HbR) → SCI QC → 短距通道回归 → GLM(AR-IRLS) → 混合效应 → 可视化
```

我们的管线：
```
.nirs (raw intensity) → 小波 → MBLL → BPF → FC → NBS
                      → MBLL → BPF → GLM(OLS) → pre-post subtraction
```

**MNE-NIRS 强项**: 统计推断（AR-IRLS, 混合效应）、短距通道、可视化
**我们强项**: 数据格式兼容（.nirs）、运动校正（小波）、网络分析（NBS）
**互不覆盖**: 我们的 MBLL + 小波 + NBS 是 MNE-NIRS 没有的；MNE-NIRS 的 AR-IRLS + 混合效应 + 可视化是我们没有的

## 1. I/O 层：为什么只支持 SNIRF？

### 架构：读/写分离

```
读取: MNE core (mne.io.read_raw_snirf)  ← 不在 MNE-NIRS 里
      └─ RawSNIRF 类: HDF5 → MNE Raw 对象
      
写入: MNE-NIRS (mne_nirs.io.write_raw_snirf)
      └─ MNE Raw 对象 → HDF5 (SNIRF v1.1 spec)
```

MNE-NIRS 只负责**写**，读是 MNE 核心层的事。也支持 NIRx 格式（`mne.io.read_raw_nirx`），但不支持 Homer2 `.nirs`。

### 写入流程（`_snirf.py:19-68`）

```
write_raw_snirf(raw, fname)
  ├── /nirs/formatVersion        → "1.1"
  ├── /nirs/metaDataTags/        → 日期/时间/被试ID/单位/生日/性别
  ├── /nirs/data1/
  │   ├── dataTimeSeries         → raw.get_data().T
  │   ├── time                   → raw.times
  │   └── measurementList1..N/   → 每个通道: sourceIndex, detectorIndex,
  │                                  wavelengthIndex, dataType, dataTypeLabel
  ├── /nirs/probe/
  │   ├── sourceLabels/ detectorLabels/ wavelengths
  │   ├── sourcePos3D/ detectorPos3D   ← raw.info['chs'][i]['loc'][3:9]
  │   └── landmarkPos3D/ landmarkLabels  ← dig points + 可选 10-20 montage
  └── /nirs/stim1..N/            → MNE Annotations → SNIRF stim 表
```

### 通道命名约定

```python
# 正则: r"^S(?P<source>\d+)_D(?P<detector>\d+) (?P<wavelength_type>[\w]+)$"
"S1_D3 760"     # 波长模式
"S1_D3 850"     
"S1_D3 hbo"     # 色团模式（已做 MBLL）
"S1_D3 hbr"
```

和我们的 `.nirs` 格式差异：我们用 `S1-D3`（短横线），MNE 用 `S1_D3`（下划线）；我们没有波长后缀，因为两个波长存在同一通道的不同列里。

### MBLL 在哪里完成？

MNE-NIRS **不包含 MBLL 代码**。MBLL 由 MNE core 的两个方法完成：

```python
raw_od = raw.optical_density()    # raw intensity → ΔOD
raw_hb = raw.haemoglobin()        # ΔOD → HbO/HbR
```

这两个方法内置于 `mne.io.BaseRaw`，所有 fNIRS 格式（SNIRF、NIRx）共享。这和我们的管线不同——我们是手写 MBLL，MNE 把它做成了 Raw 对象的方法。

### 对照我们手写管线

| | MNE-NIRS/MNE core | 我们的 fnirs_loader |
|------|------|------|
| 读 .nirs | ❌ 不支持 | ✅ 专用解析器 |
| 读 SNIRF | ✅ `read_raw_snirf` | ❌ |
| 写 SNIRF | ✅ `write_raw_snirf` | ❌ |
| MBLL | `raw.optical_density()` → `raw.haemoglobin()` | 手写 `to_hbo_hbr()` |
| probe 坐标 | MNE `info['chs'][i]['loc']` 自动管理 | 手写 `ProbeLayout` dataclass |
| stim 标记 | MNE `Annotations` 对象 | `s` 数组直接拿 |
| 消光系数 | MNE 内置（和我们的值一致） | 硬编码 |
| DPF | MNE 用默认值（不做年龄校正） | 硬编码 6.0/5.5 |

### 数据流差异

```
我们的管线:
  .nirs → scipy.loadmat → raw intensity → 小波 → MBLL → HbO/HbR

MNE 管线:
  SNIRF → RawSNIRF → raw intensity → optical_density() → haemoglobin() → HbO/HbR
  (运动校正需要在 MBLL 之前或之后手动调用，MNE-NIRS 不提供)
```

### 依赖关系（2026-07-09 核查）

```
mne-nirs 0.7.3
  ├─ mne 1.9.0         ← 核心 I/O、MBLL (optical_density / haemoglobin)
  ├─ nilearn            ← GLM 引擎
  ├─ h5io               ← SNIRF HDF5 读写
  ├─ seaborn            ← 可视化
  └─ numpy, scipy
```

MNE 1.9.0 完整 fNIRS 能力已就绪，但卡在入口：**MNE 无 Homer2 `.nirs` reader**。数据要进 MNE 只有两条路：
1. MATLAB Homer3 → `.nirs` 转 `.snirf` → MNE `read_raw_snirf`
2. 自写转换器：`fnirs_loader` 输出 → MNE `Raw` 对象（手动构建 `info` + `data`）

## 2. preprocessing 模块：名义是"预处理"，实际全是 QC

### 模块内容

```python
# __init__.py 导出三个函数，全是质量评估，没有一个做信号处理
from ._peak_power import peak_power
from ._scalp_coupling_segmented import scalp_coupling_index_windowed
from ._mayer import quantify_mayer_fooof
```

**MNE-NIRS 没有运动校正、没有带通滤波、没有 MBLL。**
这些在 MNE 生态中的分工是：
- MBLL → MNE core `raw.optical_density()` / `raw.haemoglobin()`
- 滤波 → MNE core `raw.filter()` 或 `mne.filter.filter_data()`
- 运动校正 → **MNE 生态完全缺失**

### 2.1 SCI（头皮耦合指数）— `_scalp_coupling_segmented.py:13-111`

**原理**（Pollonini et al. 2016, PHOEBE）：
```
同一个 S-D 对的两个波长（730/850 nm），如果探头与头皮耦合良好，
两个波长都应该测到心跳信号 → 在心脏频段（0.7-1.5 Hz）相关性高。
如果探头脱开 → 只有噪声 → 相关性低。
```

**算法**：
1. 带通 0.7-1.5 Hz → 提取心跳频段
2. 对每个 10s 窗口，计算同 S-D 对两波长信号的 Pearson r
3. r < threshold（默认 0.1）→ 标记 `BAD_SCI` 注释
4. 只在 OD 数据上运行（fnirs_od）

**对我们的适用性**：
- 理论上可以移植到我们的管线（`fnirs_loader` 出来的 raw intensity 转 OD 后就能跑）
- 但 11 Hz 采样率下 0.7-1.5 Hz 心跳频段刚好在 Nyquist 边缘 → 信号失真 → 相关系数可能被低估
- PHOEBE 原始论文用的是更高采样率的设备 → 11 Hz 可能不够可靠

### 2.2 峰值功率 — `_peak_power.py:14-120`

**原理**：SCI 的变体。不用相关系数，改算两个波长信号互相关后的功率谱峰值。峰值高 → 心跳信号强 → 耦合好。

算法差异：SCI 用 `np.corrcoef`（时域），PeakPower 用 `periodogram(cross-correlate(c1,c2))`（频域）。

### 2.3 Mayer 波量化 — `_mayer.py:12-178`

**原理**（Luke et al. 2021, Donoghue et al. 2020 FOOOF）：
```
Mayer 波 = ~0.1 Hz 的血压振荡，是 fNIRS 中主要的系统性生理噪声之一。
FOOOF 算法把功率谱分解为"非周期性 1/f 背景" + "周期性振荡峰"，
识别最接近 0.1 Hz 的峰 → 提取中心频率、带宽、功率。
```

这个函数不是 QC，是**研究工具**——用于对比不同组/条件的 Mayer 波参数，回答"干预是否改变了自主神经调节"这类次生问题。

### 对照我们手写管线

| 功能 | MNE-NIRS | 我们的手写管线 |
|------|------|------|
| 运动校正 | ❌ | ✅ 小波 db4 IQR=1.5 |
| 带通滤波 | ❌（MNE core `filter`） | ✅ `scipy.signal.butter + filtfilt` |
| MBLL | ❌（MNE core `haemoglobin`） | ✅ 手写 `to_hbo_hbr()` |
| SCI QC | ✅ `scalp_coupling_index_windowed` | ❌ 未实现 |
| 峰值功率 QC | ✅ `peak_power` | ❌ 未实现 |
| Mayer 波 FOOOF | ✅ `quantify_mayer_fooof` | ❌ 未实现 |
| 坏段标记 | ✅ `raw.annotations.append` | ❌ 手动处理 |

**核心差异**：MNE-NIRS 的 preprocessing 是"数据质量评估层"，我们的 preprocessing 是"信号清洗层"。两者不矛盾——理想管线应该两层都有。

## 3. signal_enhancement：两种 SPN 处理方法

### 3.1 短距通道回归 — `_short_channel_correction.py`

**原理**（Scholkmann et al. 2014; Saager et al. 2005）：

```
短距通道 (SD < 1 cm): 光只穿透头皮/颅骨，未到皮层
长距通道 (SD > 1 cm): 光穿透头皮+颅骨+皮层

对每个长距通道：
  alfa = (A_short · A_long) / (A_short · A_short)   # 式 (27)
  A_clean = A_long - alfa × A_short                  # 式 (26)
```

**实现细节**（55-67 行）：
- 要求输入为 OD 数据（`fnirs_od`），不处理 HbO/HbR
- `max_dist` 默认 0.01 m = 1 cm
- `_find_nearest_short`：欧几里得距离找最近短距通道（不一定是同源-检对）
- 回归模型极简：单变量线性回归，无截距，无延迟

**与我们管线的关联**：❌ 无法使用 — 无短距通道硬件。但代码简洁（67 行），如果未来换设备很容易移植。

### 3.2 负相关增强 — `_negative_correlation.py`

**原理**（Cui et al. 2010, *NeuroImage*）：

```
神经激活:  HbO ↑ + HbR ↓ → 负相关（真信号）
系统性噪声: HbO ↑ + HbR ↑ → 正相关（噪声）

利用这一差异增强真信号:
  alpha = std(HbO) / std(HbR)
  HbO_clean = 0.5 × (HbO - alpha × HbR)
  HbR_clean = -(1/alpha) × HbO_clean
```

**实现细节**（59-72 行）：
- 要求输入为 HbO/HbR 数据（post-MBLL）
- 通道必须交替排列（S1_D1 hbo, S1_D1 hbr, S1_D2 hbo, S1_D2 hbr, ...）
- 逐个 S-D 对处理：减均值 → 算 alpha → 增强
- 两个关键假设：a) HbO 和 HbR 通道一一匹配；b) 系统性噪声在两色团中同号

**对我们管线的意义**：

| | GSR | PCA | 负相关增强 |
|------|------|------|------|
| 额外硬件 | ❌ | ❌ | ❌ |
| 核心理念 | 全局信号=噪声 | PC1=最大方差=噪声 | HbO✗HbR=噪声, HbO↘HbR=信号 |
| 可能去掉脑信号 | 会 | 会 | **不太会**（基于生理约束） |
| 需要什么 | 38ch 均值 | 38ch SVD | 匹配的 HbO/HbR 对 |
| 审稿人接受度 | 争议大 | 中等 | 较好（Cui 2010 引用高） |
| 我们能用吗 | 能但审稿人会问 | 能 | ✅ **能用且不需要额外硬件！** |

**这可能是我们 SPN 问题的最佳答案**：基于 HbO/HbR 的反相关生理特性来增强信号，不假设全局=噪声，不需要短距通道，有明确的文献支撑（Cui 2010, *NeuroImage*, ~600+ citations）。

### 3.3 对照我们手写管线

| 功能 | MNE-NIRS | 我们的手写管线 |
|------|------|------|
| 短距通道回归 | ✅ 需硬件 | ❌ 无硬件 |
| 负相关增强 | ✅ `enhance_negative_correlation` | ❌ 未实现 |
| GSR | ❌ | ❌（有意不做） |
| PCA 去噪 | ❌ | 原理已讨论，代码未实现 |

**建议**：`enhance_negative_correlation` 值得移植到我们的管线。放在 MBLL 之后、BPF 之前。作为 SPN 去除的 defend 比 GSR/PCA 更有生理学依据。

## 4. statistics 模块：GLM + 混合效应 + ROI 聚合

### 4.1 架构全景

```
MNE-NIRS statistics/
├── _glm_level_first.py     ← 一阶 GLM（第 7-796 行）
│   ├── RegressionResults   ← 包装 nilearn RegressionResults，暴露 theta/MSE/vcov/per_channel
│   ├── ContrastResults     ← 包装 nilearn Contrast，暴露 effect/p_value/stat/z_score
│   ├── run_glm()           ← 主入口：逐通道调用 nilearn.glm.first_level.run_glm
│   └── _BaseGLM            ← 共享父类：to_dataframe/plot_topo/surface_projection/save/read_glm
├── _statsmodels.py         ← 组分析：statsmodels → tidy DataFrame
│   ├── statsmodels_to_results()    ← MixedLM/RLM/OLS 输出 → df
│   ├── summary_to_dataframe()      ← HTML summary 解析 + 数值精度修复
│   └── expand_summary_dataframe()  ← 索引展开：Group[T.MT]:Time[T.Post] → Group + Time
├── _roi.py                 ← ROI 级聚合
│   └── _glm_region_of_interest()   ← 加权聚合（1/SE 默认）
└── ../utils/_io.py         ← tidy 导出
    ├── glm_to_tidy()                ← RegressionResults/Contrast → tidy df
    └── _tidy_long_to_wide()         ← 宽表转换 + Source/Detector/Chroma 展开
```

**核心洞察：MNE-NIRS 的统计层是 Nilearn + statsmodels 的薄封装。**
不实现自己的 GLM 引擎，而是在两个成熟统计库之上做 fNIRS 适配。

### 4.2 一阶 GLM（`_glm_level_first.py`）

#### run_glm() — 主入口（736-796 行）

```python
def run_glm(raw, design_matrix, noise_model="ar1", bins=0, ...):
    """
    raw: MNE Raw (fnirs_cw_amplitude / fnirs_od / hbo / hbr)
    design_matrix: Nilearn DataFrame (frame_times + regressors)
    noise_model: 'ols' | 'ar1' | 'arN' | 'auto'
    bins: 0=逐通道独立 AR, >0=按 S-D 距离分 bin 共享 AR 系数
    """
    for pick in picks:
        labels, glm_estimates = nilearn_glm(
            raw.get_data(pick).T,          # (time,) → (time, 1)
            design_matrix.values,           # Nilearn DM
            noise_model=noise_model,        # 默认 'ar1'
            bins=bins,                      # 0 = per-channel
        )
        results[pick_name] = glm_estimates[labels[0]]
    return RegressionResults(raw.info, results, design_matrix)
```

**关键设计选择：逐通道 GLM，不是逐体素。** 每个 fNIRS 通道独立跑 GLM，无空间协方差。

#### RegressionResults 类（52-568 行）

包装 `dict[ch_name → nilearn RegressionResults]`，提供：

| 方法 | 返回 | 说明 |
|------|------|------|
| `.theta()` | `(n_ch, n_reg)` | beta 系数 |
| `.MSE()` | `(n_ch, 1)` | 均方误差 |
| `.vcov()` | per-channel | 方差-协方差矩阵 |
| `.t(column)` | per-channel | t 值 |
| `.compute_contrast(con)` | `ContrastResults` | 调用 nilearn `compute_contrast` |
| `.to_dataframe()` | pandas | tidy 格式 |
| `.to_dataframe_region_of_interest(group_by, weighted)` | pandas | ROI 级聚合 |
| `.plot_topo()` | matplotlib | 地形图子图 |
| `.surface_projection()` | PyVista | 3D 脑表面 |
| `.save()` / `.read_glm()` | HDF5 | 序列化 |

#### ContrastResults 类（570-697 行）

包装单个 nilearn `Contrast` 对象：
- `.effect` / `.p_value()` / `.stat()` / `.z_score()` → 标量或逐通道向量
- `.plot_topo()` / `.surface_projection()` → 对比图
- `.to_dataframe()` → tidy

### 4.3 AR-IRLS 噪声模型 — 与我们 OLS 的核心差异

这是 MNE-NIRS 统计模块最重要的设计选择。

**为什么 fNIRS 残差不是白噪声：**

```
测量: Y(t) = β₀ + β₁ × HRF_convolved(t) + ε(t)

OLS 假设: ε(t) 独立同分布 (i.i.d.) → Cov(ε(t), ε(t-k)) = 0

fNIRS 现实: 滤波后 (0.01-0.5 Hz) 的信号有强自相关
            ε(t) 和 ε(t-1) 高度相关 → Cov(ε(t), ε(t-1)) ≠ 0
```

OLS 假设违反 → 标准误被低估 → t 值虚高 → **假阳性率 > 名义 α**。

**AR(1) 模型怎么修：**

```
ε(t) = ρ × ε(t-1) + η(t)    # η(t) 才是白噪声

Nilearn 做法: IRLS (迭代重加权最小二乘)
  1. OLS 初拟合 → 估计残差 ε̂(t)
  2. ε̂(t) = AR(1) 拟合 → 估计 ρ → 构造协方差矩阵 Σ  
  3. 用 Σ 做广义最小二乘 → 新 β
  4. 重复 2-3 直到 β 收敛
```

MNE-NIRS 默认 `noise_model="ar1"`, `bins=0`（逐通道独立估计 ρ）。

其他选项：
- `ols`：纯 OLS（等同我们现在的做法）
- `arN`：用户指定 AR 阶数
- `auto`：自动确定 AR 阶数 = `4 × sample_rate`（11 Hz 时 = 44 阶 — 非常高，可能过度参数化）

**对我们的影响和行动建议：**

我们 `task_glm.py` 用 `np.linalg.lstsq` → 纯 OLS，无自相关校正。
有两种修复路径：
1. **引入 nilearn 一阶 GLM**（推荐）→ 直接用 `noise_model="ar1"`，改动最小
2. 手动在 OLS 后用 Newey-West 或 cluster-based permutation 校正标准误

路径 1 的代码成本：~30 行改动，去掉 `np.linalg.lstsq` 那块，换成 nilearn 调用。Nilearn 已在 MNE-NIRS 依赖中，不需要新装。

### 4.4 Group-level：statsmodels 混合效应（`_statsmodels.py`）

`statsmodels_to_results()`（79-145 行）的核心逻辑：

```
statsmodels mixedlm/rlm/ols 输出
  └─ model.summary() → HTML 表格 → pandas DataFrame
       ├─ 数值精度修复  ← summary 表里数字被截断
       │   ├─ 从 model.pvalues/model.tvalues 直接取精确值
       │   └─ MixedLM: 从 model.cov_params() 重算 SE 和 CI
       └─ expand_summary_dataframe()
            └─ "Group[T.MT]:Time[T.Post]" → Group=MT, Time=Post
```

**MixedLM 特殊处理（105-132 行）：**
statsmodels 的 `MixedLMResultsWrapper.summary()` 有数值精度 bug — HTML 表格中的 SE 和 CI 被截断。MNE-NIRS 绕过它：
- 从 `model.cov_params()` 取精确协方差矩阵
- `SE = sqrt(diag(cov_params))`
- `CI = mu ± qnorm(0.975) × SE`
- 这是代码评审中的好习惯 — 不信任格式化输出，从数值对象直接取。

**这与我们 §3.2d.1 讨论的推荐方案一致：mixed effects 做估计 + permutation 做多重比较校正。**

### 4.5 ROI 聚合：加权 vs 等权（`_roi.py`）

`_glm_region_of_interest()`（51-153 行）：

```python
# 对每个 ROI 的每个 chroma：
thetas = [通道1的beta, 通道2的beta, ...]
ses = [通道1的SE, 通道2的SE, ...]

if weighted:
    weights = 1.0 / ses      # 反比标准误 → 噪声小的通道贡献大
else:
    weights = np.ones(...)   # 等权

theta_roi = Σ(weights × thetas) / Σ(weights)
se_roi = Σ(weights × ses) / Σ(weights)
t_roi = theta_roi / se_roi
p = 2 × t.cdf(-|t_roi|, df)
```

**为什么加权重要：**
- 同一个 ROI 内，通道 1 接触良好（低 SE），通道 3 接触边缘（高 SE）
- 等权均值 → 通道 3 拉低整个 ROI 的 t 值
- 1/SE 加权 → 通道 3 权重自动降低
- 这在卒中患者的实际数据中特别重要（探头贴合常不完美）

**我们目前的做法：`hbo[:, ch_indices].mean(axis=1)` 在时域平均 → 再跑 GLM。**
这和 ROI-level 的 weighting 不一样 — 我们是时域平均后再建模，MNE-NIRS 是先逐通道建模再聚合。两种做法各有利弊。

### 4.6 我们的 `task_glm.py` vs MNE-NIRS：完整对比

| 维度 | task_glm.py (我们) | MNE-NIRS |
|------|------|------|
| GLM 引擎 | `np.linalg.lstsq` | Nilearn `run_glm` |
| **自相关校正** | ❌ OLS 假设 i.i.d. | ✅ **AR(1)-IRLS 默认** |
| 设计矩阵 | 手写 boxcar→SPM HRF 卷积 | Nilearn `make_first_level_design_matrix` |
| HRF | 手写双 Gamma (SPM 参数) | Nilearn 内置 SPM/Glever |
| 漂移建模 | 手写 3 阶多项式 (Legendre) | Nilearn drift 项 |
| Contrast | ❌ 只取 col 0 的 t-stat | ✅ Nilearn `compute_contrast` |
| 输出格式 | 手动 dict/array | tidy DataFrame + HDF5 持久化 |
| ROI 聚合 | 时域均值→GLM | 逐通道 GLM→加权聚合(1/SE) |
| 组分析 | pre-post 减法 + t 检验 | statsmodels mixedlm |
| 可视化 | 手动 matplotlib | plot_topo + surface_projection |

**差距最大的三项：**
1. **自相关校正**（AR-IRLS vs OLS）— 可能系统性高估显著性
2. **Contrast 框架** — 多条件比较没有正式框架
3. **输出生态** — tidy DataFrame + 持久化 vs 手动管理 array

### 4.7 设计模式评注

MNE-NIRS 的统计模块体现了几个好的软件工程实践：

1. **薄封装而非重造轮子**：GLM 委托 Nilearn，组分析委托 statsmodels。代码量小但功能覆盖全。

2. **tidy data 原则**：所有输出统一为 tidy DataFrame（Hadley Wickham 的 "tidy data" 概念 — 每列一个变量，每行一个观测，每种观测表一个类型）。`glm_to_tidy()` 是统一的出口。

3. **防御数值精度问题**：不信任格式化输出（summary HTML），从数值对象直接取精确值。见 `_statsmodels.py` 中对 MixedLM 的处理。

4. **灵活的 ROI 权重**：支持 bool/dict 两种输入，扩展性好。

## 5. visualisation 模块（快速扫描）

### 模块清单

```
visualisation/
├── _plot_3d_montage.py          ← 3D probe 地形图（plotly/mayavi/PyVista）
├── _plot_GLM_projection.py      ← GLM 结果投到脑表面
├── _plot_GLM_topo.py            ← GLM 结果地形图子图
├── _plot_NIRS_topo.py           ← 通用 fNIRS 地形图
└── _plot_oktopus_probe.py       ← probe layout 示意图
```

### 关键函数

**`plot_glm_group_topo`**（`_plot_GLM_topo.py`）：
- 输入：多被试 GLM 结果 list → 组级地形图
- 按 ROI 分组，每个 ROI 区域内的通道用均值
- 支持 `vmin/vmax` 控制色阶

**`plot_glm_surface_projection`**（`_plot_GLM_projection.py`）：
- 投到 FreeSurfer 模板脑表面（`fsaverage`）
- 用的是 mne `SourceEstimate` → `plot_source_estimates`
- 需要 FreeSurfer `fsaverage` 数据（`mne.datasets.fetch_fsaverage()`）

**`plot_3d_montage`**（`_plot_3d_montage.py`）：
- 交互式 3D 探头布局
- 支持 plotly/PyVista/mayavi 三种后端
- 自动按 S-D 距离标色

### 对我们的实用性评估

| 可视化 | 能直接用？ | 条件 |
|------|------|------|
| `plot_glm_topo` | ⚠️ 部分 | 数据进 MNE Raw 对象就行 |
| `plot_glm_surface_projection` | ❌ | 需要 FreeSurfer + 通道→MNI 坐标映射 |
| `plot_3d_montage` | ⚠️ 部分 | 需要 probe 坐标在 MNE info 里 |
| `plot_glm_group_topo` | ✅ | 如果我们用 MNE-NIRS GLM 输出 |
| `plot_nirs_topo` | ⚠️ 部分 | 通用函数，数据格式要求低 |

**建议：发表级可视化，如果可以转换数据格式到 MNE Raw 对象就尽量用 MNE-NIRS。否则自己写 matplotlib。**

## MNE-NIRS 审计总结

### 我们有但没有的

| MNE-NIRS 强项 | 我们应该引入？ | 理由 |
|------|------|------|
| SCI QC | ⚠️ 可选 | 11 Hz 采样率限制可靠性 |
| **负相关增强** | ✅ **推荐** | 不需硬件，有生理学依据，写进方法有 defense |
| **AR-IRLS GLM** | ✅ **推荐** | 自相关校正是当前 OLS 的明显短板 |
| 混合效应组分析 | ✅ 已在计划 | 与 §3.2d.1 讨论一致 |
| 可视化 | ⚠️ 可选 | 数据转换成本 vs 手写 matplotlib |

### MNE-NIRS 没有但我们有的

| 我们强项 | MNE-NIRS 状态 |
|------|------|
| .nirs 格式解析 | ❌ |
| 小波运动校正 | ❌ MNE 生态完全缺失 |
| NBS 脑网络 | ❌ |
| 手写 MBLL (透明度高) | ⚠️ 委托 MNE core，黑盒 |

### 一句话总结

MNE-NIRS 是"统计推断层 + QC 层"的优质薄封装，信噪分离（运动校正）和原始 I/O（.nirs）是我们自己补的。**推荐走"适配层"方案：我们的数据进 MNE Raw 对象，复用 MNE-NIRS 的 AR-IRLS GLM + 混合效应 + 可视化，保留我们自己的小波 + MBLL + NBS。**

## 下一步

- [x] 架构总览
- [x] I/O 层
- [x] preprocessing: 三个 QC 函数
- [x] signal_enhancement: 短距通道回归（无法用）+ 负相关增强（可用！）
- [x] statistics: GLM (AR-IRLS) + 混合效应 + ROI 加权
- [x] visualisation: 3D montage + GLM topo + surface projection
- [ ] **下一工具：Homer3 / SPM-fNIRS / NIRS-KIT / Cedalion**
