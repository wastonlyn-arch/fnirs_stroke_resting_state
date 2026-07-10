# Cedalion 25.1.0 源码审计

> 日期: 2026-07-09 | 状态: 完成

## 项目背景

Cedalion 是 TU Berlin 的 fNIRS Python 工具箱（https://github.com/ibs-tu-berlin/cedalion）。不同于 MNE-NIRS 的"薄封装"策略，Cedalion 采用**自实现 + xarray 生态**的路线。

## 架构概览

```
cedalion/  (74 .py files)
├── __init__.py              ← pint Quantity, units → xarray DataArray
├── dataclasses/             ← Recording, PointCloud, NDTimeSeries schema
├── io/                      ← SNIRF, BIDS, anatomy, probe_geometry
├── nirs.py                  ← int2od → od2conc (MBLL 全链)
├── sigproc/                 ← 信号处理核心
│   ├── motion_correct.py    ←   ★ 6 种运动校正算法
│   ├── quality.py           ←   ★ SCI, PSP, GVTD, SNR, id_motion
│   ├── frequency.py         ←   滤波
│   ├── epochs.py            ←   epoch 切分
│   └── physio.py            ←   生理噪声
├── sigdecomp/               ← ICA_EBM (Entropy Bound Minimization)
├── math/                    ← AR-IRLS GLM, AR model, stats helpers
├── imagereco/               ← 图像重建 (forward model, solver)
├── models/                  ← (空目录)
├── sim/                     ← 仿真数据生成 (HRF, artifact)
└── vis/                     ← probe plot, time_series
```

### 核心设计理念

**数据模型**: `xarray.DataArray` + `pint` 物理单位
```
每个 DataArray 携带:
  - dims: channel, wavelength (或 chromo), time
  - coords: source, detector, wavelength values
  - attrs: pint units → 单位自动追踪/转换
```

与 MNE-NIRS 的 `mne.io.Raw` 对象完全不同。xarray 更接近 numpy，但保留了 labeled dimensions。

## 与 MNE-NIRS 的关键差异

| 维度 | Cedalion | MNE-NIRS |
|------|------|------|
| 数据模型 | xarray.DataArray + pint | mne.io.Raw |
| 设计哲学 | 自实现核心算法 | 薄封装 nilearn/statsmodels |
| I/O | SNIRF (h5py 直接) | SNIRF (MNE core + h5io) |
| MBLL | `int2od` → `od2conc` (Homer3 Prahl) | MNE `optical_density()` → `haemoglobin()` |
| 运动校正 | **6 种算法** + QC | **无** |
| GLM | **自实现 AR-IRLS + RLM** | nilearn AR(p)-IRLS |
| 信号分解 | **ICA-EBM** | 无 |
| 图像重建 | ✅ forward model + solver | ❌ |
| 依赖 | xarray, pint, statsmodels, pywt | mne, nilearn, h5io, seaborn |
| 成熟度 | 学术项目 (TU Berlin) | MNE 官方子项目 |
| 文档 | Sphinx + 引用标注 | MNE 风格 |

---

## 1. 数据模型：xarray + pint = 单位安全的 fNIRS

### 核心类型

```python
# cedalion/typing.py
NDTimeSeries = xr.DataArray  # dims: (channel, wavelength, time)
                             # coords: source, detector
                             # attrs: pint units

LabeledPointCloud = xr.DataArray  # dims: (label, crs)
                                  # coords: 3D positions
```

### 单位系统

```python
from cedalion import units, Quantity

# 所有物理量都带单位
fs = sampling_rate(data)       # → Quantity: 11 Hz
window = 10 * units.s          # → Quantity: 10 s
cardiac = (0.5, 2.5) * units.Hz

# 自动单位转换
dists = dists.pint.to("mm")    # cm → mm 自动转换
conc = conc.pint.to("micromolar")
```

**优势**: 杜绝单位错误（比如忘记把 cm 转 mm 导致 MBLL 结果差一个数量级）。

**劣势**: 学习曲线高，调试困难（pint 的错误信息不友好）。

---

## 2. 运动校正：Cedalion 的最强项

### 算法清单

| 函数 | 算法 | 来源 | 行数 |
|------|------|------|------|
| `motion_correct_spline` | Spline 插值 | Homer3 `hmrR_tInc_baselineshift_Ch_Nirs` | ~168 |
| `motion_correct_splineSG` | Spline + Savitzky-Golay | Homer3 + SG 平滑 | ~60 |
| `motion_correct_PCA` | PCA 去运动伪影 | Homer3 `hmrR_MotionCorrectPCA` | ~135 |
| `motion_correct_PCA_recurse` | 递归 PCA | Homer3 迭代法 | ~60 |
| `tddr` | 时间导数分布修复 | Fishburn 2019 (NeuroImage) | ~100 |
| `motion_correct_wavelet` | 小波阈值去尖峰 | Molavi 2012, Homer3 | ~70 |

### 与我们的管线对比

**`motion_correct_wavelet`** 是和我们管线最接近的：
```
我们的实现:
  pywt.swt → IQR 阈值 → pywt.iswt

Cedalion:
  pad_to_power_2 → dc removal → normalize (MAD) → pywt.swt
  → process_coefficients (block-IQR) → pywt.iswt → denormalize → trim
```

差异：
1. Cedalion 加了 **MAD 归一化**（Homer3 `NormalizationNoise` 函数）— 这个我们没做
2. Cedalion 的 IQR 是按 block 分的（`process_coefficients` 里分块），我们是全局 IQR
3. 都用了 `pywt.swt`（平稳小波变换），不是 DWT

**TDDR** 是完全不同的思路：
- 不检测/修复特定伪影段
- 假设运动伪影是"大波动"→ 用 Tukey biweight 迭代降权
- 优点：无参数调优（tune=4.685 = 95% 统计效率）
- 缺点：可能把真正的 strong activation 也当成伪影

### 运动检测（`quality.py`）

| 函数 | 算法 | 来源 |
|------|------|------|
| `id_motion` | stdev/amp 阈值检测 | Homer3 `hmR_MotionArtifact` |
| `id_motion_refine` | 按通道/全局聚合 | Homer3 |
| `detect_outliers` | std + grad IQR | Homer3 `hmR_tInc_baselineshift` |
| `detect_baselineshift` | 段 delta 阈值 | Homer3 |

---

## 3. GLM：自实现 AR-IRLS（不用 Nilearn）

### 核心函数：`ar_irls_GLM` (`math/ar_irls.py:10-113`)

```
算法流程:
  1. OLS 初拟合 (statsmodels RLM with HuberT)
  2. 残差 → BIC AR 模型选择 (pmax=40)
  3. AR 系数 → prewhitening filter
  4. 对 y 和 X 都应用 prewhitening filter
  5. 重拟合 RLM → 新残差
  6. 重复 2-5 (4 次迭代)
```

**与 MNE-NIRS (nilearn) 的关键区别：**

| | Cedalion | MNE-NIRS |
|------|------|------|
| 回归方法 | **Robust LM (HuberT)** | OLS/GLS |
| AR 估计 | BIC 选阶 (pmax=40) | 用户指定 ar1/arN |
| Prewhitening | **对 y 和 X 都滤波** | nilearn 内部处理 |
| 异常值处理 | **HuberT biweight** 内建 | 无 |
| 迭代次数 | 固定 4 次 | IRLS 至收敛 |
| 实现方式 | 自己写的 | nilearn 黑盒 |

**重要提示**（源码注释）：
```
"DO NOT preprocess your data with a low pass filter.
The algorithm is trying to transform the residual to create a
white spectrum. If part of the spectrum is missing due to low pass
filtering, the AR coefficients will be unstable."
```

→ 如果用 Cedalion 的 AR-IRLS，**不能用我们现在的 BPF (0.01-0.5 Hz)**。应该把低频漂移建模进设计矩阵（Legendre polynomials），让 AR 模型自己处理残差的白化。

这跟 Nilearn 的 AR-IRLS 有同样要求 — 但我们在 §4.3 MNE-NIRS 审计中没提这一点。这是个重要的修正。

### AR 模型选择：`bic_arfit` (`math/ar_model.py:7-31`)

```python
for p in range(pmax + 1):
    model = sm.tsa.AutoReg(dd, lags=p).fit()
    if model.bic > last_bic:   # BIC 开始上升 → 前一个 p 最优
        break
return sm.tsa.AutoReg(dd, lags=p-1).fit()
```

BIC 自动选择最优 AR 阶数 — 不需要用户指定 `ar1` vs `arN`。

---

## 4. 信号质量评估（`sigproc/quality.py`）

Cedalion 的 QC 模块比 MNE-NIRS 丰富得多：

| 指标 | 函数 | 原理 | 对应 MNE-NIRS |
|------|------|------|------|
| SCI | `sci()` | 心跳频段两波长相关 | `scalp_coupling_index_windowed` |
| PSP | `psp()` | 互相关功率谱峰值 | `peak_power` |
| GVTD | `gvtd()` | 全局时间导数 RMS | ❌ 无 |
| SNR | `snr()` | mean/std | ❌ 无 |
| 振幅 | `mean_amp()` | 均值范围 | ❌ 无 |
| S-D 距离 | `sd_dist()` | 距离范围 | ❌ (channels 模块有类似) |
| 通道修剪 | `prune_ch()` | 多指标组合 | ❌ 无 |

GVTD（Global Variance of Temporal Derivatives）是 Sherafati 2020 提出的指标 — 用所有通道的时间导数的 RMS 来判断是否有全局运动。比逐通道检测更稳健。

---

## 5. 信号分解：ICA-EBM

`ICA_EBM()` (`sigdecomp/ICA_EBM.py`) — 955 行的独立组件分析实现。

**与 FastICA 的区别：**
- FastICA：最大化非高斯性（negentropy 近似用 `log cosh` 或 `exp(-x²/2)`）
- ICA-EBM：使用 **熵界最小化**（Entropy Bound Minimization）— 4 种非线性函数（x⁴, |x|/(1+|x|), x|x|/(10+|x|), x/(1+x²)）同时评估，选择 tightest entropy bound

EBM 据说对 fNIRS 中的多种分布类型（超高斯/亚高斯混合）更稳健。但代码量大（~900 行 vs FastICA ~100 行），计算更重。

### 算法阶段：
1. **Part 0**: FastICA 初始化
2. **Part 1**: 正交 ICA（随机梯度）
3. **Part 2**: 鞍点检测 → 旋转成对组件
4. **Part 3**: 如有鞍点 → 正交 ICA 精炼
5. **Part 4**: 非正交 ICA（固定步长梯度）

---

## 6. MBLL 转换链（`nirs.py`）

完全自实现，不使用 MNE core：

```python
beer_lambert(amplitudes, geo3d, dpf):
  od = int2od(amplitudes)       # I → -ln(I/I₀)
  conc = od2conc(od, geo3d, dpf)# ΔOD → Δ[HbO]/Δ[HbR]
    └─ E = get_extinction_coefficients("prahl", wavelengths)
       Einv = pinv(E)
       dists = channel_distances(od, geo3d)  # S-D 距离
       conc = Einv @ (od / (dists * dpf))    # MBLL
```

消光系数来自 **Prahl 1998**（和 Homer2/3 一致）。与 MNE core 的区别：MNE 用 Moaveni 或 Gratzer 数据，值近似但不完全相同。

DPF 作为独立参数传入 — 不做默认假设，这比 MNE 更好（MNE 默认 DPF=1 或固定值）。

---

## 7. 与我们管线的对比

### Cedalion 有但我们没有的

| 能力 | 可用性 | 行动建议 |
|------|------|------|
| TDDR 运动校正 | ✅ | 可以作为小波外的替代方法 |
| Spline + PCA 运动校正 | ✅ | 多算法比较: "哪个运动校正对结果最稳健?" |
| GVTD 全局运动检测 | ✅ | 比我们现在的逐段人工检查好 |
| BIC AR-IRLS GLM | ✅ | **不同思路**: 不用 nilearn，用自实现 |
| ICA-EBM | ⚠️ 代码量大 | 审稿人可能要求 ICA 去噪时可用 |
| pint 单位安全 | ⚠️ 需改造管线 | 重构成本高，但论文发表级 |
| 图像重建 (imagereco) | ✅ | 可以把通道级激活投到脑表面 |

### Cedalion 没有但我们有的

| 能力 | Cedalion 状态 |
|------|------|
| .nirs 格式解析 | ❌ 只支持 SNIRF |
| NBS 脑网络分析 | ❌ |
| 负相关增强 SPN | ❌ 但可移植 |
| 小波 IQR=1.5 | ✅ 有（几乎是 1:1 移植） |

### Cedalion 也没有的（和 MNE-NIRS 一样）

| 缺失能力 | 说明 |
|------|------|
| Homer2 .nirs 直接读取 | 只支持 SNIRF |
| 完整 GLM pipeline | 有 AR-IRLS GLM 函数但无完整统计框架 |
| NBS/图论分析 | 无网络分析模块 |

---

## 8. Cedalion vs MNE-NIRS：选型对比

| 考量 | Cedalion | MNE-NIRS |
|------|------|------|
| 学习曲线 | 陡峭 (xarray+pint) | 中等 (MNE 生态) |
| 运动校正 | **完胜** (6 种算法) | 无 |
| GLM 统计 | 自实现 AR-IRLS | Nilearn 成熟框架 |
| 代码质量 | 学术代码，有 TODO/FIXME | 工程级，测试覆盖好 |
| 社区支持 | 小 (TU Berlin) | MNE 社区大 |
| 长期维护风险 | 高（学术组，可能停更） | 低（MNE 官方维护） |
| fNIRS 领域采用 | 增长中 | 已成为 Python 标准 |
| 对我们项目的价值 | 运动校正方法库 + ICA | 统计推断管线 |

## Cedalion 审计结论

Cedalion 是 **运动校正方法的最佳参考实现** — 6 种算法的代码可以作为我们的"算法对比"章节的基线。其 AR-IRLS GLM 提供了不同于 nilearn 的思路（Robust + prewhitening vs GLS），但没有完整的统计输出框架。**对于我们的管线，推荐：运动校正算法参考 Cedalion，统计推断用 MNE-NIRS/nilearn。**
