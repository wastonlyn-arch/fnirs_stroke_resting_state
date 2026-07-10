# fNIRS 分析工具源码审计

> 目的：系统比较主流 fNIRS 分析工具的架构、算法实现和设计取舍，为手写管线提供参考基准。
> 日期: 2026-07-09 | 状态: ✅ 完成

## 工具审计进度

| 工具 | 语言 | 可运行 | 审计状态 | 笔记 |
|------|------|------|------|------|
| **MNE-NIRS** 0.7.3 | Python | ✅ | ✅ 完成 | `notes/mne-nirs.md` |
| **Cedalion** 25.1.0 | Python | ✅ | ✅ 完成 | `notes/cedalion.md` |
| **Homer3** | MATLAB | ❌ | ✅ 概要 | `notes/matlab-tools-summary.md` |
| **SPM-fNIRS** | MATLAB | ❌ | ❌ 无法获取 | — |
| **NIRS-KIT** | MATLAB | ❌ | ❌ 无法获取 | — |
| **手写管线** | Python | ✅ | 参考基准 | `code/` |

## 核心发现：三工具对比

| 维度 | MNE-NIRS | Cedalion | Homer3 |
|------|------|------|------|
| **设计哲学** | 薄封装 (nilearn/statsmodels) | 自实现 + xarray/pint | GUI 驱动的完整管线 |
| **运动校正** | ❌ 完全缺失 | ✅ **6 种** (移植自 Homer3) | ✅ Spline/PCA/Wavelet |
| **GLM** | Nilearn AR-IRLS | 自实现 AR-IRLS (RLM) | ❌ 实验性 |
| **组分析** | **statsmodels 混合效应** | ❌ | 配对 t 检验 |
| **SPN 去除** | ✅ 短距通道 + 负相关增强 | ⚠️ 距离阈值分离 | ✅ 短距通道 |
| **QC** | SCI + PeakPower + Mayer | **SCI + PSP + GVTD + SNR + amp** | 基础 |
| **.nirs 读取** | ❌ | ❌ | ✅ 原生 |
| **NBS** | ❌ | ❌ | ❌ |
| **ICA** | ❌ | ✅ ICA-EBM | ❌ |
| **可视化** | ✅ 3D + topo + surface | ⚠️ 基础 | ✅ 丰富 |

## 对我们管线的关键建议

### 1. 运动校正：参考 Cedalion
- Cedalion 的 6 种算法是最佳参考实现（1:1 移植自 Homer3）
- 建议在方法学中报告多种运动校正方法的结果比较（稳健性分析）
- 我们的 `motion_correct_wavelet` 与 Cedalion/Homer3 一致，但缺少 MAD 归一化

### 2. GLM 统计：推荐 MNE-NIRS (Nilearn)
- AR-IRLS 自相关校正是我们必须加的（当前 OLS 有假阳性风险）
- Nilearn 比 Cedalion 自实现的 AR-IRLS 更成熟
- 组分析推荐 statsmodels 混合效应（与 §3.2d.1 讨论一致）

### 3. SPN 去除：用负相关增强
- 不依赖短距通道硬件（我们没有）
- Cui 2010 方法有明确生理学依据
- 代码可以从 `mne_nirs.signal_enhancement.enhance_negative_correlation` 移植

### 4. 数据适配层
- 需要 `.nirs` → MNE Raw 或 `.nirs` → xarray 的转换器
- 这是复用 MNE-NIRS/Cedalion 功能的前提

## 输出

- `notes/mne-nirs.md` — MNE-NIRS 0.7.3 完整审计
- `notes/cedalion.md` — Cedalion 25.1.0 完整审计
- `notes/matlab-tools-summary.md` — Homer3/SPM-fNIRS/NIRS-KIT 概要
- `src/` — 源码副本 (mne-nirs, cedalion, homer3)
- `README.md` — 本文件（最终对比矩阵）
