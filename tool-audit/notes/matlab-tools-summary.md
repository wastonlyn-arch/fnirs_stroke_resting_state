# Homer3 + SPM-fNIRS + NIRS-KIT 源码审计（概要）

> 日期: 2026-07-09 | 状态: 概要（无可运行 MATLAB license）

## Homer3（482 .m 文件，GitHub BUNPC）

### 架构

```
Homer3/
├── Homer3.m                   ← 主入口（GUI + 脚本模式）
├── FuncRegistry/              ← 插件式函数注册系统
│   └── UserFunctions/
│       ├── hmrR_*             ← Run-level 处理（逐被试逐 run）
│       ├── hmrE_*             ← Epoch-level
│       ├── hmrG_*             ← Group-level (第 2 层)
│       └── hmrS_*             ← Session-level
├── DataFiles/                 ← DataClass: SNIRF .snirf / .nirs 读写
├── ProcStreamEditGUI/         ← GUI 处理流编辑器
├── PlotProbe2GUI/             ← 探头布局可视化
└── Utils/                     ← 共享工具 (wavelet, filter, stats)
```

### 核心设计模式

1. **DataClass 对象**：MATLAB class，封装 SNIRF 兼容的数据结构
   - `.d` = 数据矩阵 (time × channels)
   - `.t` = 时间向量
   - `.ml` = MeasurementList (source/detector/wavelength 元数据)
   - `.s` = stim 标记向量 (time × conditions)

2. **处理流（Processing Stream）**：GUI 或脚本定义的函数链
   ```
   hmrR_Intensity2OD → hmrR_MotionArtifact → hmrR_MotionCorrectWavelet
   → hmrR_BandpassFilt → hmrR_BlockAvg → hmrG_SubjAvg → hmrG_t_paired_contrast
   ```

3. **命名约定**：
   - `R` = Run-level (单文件)
   - `E` = Epoch-level
   - `G` = Group-level
   - `S` = Session-level

### 关键函数与 Cedalion 对比

| 函数 | 行数 | Cedalion 移植 | 移植质量 |
|------|------|------|------|
| `hmrR_Intensity2OD` | 32 | `int2od()` | 1:1（但有差异：Homer3 用 `mean(abs(d))`，Cedalion 用 `mean(d)`）|
| `hmrR_BandpassFilt` | 91 | `freq_filter()` | API 不同 |
| `hmrR_MotionArtifact` | 142 | `id_motion()` | 1:1 |
| `hmrR_MotionCorrectWavelet` | 174 | `motion_correct_wavelet()` | 1:1 |
| `hmrR_MotionCorrectPCA` | 166 | `motion_correct_PCA()` | 1:1 |
| `hmrR_MotionCorrectSpline` | ? | `motion_correct_spline()` | 1:1 |
| `hmrR_BlockAvg` | 170 | 无直接对应 | — |
| `hmrG_SubjAvg` | ? | 无直接对应 | — |
| `hmrG_t_paired_contrast` | ? | 无直接对应 | — |

### Homer3 的设计优势

1. **GUI 降低门槛**：非编程用户也能用 ProcStreamEditGUI 构建处理流
2. **.nirs 原生支持**：我们的数据格式就是 Homer2 的 `.nirs`，Homer3 向后兼容
3. **完整管线**：I/O → 预处理 → epoch → 个体统计 → 组统计，一站式
4. **成熟生态**：20+ 年历史，论文引用最多的 fNIRS 工具

### Homer3 的设计劣势

1. MATLAB license 依赖：¥ 数万/年
2. 无程序化接口友好：GUI 导向，批量处理/自动化困难
3. 统计方法基础：BlockAvg + t-test，无 GLM（Homer3 的 GLM 是实验性的 `hmrR_GLM_sim.m`，仅仿真模式）
4. 无网络分析（NBS）：完全缺失
5. 无混合效应模型：只有 paired t-test

---

## NIRS-KIT（无法获取源码）

- GitHub: `nirstorm/nirs-kit`，专为任务态和静息态 fNIRS 设计的 MATLAB 工具箱
- 特点：内置 NBS、度中心性、FDR 校正、多种 FC 方法
- 状态：❌ HTTPS clone 失败（需认证），安装包可通过官网下载
- 无法审计

---

## SPM-fNIRS（未获取）

- SPM12 的 fNIRS 插件，基于 MATLAB
- 特点：GLM mass-univariate 方法（像 fMRI）
- 使用 SPM 的 canonical HRF + 设计矩阵框架
- 无独立发布 — 是 SPM12 的 Toolbox
- 无法审计（需 SPM12 源码 + MATLAB）

---

## MATLAB 工具 vs Python 工具：总结

| 维度 | Homer3 | MNE-NIRS | Cedalion |
|------|------|------|------|
| 语言 | MATLAB | Python | Python |
| 运动校正 | ✅ Spline/PCA/Wavelet | ❌ | ✅ 6 种（移植自 Homer3）|
| GLM | ❌（实验性）| ✅ Nilearn AR-IRLS | ✅ 自实现 AR-IRLS |
| 组分析 | 配对 t 检验 | 混合效应 (statsmodels) | ❌（需自行拼接）|
| NBS | ❌ | ❌ | ❌ |
| .nirs 读取 | ✅ 原生 | ❌ | ❌ |
| 短距通道 | ✅ | ✅ | ✅（距离阈值分离后手动回归）|
| GUI | ✅ | ❌ | ❌ |
| 自动管线 | ⚠️ GUI 流 | ✅ 脚本 | ✅ 脚本 |
| MBLL | `mean(abs(d))` | MNE core | `mean(d)`（不含 abs）|

**MBLL 基线差异**：Homer3 用 `mean(abs(d))`，Cedalion 用 `mean(d)`。这个差异在原始强度数据为负值时（AD 转换器可能输出 signed int）会产生不同的 OD 值。我们的 `.nirs` 数据用的是 `mean(d)`（`fnirs_loader` 第 137 行：`baseline = np.mean(intensity, axis=-1)`），与 Cedalion 一致。
