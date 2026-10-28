# 原三角曲面流体控制体实验库

本库将v58已测完整MATLAB模块整理为独立目录，沿用原X轴/SI和全部表面网格。默认工程链路仍使用既有区域方法；本库没有自动取得生产物理解覆盖资格。常gamma/mu/k层流近似、未收敛派生壁热及数值稳态的身份保持明确。

从项目根目录设置路径：

```matlab
addpath('internalFlow');
setupCutCellFlow;
```

- `partitionOrientedTriangles`：按原三角面与Cartesian格裁墙片，保留原面编号、面积/矩和共面唯一几何归属。
- `buildOriginalSurfaceCutCells`：非凸原实体的真实流体体积、质心、唯一共享开放面和真正墙所有者；逐CV向量/九矩闭合。
- `agglomerateFluidControlVolumes`：只沿正开放共享面有界合并，严格保持ΣV与ΣVc及全部墙片。无可行邻居保留原几何并明确报告。
- `cutCellNavierStokesResidual`：flat N×5五量、共享HLLC、Newton应力与Fourier曲面热流。固体NaN；墙温给定、p/q当前态派生。
- `cutCellBoundaryStencilState`、`buildCutCellDiffusionStencil`：与空间内核同序的primitive/特征分支及批量稀疏扩散模板。绝对行和仅是刚度估计，冻结谱不证明完整非线性NS稳定。
- `advanceCutCellNavierStokes`：变体积共享全局SSP-RK2，候选失败整步撤回，五量两阶段边界账；目标时间/工作上限/持续残差稳态分别报告。
- `mapCutCellWallFields`：按Σq*A与Σp*n_fluid*A积分回原三角面；每变量独立未知，最小未知正片也不会被面积容差忽略。给定Tw不标成解。

上述8个模块算法源自已测v58快照；除去已失效的草稿目录addpath，计算语句保持。迁移与函数解析测试位于`tests/testCutCellLibrary.m`，记录源/目标SHA和实际数值对照，不重复基础六步长回归。

v60另外提供4个标量空间实验模块：`buildCutCellVolumeMoments`在MATLAB内部计算真实二阶矩；`buildQuadraticTemperatureStencil`、`evaluateQuadraticTemperature`及`quadraticTemperatureFourierFlux`在真实温度均值/Cov合同下重构曲面热流。它们尚未与上面的NS推进耦合；不能把一般EOS primitive温度直接当温度体积均值。逐模块用途、真三次制造场误差及可移植测试命令见[二次温度重构说明](../../docs/QUADRATIC_TEMPERATURE_RECONSTRUCTION.md)。完整几何身份、秩不足及非有限通量均明确检查；绝热零热流不求墙温。空间RMS改善不代表正性、完整NS稳定或新增物理解面积。

当前支持的是数值实验。尚需非仿射近壁/网格/外域/特征边界、真实空气物性/湍流适用性和实际区域收敛/实验对照。`physicalSolution=false`和`eligibleForSurfaceCoverage=false`继续保留。未解区不能补零，Tw/储库始终是输入，现有科研五步骤工作台不受本库设置路径影响。

v80新增显式`solidWall.thermalGradient='dirichlet-tangent'`：由已给两切向壁温导数和真实正壁距形成墙面温度梯度，并在扩散模板中采用一致owner导数。常Tw两切向导数为零；函数Tw须明确给切向梯度。默认legacy完整rate/info保持原数值；显式新模式不改质量/动量、内面通量、WLS或墙压。约束模板改为按原顺序预分配，原动态完整实现保留为小几何精确参考。实际20°箭体217272CV、墙热扰动及三步完整NS已测试，仍未收敛或取得物理资格；说明见[真实20°三维研究](../../docs/REAL_20DEG_SHELL_STARTUP.md)。上面的v58计算语句保持是当时的迁移身份，不应当成当前v80源码与v58逐字节相同的声明。

v81显式`boundary.viscousWorkVelocity='arithmetic'`只改共享内部面黏性功速度；默认`legacy-wls`完整数值保持，牵引、墙热、Fourier和边界不变，原无限制WLS迹及负配对功保留。新增保守态导数、两点面线性化、原面法向线装配/批量块Thomas及完整残差方向导数/右预条件Krylov。线近似不替换原三维HLLC/应力/导热方程，伪步不算物理时间或BE。真实六次修正及独立字段/空间五量账已测试，末次仍受HLLC星压限制、尚未稳态；说明见[三维功与Krylov研究](../../docs/REAL_20DEG_LINE_KRYLOV_DRAFT.md)。批量线求解需要当前MATLAB支持`pagemldivide`/`pagemtimes`，本机已实际运行，不依赖外部CFD或并行工具箱。
