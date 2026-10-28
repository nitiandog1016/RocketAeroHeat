# 内部黏性区域实验数值库

本库在MATLAB中实际计算三维定比热理想气体的对流、黏性应力与导热，并按共享全局步长SSP-RK2推进五个守恒量。是后续程序内部回流、横流和区域能量求解的数值基础，未默认替代箭体工程外流，当前不计入箭体热流覆盖。

- `compressibleNavierStokesResidual3D.m`：均匀笛卡尔共享面HLLC、Newton应力与Fourier热流；显式墙速度/壁温、固体掩码和边界积分。
- `resolveNSCharacteristicBoundary3D.m`：显式储库的法向Riemann远场、压力出口及回流。给定出口压力与回流入射数据分别使用；声速转换/不自洽分支拒绝。
- `advanceCompressibleNavierStokes3D.m`：按对流/声学及扩散频率选步，真实RK2与整步拒绝/减半；累计五量改变及两阶段边界通量。目标时间、步数上限和稳态收敛分别报告。

入口示例及全部输入构造见本目录测试。项目根目录运行：

```matlab
maxNumCompThreads(1);
addpath('internalFlow','-begin');
folder=fullfile('outputs','v57_validation','solver');
if ~exist(folder,'dir'),mkdir(folder);end
testCompressibleNavierStokesResidual3D(folder);
testNSCharacteristicBoundary3D(folder);
testNSOpenRegionAdvance(folder);
report=testAdvanceCompressibleNavierStokes3D(folder);
```

检查包括法向不变量/旋转/回流、真实开放输运、两墙导热和Couette黏性发热稳态、闭周期五量守恒、整步回退。使用人工常物性气体的校核明确记录参数，不能称为空气或实验精度验证。Couette的温度二阶网格偏移与停止误差分开列示。

`steadyConverged`只由连续空间残差停止条件取得；`physicalSolution`及`eligibleForSurfaceCoverage`仍为false。给定Tw与储库/出口压力始终标为输入；未收敛墙热流只是当前状态的派生通量。详细结果、来源、适用范围及后续区域接入要求见 `docs/INTERNAL_VISCOUS_REGION.md`。

还缺真实STL/曲面映射、实际区域入口外流及热边界层衔接、网格/域长度/边界收敛、变物性及湍流闭合、实验对照。固体墙采用阶梯掩码，不能把其粗网格热流直接涂回原始三角表面。主轴/鼻尖/来流与SI约定保持原项目要求，未引入外部CFD数据前提。
