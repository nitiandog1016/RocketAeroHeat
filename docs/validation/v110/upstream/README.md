# 官方输出契约

2026-10-31通过GitHub连接读取固定SU2 8.6.0源码。Copyright 2012–2026 SU2 Contributors；LGPL 2.1或后续版本，完整许可证见[third_party/su2/COPYING](../../../../third_party/su2/COPYING)。文件只用于核对输出契约，不编译进本项目。

- [CFlowCompOutput.cpp](https://github.com/su2code/SU2/blob/v8.6.0/SU2_CFD/src/output/CFlowCompOutput.cpp)：完整原文，节点`RESIDUAL`来自`LinSysRes`，同组包含密度、三动量及能量。
- [CFVMOutput.cpp](https://github.com/su2code/SU2/blob/v8.6.0/SU2_CFD/src/output/CFVMOutput.cpp)：完整原文，原生正交角、面面积比及子体积比字段。
- [CConfig.cpp](https://github.com/su2code/SU2/blob/v8.6.0/Common/src/CConfig.cpp#L3808)：本目录`CConfig_mesh_quality.cpp`保留版权及3808—3815行，`MESH_QUALITY`组自动开启内部开关；不存在独立`WRT_MESH_QUALITY`配置项。
- 迭代顺序另对照已交付的[v105 CFVMFlowSolverBase.inl](../../v105/upstream/CFVMFlowSolverBase.inl)，588—640行：求线性增量后更新解，输出残差和更新后解须区分。

这些代码仅解释接口/语义；实际本地检查见父目录，不能凭成熟引擎名称推定物理资格。
