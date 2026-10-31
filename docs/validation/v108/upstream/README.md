# SU2 8.6.0 原生方向与重启动依据

2026-10-31通过官方仓库读取固定版本，摘录保留原版权与许可证说明。Copyright 2012–2026 SU2 Contributors；LGPL 2.1或后续版本，完整许可证见 [`third_party/su2/COPYING`](../../../../third_party/su2/COPYING)。摘录用于接口核对，不编译进本项目。

- [`CPhysicalGeometry_orientation.cpp`](CPhysicalGeometry_orientation.cpp)：[官方 CPhysicalGeometry.cpp](https://github.com/su2code/SU2/blob/v8.6.0/Common/src/geometry/CPhysicalGeometry.cpp#L3995) 原3995–4278行，检查棱柱/六面体及边界朝向。SU2边界法向须指向流体域；数学体积审计的外向边界仍单独保留。
- [`CSolver_restart.cpp`](CSolver_restart.cpp)：[官方 CSolver.cpp](https://github.com/su2code/SU2/blob/v8.6.0/SU2_CFD/src/solvers/CSolver.cpp#L2979) 原2979–3112行，ASCII重启动列与节点读取。
- `SOLUTION_FILENAME`是基名，稳态ASCII读取追加`.csv`，见已保留的[CFVMFlowSolverBase.inl](../../v105/upstream/CFVMFlowSolverBase.inl)及[官方来源](https://github.com/su2code/SU2/blob/v8.6.0/SU2_CFD/include/solvers/CFVMFlowSolverBase.inl#L854)。监督器现在锁定实际初场CSV，检查点不得覆写该文件。
- 棱柱换序另对照[官方 CPrism.cpp](https://github.com/su2code/SU2/blob/v8.6.0/Common/src/geometry/primal_grid/CPrism.cpp)。本项目导出时交换0/1及3/4节点，完整数学几何不改变。

这些来源解释格式/方向契约；不能用上游代码存在证明本地物理精度。实际导入日志、独立几何审计和本地源版本在上一级目录。
