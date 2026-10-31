# 官方实现依据

均为SU2 **v8.6.0**原源码摘录，保留版权头、URL及原行号；Copyright 2012–2026 SU2 Contributors，LGPL-2.1-or-later，完整许可证见[third_party/su2/COPYING](../../../../third_party/su2/COPYING)。不是可独立编译的源码副本。

- [CSysSolve.cpp](https://github.com/su2code/SU2/blob/v8.6.0/Common/src/linear_algebra/CSysSolve.cpp)：BCGSTAB先更新再判断、直接返回从0开始的循环索引，返回残差为相对范数。
- [CEulerSolver.cpp](https://github.com/su2code/SU2/blob/v8.6.0/SU2_CFD/src/solvers/CEulerSolver.cpp)：MUSCL负重构状态的局部退回。
- [CNSSolver.cpp](https://github.com/su2code/SU2/blob/v8.6.0/SU2_CFD/src/solvers/CNSSolver.cpp)：等温壁能量弱条件与速度强条件。
- [CFVMFlowSolverBase.inl](https://github.com/su2code/SU2/blob/v8.6.0/SU2_CFD/include/solvers/CFVMFlowSolverBase.inl)：隐式线性求解、欠松弛更新和历史计数。
- [CFVMFlowSolverBase.hpp](https://github.com/su2code/SU2/blob/v8.6.0/SU2_CFD/include/solvers/CFVMFlowSolverBase.hpp)：对流／黏性时间步取限。
- [CMultiGridIntegration.cpp](https://github.com/su2code/SU2/blob/v8.6.0/SU2_CFD/src/integration/CMultiGridIntegration.cpp)：粗层残差修正与限制／延拓，不将粗层流程与原单层的一次迭代等同运算量。
- [CConfig.cpp](https://github.com/su2code/SU2/blob/v8.6.0/Common/src/CConfig.cpp)：多重网格配置。阻尼的实际代码默认0.5与旧行内描述0.75不同，本次明确写0.5，不依赖描述推测。
