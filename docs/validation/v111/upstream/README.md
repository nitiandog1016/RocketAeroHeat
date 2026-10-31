# SU2 8.6.0 官方源码摘录

用于核对几何定义及输出文件语义；这些文件是带原版权头、标注原行号的摘录，不能编译成完整SU2。

- [CPhysicalGeometry.cpp](https://github.com/su2code/SU2/blob/v8.6.0/Common/src/geometry/CPhysicalGeometry.cpp)：`SetControlVolume`，边／面／单元中心及节点控制体累加。
- [CEdge.cpp](https://github.com/su2code/SU2/blob/v8.6.0/Common/src/geometry/dual_grid/CEdge.cpp)：`GetVolume`，三重积绝对值除6。
- [COutput.cpp](https://github.com/su2code/SU2/blob/v8.6.0/SU2_CFD/src/output/COutput.cpp)：稳态基础文件和具名文件均调用`WriteData`，而非只生成具名文件。
- [CConfig.cpp](https://github.com/su2code/SU2/blob/v8.6.0/Common/src/CConfig.cpp)：`GetFilename_Iter`使用局部迭代的至少六位后缀。

Copyright 2012–2026 SU2 Contributors。原代码采用LGPL-2.1-or-later，许可证全文见项目[third_party/su2/COPYING](../../../../third_party/su2/COPYING)。Python几何重构是本项目独立实现，核对方法遵循上述定义；解析单元与分账测试不等于直接比对引擎内部体积数组。
