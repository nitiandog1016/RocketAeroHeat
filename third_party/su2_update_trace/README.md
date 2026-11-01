# SU2 原生隐式更新只读记录

固定上游：SU2 8.6.0，commit `7048fd35f56d4bf52e586429d6746a36a10255f8`，LGPL-2.1-or-later。上游版权声明保留；此目录提供明确标记的补丁及记录头文件，大型运行包仅本地保留。

在流动和守恒SST各自更新前后读取实际 Solution、Solution_Old、LinSysSol、共享松弛及采用的密度。所有线程同步后由主线程写出，代码不写数值数组，不修改物理下限、方程、Jacobian或时间推进。环境开关仅由任务合同显式开启；文件存在时拒绝覆盖。当前仅支持固定单区、无MPI、定比热SST的一次更新诊断。

v1曾错误输出RCM重排后的局部编号，保存状态对应审核失败；该运行包和失败结果保留，不能用于真实空间诊断。v2使用 `geometry->nodes->GetGlobalIndex(i)`，读取时检查全局编号完整唯一并排序。源与本段局部步号分开；流动后态和保存primitive的潜在回退差异单列，不能将后更新k混入流动热状态。

构建与检查入口：`prepareSU2UpdateTracePatch.py`、`buildSU2UpdateTraceRuntime.py`、`prepareSU2UpdateTraceCase.py`、`reviewSU2UpdateTrace.py`、`diagnoseSU2NativeUpdate.py`。生成器验证移除观测块可逐字还原固定数值源码；构建器独立安装、审核所有CFD对象和归档成员，最后恢复共享源码。禁止覆盖已有运行记录。

记录没有各项源/输运贡献，不得由总RHS断言某一个源项导致失败。观测兼容与物理资格分开；实际结果以 `docs/validation/v147/` 为准。
