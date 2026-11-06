# v173 下一路线决定：关闭当前节点后处理限步支线

2026-11-06。决定是关闭当前 `spatial-step` 候选及其节点后处理限步支线，不登记新的同源试算，不直接开启 `NEWTON_KRYLOV`，也不进入自研全残差线搜索内核。该决定针对目前保护范围和开发成本，**不是对 exit 134 原因的诊断**；最终终止原因仍按 `FAILURE_REVIEW.md` 保持未知。

## 有限准确源证据

局部读取执行约 3 秒，涉及文件先与 v173 准确 `outputs/v173_validation/build01/source_tree.json` 的 SHA 核对相同；没有大规模源码搜索、构建、CFD或大场读取。本页只整理已经取得的检查结果，没有为文档再次核验源文件。

| 准确源相对路径 | 已检查位置 | 结论 |
|---|---|---|
| `SU2_CFD/src/solvers/CEulerSolver.cpp` | 2567–2602 `ComputeUnderRelaxationFactor` | 根据密度和静内能相对变化计算原生局部松弛；不是完整 60–2000 K 物性域与边界残差接受测试。 |
| `SU2_CFD/include/solvers/CFVMFlowSolverBase.inl` | 588–613 `CompleteImplicitIteration` | 直接把原生 `UnderRelaxation × LinSysSol` 加到全部流动守恒分量，然后通信；没有试探全残差后统一接受/拒绝循环。 |
| `SU2_CFD/src/solvers/CSolver.cpp` | 2065–2086 | CFL 反馈读取原生流动/湍流松弛因子；不读取本轮 `CNSVariable` 后处理得到的 `theta`。这说明覆盖差别，不证明它导致终止。 |
| `Common/src/CConfig.cpp` | 1905–1906 | `NEWTON_KRYLOV` 确为真实配置项，默认 false；不能虚构配置，但也不能把真实开关解释成它未实现的功能。 |
| `SU2_CFD/src/integration/CNewtonIntegration.cpp` | 123–158、193–208、298–317 | 具备扰动和完整残差计算，用于矩阵自由线性求解；309 行明确保留 `TODO: Clever back-tracking and CFL adaptation based on residual reduction.`。313 行仍调用 `CompleteImplicitIteration`，317 行再预处理。没有可直接启用的全域试探接受/拒绝机制。 |
| `SU2_CFD/include/integration/CNewtonIntegration.hpp` | 197–210 | 声明扰动、残差计算和有限差分步长接口；这些接口本身不是可接受域回滚合同。 |

准确 SHA256：

```text
SU2_CFD/src/solvers/CEulerSolver.cpp
d85cd9efa8372b7a3bc1737c73245538affc4445339421355ee6d2d942214ee8
SU2_CFD/include/solvers/CFVMFlowSolverBase.inl
ef228099073057584c196405b674ff8dbf246a3fa325e636e36eff2a423c4e68
SU2_CFD/src/solvers/CSolver.cpp
b9393d14a1c651d6cceeac7177507b13b2f23e3bbb576171aa65b6840514b400
Common/src/CConfig.cpp
bb137ca80063b9c9d775c8ff576b155ff3a0ad48bb8998572fd34afb3b54ec47
SU2_CFD/src/integration/CNewtonIntegration.cpp
77aa24c092bc48184eac3f3bd10a1e4dc987e1876bb107b7c9cf1a21b87f5d58
SU2_CFD/include/integration/CNewtonIntegration.hpp
63307e7077e6353a0d32f5a256a2c46a71990ce18144428d06b8e800e3c74c8c
```

CEulerSolver 和 CConfig 的准确文本取自 `outputs/v173_validation/build01/executed_sources`。其余文件在核对时共享路径 `tools/runtime/su2-thermal-build/source` 与准确源清单一致；后续共享树变化时应按 SHA 寻找准确源，不能直接沿用当前共享内容。

## 为什么停止，而不是再加一个开关

当前补丁仅在节点物性查询之后回溯守恒态，MUSCL 原配对回退保留，静压组合和特征混合边界继续严格拒绝。把真正的全域试探残差接受/拒绝补进引擎，至少需要定义守恒态、旧态、梯度/限制器、边界缓存和残差的完整恢复，以及 OpenMP 各线程失败的统一汇总；还需要把接受步长和推进历史的语义闭合。这是新的求解器级开发，不能作为“小适配”自然延长当前候选。

原生 Newton–Krylov 也不能靠一个开关获得这些缺失能力。本次没有足够证据支持唯一确定的网格或出口边界改动，因此不因最后几条恢复事件的位置而改出口，不扩大网格碰运气，不扫描 CFL、线性参数或初猜。

停止范围是这一条实现与试算支线，不是项目总目标。后续主任务仍应优先回到具名适用条件下的真实20°全面积、分离再附、底面/突起和科研工作台；任何替代成熟引擎/受支持物理路线须另有明确依据、来源和有限验收计划，不把本轮机制检查称作新增求解能力或物理资格。
