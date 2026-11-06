# v174 原生定比热路径：准确运行包与接口评审

2026-11-06，按本轮 `PLAN.md` 第一阶段执行。**源码和对象证据支持复用 v171 空间 Tw 高精度包，下一步仅进入另行登记的最小定比热接口任务；当前没有完成任何新的气体 CFD 验证。** 不使用 v173 限步包，不改变旧 Lemmon/OTH 合同，不把本次模型改变解释成旧异常的原因。

## 唯一有限审核

脚本：`scripts/reviewNACAConstantCpRuntime.py`。输出：`outputs/v174_validation/runtime_review01`，保存准确审核源、计划快照、22 份准确执行源、具名摘录、对象绑定与报告。

监督上限 40 秒，实际 **0.234 秒、exit 0**；脚本内部 0.109 秒。22 源、14 对象、26 项具名断言通过，46 项输入 SHA 稳定。没有 CFD、构建、大场读取或共享源改动。

- 实际包：`tools/runtime/su2/8.6.0-rah-ideal-air-spatial-wall-v1`。
- 可执行文件 SHA256：`3309e77ed0fdfdb465e6d774f2157b7ef8d2cb2252a3d4abbfab7ee49fd7a217`。
- v171 准确源清单 SHA256：`c600dc43dc777e40a344f253b1adf89e7c51aa95caa98c6a58f57db1c65f0ee8`。
- 构建记录含 298 对象；本次逐字绑定其中 14 个相关对象。原生 `CIdealGas` 对象 SHA 为 `b3d0ccdee27584ccf30a66ab14136c10220ef5428d472ee4bab945c5c0eb0ea4`，SIMD 数值对象为 `1ed70baf0cb2baec57802cf1491a1f4994b0bc70b91210048f7e944685448ee7`。

包内仍存在项目变比热/空间温度/输出适配，不能称“未修改官方发行版”。`runtime.json` 的原 Lemmon 元数据保持原文；**包身份和实际 case 物理模型必须分开**，由新的 `ideal_gas_contract.json` 说明此次实际使用 `IDEAL_GAS`。

## 已确认的实际路径

以下行号指准确 v171 源。每份源的完整 SHA、选中路径和本轮准确副本见报告，不依据未经核对的共享树。

| 路径 | 源码依据 | 对本轮的含义 |
|---|---|---|
| 原生气体与 gamma/R | `CIdealGas.cpp:31–101`；`CEulerSolver.cpp:888–891,1145–1147`；`CConfig.cpp:1326–1329` | `IDEAL_GAS` 创建 `CIdealGas(gamma,R)`；应显式给定 gamma=1.4、R=8.31451/0.02896546。不要用会覆盖 R 的 `STANDARD_AIR` 代替。 |
| 守恒能量 | `CIdealGas.cpp:44–67`；`CNSVariable.cpp:140–206` | `cv=R/(gamma−1)`，`e=cv*T`，`rhoE=rho*(cv*T+|u|²/2)`，层流 k=0。没有 Lemmon 的 OTH 常量参考。原 OTH 初场/冷重启的 `rhoE` 不可直接沿用。 |
| 单元原生恢复 | `CNSVariable.cpp:153–185` | 保留原生密度/压力/声速/温度检查与旧态恢复，没有 v173 同方向回溯或查询上下文补丁。原生气体模型没有项目 60–2000 K 域拒绝，不能据此称无条件稳定。 |
| 主内部对流/黏性 | `CEulerSolver.cpp:1838–1849`；`CNumericsSIMD.cpp:38–55,97–119` | 理想气体 Roe 走 `EdgeFluxResidual` 和原生 SIMD `CRoeScheme<CCompressibleViscousFlux>`，没有进入项目一般 EOS 重构异常捕获路径。 |
| 原生 MUSCL | `numerics_simd/flow/convection/common.hpp:74–130` | 原生重构 T/速度/p，按 `p/(RT)` 求 rho，按定 cp 与动能增量更新 H；原生非物理边检查/一阶回退保留。不是项目 variable-cp 回退，也不是全残差线搜索。 |
| 标量黏性/边界 | `CDriver.cpp:1667–1682,1804–1828`；`flow_diffusion.cpp:393–446` | 创建原生理想气体 Roe / `CAvgGrad_Flow`，按原 primitive 索引读取 μ、κ、cp；不走一般 EOS 黏性路径。 |
| 导热率 | `CFluidModel.cpp:102–108`；`CConstantPrandtl.hpp:49–55`；`CConfig.cpp:5772–5781` | 原生层流 `κ=μ cp/Pr`。NS 时关闭湍流导热项，`μ_t=0`。项目 variable-cp 的 `dcp/dT` 导热修正不激活。 |
| 弱等温壁 | `CNSSolver.cpp:690–750` | 保留原生无滑移和 `κ*[-(T_neighbor−Tw)/distance]*area` 能量残差。variable-cp 的 secondary 导数分支仅受相应模型枚举激活；IDEAL_GAS 使用原生 gamma/R 近似隐式对角。不能强制覆盖求解壁流体温度。 |
| Riemann 边界 | `CEulerSolver.cpp:5331–5348,5399–5427,5469–5488,5527–5558` | 保留 `outer` 静 PT 和 `outlet` 静压标记。外态按原生气体求取；特征混合仍按原式形成。共享 `thermal_energy` 辅助函数只做 E−速度平方/2−k，不引入 OTH。 |
| 空间 Tw | `CGeometry.cpp:2752–2817` 与准确 `RAHAxialWallTemperature.hpp` | 表加载不按流体模型门控；仍要求单区、串行 MPI、固定 SI 网格、MGLEVEL=0、唯一 custom `wall`。逐节点输入回执17位输出，不外推。 |
| 保存 | `CSU2FileWriter.cpp:44–68`；`CFlowCompOutput.cpp:220–276,309–327` | Restart ASCII17位精度对所有模型有效；`Flow_Primitive_K` 只在 variable-cp 模型注册，因此本轮不能要求该列。表面 VTK 仍须沿原文件精度审核，不误称所有格式都17位。 |

## 三维数组和气体输出合同

原生三维 primitive 共 12 项，**零基索引**为：

| 量 | T | u/v/w | p | rho | 总焓H | 声速a | μ | μ_t | κ | cp |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 索引 | 0 | 1/2/3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |

`EulerNSecVar` 对 IDEAL_GAS 返回 **0**，不能读取 variable-cp 的 secondary[2]/[3] 等列。当前 Roe 梯度分配为 nDim+2=5（T、三速度、p），重构容器另外保存由原生关系求取的 rho/H，不能将两种数量混淆。

限定 `VOLUME_OUTPUT=(COORDINATES,SOLUTION,PRIMITIVE,RESIDUAL)`、三维层流、无可选附加输出时，预计 Restart 为27列：PointID、xyz、五守恒量、p/T/Mach/Cp/三速度/μ/三摩擦系数/q/y+、五残差。这里压力系数 Cp 与定压比热 cp 是不同量。完整预测列名已放在审核报告，必须由唯一小例实际列头确认，不能只据源码直接授予输出资格。

原生输出没有单独保存定压比热 cp 和导热率 κ；若后处理需要，应以具名定比热方法从明确 gamma/R、实际 μ 和 Pr 推导并标为派生量，不冒称引擎输出字段。层流不需要伪造 `Flow_Primitive_K`，也不允许沿用 SST/OTH 阶段标签。

## 建议的唯一最小接口任务

保持 v171 已验证小网格、原单 wall 非均匀 Tw、`outer/outlet/sym_y/sym_z` 标记及原生弱壁公式。**保留 Riemann，不临时换 farfield。** 静 PT 外态的方向输入实际是 Mach 向量，应为 `[-1.97,0,0]`，不是单位方向 `[-1,0,0]`。外场 PT 和速度必须按新模型相容条件重算；这是新物理模型初场，不是 v171/v173 结果的连续接续。

独立接口合同至少固定：准确 v171 二进制/源/对象身份；IDEAL_GAS 与显式 gamma/R；Sutherland 和 Pr；全新 cv*T 守恒初场；准确 Tw 表与实际施加回执；上述 Riemann 标记与 Mach 向量；模型专用保存列/阶段；启动前 resolver、容量监督及 checkpoint 分派；正常与异常任务保存帧身份。不得修改原运行包元数据来冒充另一包。

每帧须独立核对有限正状态、`p=rho*R*T`、`rhoE`、速度、μ、κ、原列热流和功率，原 −8 与严格热梯度 2e−10 门槛保持。给定 Tw、求解壁流体 T 和固体材料 Tw 分开；原生层流 y+=0 不作为真实近壁资格。

本轮计划的 [150,360] K 是外部开发适用域检查，既不是原生模型硬温区，也不是已证实的全场包络。超域必须报告，不截温或自动扩域。原生压力/温度非正检查不能替代独立有限性审核；本次静态路径通过也不证明收敛、实验精度、全域边界可接受性或真实20°整箭能力。

当前结论允许主任务另行登记最小组合任务。尚未运行该任务，更未授权或开展实际大网格候选；物性误差、原文物性表/恢复关系、后6段中心、外缘和物理验收仍各自保持边界。
