# v174 原生定比热首帧输出包：构建前审查

2026-11-06。本文是执行前说明，不是已构建或数值验证的结果。必须先在 PLAN 写入唯一 `v174-output-build-registration` JSON 块，再允许构建脚本执行。本文编写时构建数、链接数和新 CFD 数均为0。

## 原因及唯一修改

small02 原生4更新退出0，但监督器因缺少第0帧退出1，失败身份保持。准确v171 `CFlowOutput.cpp:4138` 的实际重写函数只在 `Iter>0` 或模型为 `IDEAL_AIR_VARIABLE_CP` 时按频率保存。原生 `IDEAL_GAS` 因而只有2/3帧。基类 `COutput.cpp:2176` 也有类似条件，但本调用经过 `CFlowOutput`，修改基类不能修复实际入口。两源及相关继承/驱动头文件已经匹配准确v171源树，不能用当前共享树猜测支持能力。

唯一补丁在 `CFlowOutput::WriteVolumeOutput` 的原条件中增加 `config->GetKind_FluidModel() == IDEAL_GAS`。不改 `COutput.cpp`、流动求解、物性、能量参考、Roe重构、弱壁热公式、17位文本、字段列表、已有频率或强制末帧规则。第0帧仍是第1次完成更新，不是初始场。现有稳态输出频率选项无法绕过门控；不切换非定常模式，也不以事后改成[2,3]使原失败通过。

## 准确来源与对象恢复

父包为 `rah-ideal-air-spatial-wall-v1`，父runtime SHA256：
`e794985c282f73913c679eb72ad1818c9d56e97371813270ec7364865bfa3b3d`。

父源码树 SHA256：
`c600dc43dc777e40a344f253b1adf89e7c51aa95caa98c6a58f57db1c65f0ee8`。

19份父执行源全部来自 `outputs/v171_validation/build01/executed_sources`，逐项绑定父源码树。唯一修改前的 `CFlowOutput.cpp` SHA256 为 `90336701375e196b8cd18e87cd68218a080fec400f94e556aefd094d1770df78`。

共享构建目录实际有298对象。审查时下列4项不等于v171父对象，已经确认它们都能从v171 `produced_objects` 取得准确缓存；不能照搬v173脚本只恢复3项的假设。

| 对象 | 必须匹配的父SHA256 |
|---|---|
| `geometry_CGeometry.cpp.obj` | `f802fce250cf18b0430588f738403aba5a47edd3d35867a4d9bafbb10ba2a8b6` |
| `python_wrapper_structure.cpp.obj` | `01b996a23413e848184c7c23846af982845a9dce23049fdd05fafaeca47ea84e` |
| `solvers_CTurbSSTSolver.cpp.obj` | `1e0a7d570de5e4371da312cd47333d8285b6d9d5e0139cb703f9a70b7bf72105` |
| `variables_CNSVariable.cpp.obj` | `355fdf49bf2dc4b89643cb94ed09be44126f359cde97d3ba4d9c736673230291` |

执行时重新逐项检查298对象，只允许恢复上述准确缓存和重编译唯一 `output_CFlowOutput.cpp.obj`。任何其他非父对象直接拒绝，不进行额外编译。链接前要求297个继承对象全部匹配父SHA，使用父准确链接命令。两薄归档及构建索引保持。成功与失败都恢复原19源、触及的对象/依赖文件和共享exe，并逐项确认全部298对象回到本次开始状态。只有恢复全部通过才发布新包。

## 独立包及元数据

新包具名 `rah-native-ideal-gas-spatial-wall-v1`。只发布新exe、明确列出的许可文件与增量补丁，不递归复制运行包，不更改父runtime或旧small02来源。

- `numericalParentFlavor=rah-ideal-air-spatial-wall-v1`。
- `supportedModel=IDEAL_GAS`，`thermodynamicModel` 为本项目明确的原生定比热NACA模型ID；原父包模型另存 `inheritedPackageThermodynamicModel`，其温区另存 `inheritedPackageTemperatureRangeK`。这只声明本入口支持的模型，变比热实现仍存在于继承包中；不能把二进制称为未修改的官方发行包。
- `firstFrameOutputModelExtension=IDEAL_GAS`，完整输出合同见准备脚本的 `OUTPUT_CONTRACT`，发布为 `nativeIdealGasOutputContract`。
- `samePhysicsAndNumerics=false`（沿用包含输出集合的口径），`numericalUpdateSourceUnchanged=true`、`sameNumericalUpdateSources=true`、`outputScheduleChanged=true`。
- `nativeSharedFrameEqualityVerified=false`，编译通过不预先声称共同帧相同；真实验收另给证据。新模型与旧变比热案例不能因共用父包而称物理相同。
- 未声明核心存在150–360K裁切或适用域保护。该范围仍是外部开发检查，保持原阈值。

## 执行预算与验收边界

最多1对象编译60秒、1次链接30秒、辅助20秒，总110秒、4GiB；失败恢复必须执行，不能通过杀死恢复阶段使预算看似通过。构建脚本不得启动CFD。登记JSON可以由 `registration_specification()` 只读取得，必须与代码逐项相同。

若另行登记新小型验证，仍要求4更新和0/2/3完整帧，并比较与small02共同2/3帧的全部27列和共同历史；数值不变只能在该比较后声明。small02原错误监督回执不改。新包新增第0帧属于输出能力修复，并不授予物理资格。

small02已存2/3帧独立热核验通过，但温度最高分别361.939K和367.697K，超过360K开发上界，且原任务正常完成合同未通过。纯输出补丁不能改善这些温度或物理解。预期共同帧温区失败保持，不扩大温区、不放宽门槛、不自动重跑或启动实际网格候选。是否继续实际大任务须独立作出路线决定。
