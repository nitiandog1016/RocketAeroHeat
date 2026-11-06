# v171 准确observed基线的空间壁温运行包

2026-11-06。已在事前预算内构建独立 `rah-ideal-air-spatial-wall-v1`，复用v115空间Tw文件适配。**此报告证明来源、有限构建和共享状态恢复；新包的无表/常表/非均匀表数值组合检查由主任务另行记录。** 构建成功不授予实验或真实箭体物理资格。

## 为什么需要三个编译对象

v170已识别当前observed包没有文件加载。继续核对实际298对象时，另发现共享 `build-v117` 的 `CNSVariable` 和 `CTurbSSTSolver` 仍是v167数值版本，不能只编两个适配对象后直接链接。

- 原observed `CNSVariable.cpp.obj` SHA：`355fdf49bf2dc4b89643cb94ed09be44126f359cde97d3ba4d9c736673230291`。本轮按准确observed源纯恢复编译，最终逐字匹配。
- 原observed `CTurbSSTSolver.cpp.obj` SHA：`1e0a7d570de5e4371da312cd47333d8285b6d9d5e0139cb703f9a70b7bf72105`。现存 `tools/runtime/su2-coolprop-build/cfd-build` 同名对象逐字匹配，直接复用；没有重新编译或借包名推断相同。
- 两个实际适配对象为 `CGeometry.cpp`、`python_wrapper_structure.cpp`。完整对象清单对observed仅这两项改变，其余296项一致；对本轮开始时共享缓存294项未改。

这一修订在构建前写入 [PLAN.md](PLAN.md)，不是失败后扩大预算。没有全量编译，没有把v167状态约束/SST改动混入新层流运行包。

## 适配范围

构建器为 `scripts/buildSU2SpatialAirRuntime.py`。它从 `outputs/v163_validation/build03` 的准确执行源和 `source_tree.json` 恢复10871项来源清单，临时覆盖19项必要源/头文件，再接入逐字绑定的：

- `third_party/su2_thermal/RAHAxialWallTemperature.hpp`：SHA `9c09b738ea9327b3571ccbb716bafac2c5c82f723dcb6205d85cd3389f889920`。
- `third_party/su2_thermal/SetCustomBoundary.cpp.inc`：SHA `e2c0c754a668dba7b9667146d74d71f40c72664978e3d1015c0d2bde15e76dc1`。

加载器要求严格SI表头、递增X、有限正温度和全节点覆盖。只有显式 `RAH_WALL_TEMPERATURE_FILE` 才激活，按X线性插值、不外推；仅支持固定SI、单zone、串行MPI/可OpenMP、MGLEVEL=0、单一自定义等温 `wall`。每个全局壁节点的实际给定温度写入 `rah_wall_temperature_applied.csv`。marker遗漏在原调用点拒绝，不能静默回落常量。

无表时仍走官方标量初始化；空间表只写既有官方逐节点 `CustomBoundaryTemperature` 数组。`CNSSolver` 的弱等温导热公式、Lemmon物性、Riemann边界和17位输出均未修改。给定Tw与方程实际壁流体T仍需分别核验。

新manifest保留observed的 `numericalBaselineFlavor=rah-ideal-air-phase-v2`、17位和原终止观察字段；`numericalParentFlavor=rah-ideal-air-observed-v1`。全包 `samePhysicsAndNumerics=false`，另 `samePhysicsAndNumericsWithoutSpatialTable=true` 仅声明设计保留的无表路径，运行验证状态仍为 `compiled-combination-tests-pending`。它不能作为未运行的数值等价测试回执。

## 实际执行和恢复

唯一 `outputs/v171_validation/build01`：

| 步骤 | 真实耗时 | 终态 |
| --- | ---: | --- |
| CGeometry适配编译 | 9.360秒 | exit0 |
| python_wrapper适配编译 | 3.500秒 | exit0 |
| CNSVariable准确源恢复编译 | 2.593秒 | exit0、对象SHA等于observed |
| 链接 | 2.172秒 | exit0 |
| 外层监督（含准确源核对和恢复） | 71.047秒 | exit0 |

编译合计15.453秒，小于120秒总预算；单项均小于60秒。链接小于30秒，外层小于210秒。监督器4GiB限制保持，过程无失败、超时或自动重试。两个thin archive均核对成员，确保实际链接使用被核验对象。

构建后19项共享源逐字恢复；四对象、原exe以及三个依赖文件共8项共享二进制/依赖状态逐字恢复。再次读取全部298对象，与构建前v167缓存完全一致。新运行包独立路径不覆盖旧observed或其他包；新包必要四个对象单独保存在本轮 `produced_objects`，方便后续准确复用，不复制整个构建目录。

| 产物 | SHA256 |
| --- | --- |
| 新 `bin/SU2_CFD.exe` | `3309e77ed0fdfdb465e6d774f2157b7ef8d2cb2252a3d4abbfab7ee49fd7a217` |
| 新 `runtime.json` | `e794985c282f73913c679eb72ad1818c9d56e97371813270ec7364865bfa3b3d` |
| 新准确 `source_tree.json` | `c600dc43dc777e40a344f253b1adf89e7c51aa95caa98c6a58f57db1c65f0ee8` |

原包文件、适配文件和计划输入SHA在 `input_manifest.json`；实际命令、对象前后清单、准确执行源、缓存对象恢复和最终共享恢复分别保留。大型运行包和对象只本地保留；源码与小证据由主任务统一同步GitHub。

## 后续必须完成

主任务用新具名合同显式开放空间表，保留旧包拒绝；监督器清除继承环境，只传已锁定表。进行已登记的旧包标量、新包无表、新包常表、非均匀表组合检查，核对真实逐节点回执、同阶段热状态和独立热流；不把本构建报告代替那些结果。随后复用v170已过实际网格推进唯一有预算的真实候选。没有在本子任务启动任何CFD、生成新网格、运行MATLAB/UI或提交Git。
