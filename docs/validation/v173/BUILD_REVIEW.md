# v173 最小运行包构建可行性核对

本文件为只读评审，不是构建回执。2026-11-06 在准确 `rah-ideal-air-spatial-wall-v1` 基线上核对：v165 同守恒增量方向限步、CNSSolver 事件上下文，以及 CSU2IdealAir 属性查询拒绝上下文可以组合为 **三个对象重建、三个准确缓存对象恢复、一次链接**。尚未实施、编译或运行 CFD；是否形成实际候选由覆盖表和独立组合检查决定，不能据本评审认定 v172 异常已修复。

## 准确父包与共享状态

父运行包是 `tools/runtime/su2/8.6.0-rah-ideal-air-spatial-wall-v1`，不能以共享源码目录当前内容推断其能力。只读重新计算的身份如下：

| 项目 | 字节数 / SHA256 |
| --- | --- |
| 父 `bin/SU2_CFD.exe` | 19,501,188 / `3309e77ed0fdfdb465e6d774f2157b7ef8d2cb2252a3d4abbfab7ee49fd7a217` |
| 父 `runtime.json` | 5,630 / `e794985c282f73913c679eb72ad1818c9d56e97371813270ec7364865bfa3b3d` |
| `outputs/v171_validation/build01/source_tree.json` | `c600dc43dc777e40a344f253b1adf89e7c51aa95caa98c6a58f57db1c65f0ee8` |

只读逐项哈希当前 `tools/runtime/su2-thermal-build/build-v117` 的 298 个对象，结果完全等于 v171 `objects_before.json`，亦即该次构建结束后恢复的共享缓存；不能直接将这个缓存当作 spatial-wall 父包。当前缓存与 v171 `objects_after.json` 仅四项不同：CGeometry、python_wrapper_structure、CNSVariable、CTurbSSTSolver。四份 v171 `produced_objects` 全部仍在，逐项 SHA 与父包对象记录相等；此轮哈希与比较用时 0.186357 秒。

准确源来自 v171 的 19 份 `executed_sources`，其他文件依其 `source_tree.json` 身份读取。将覆盖内容仅放在内存中，扫描并逐项核对 7,705 个 C/C++ 源、头及 include 文件，SHA 差异 0，用时 2.990495 秒。共享源码中的 `CSU2IdealAir.hpp` 当前不存在，这正是先前临时覆盖后恢复的状态，不能据此否定已安装包。

父包保留 observed / phase-v2 数值身份、17 位 ASCII 和原有终止观察。空间 Tw 采用 single-zone / single custom wall / MGLEVEL=0 约束，弱热通量公式未改。父 manifest 的 `compiled-combination-tests-pending` 是构建时状态，后续组合检查见 v171 独立证据，不应据此抹掉已完成检查。

## 头文件影响与最小对象集

在上述准确源集合内，`CSU2IdealAir.hpp` 唯一现有直接包含者是 `SU2_CFD/src/solvers/CEulerSolver.cpp:37`；`IdealAir.hpp` 由 CSU2IdealAir 包含。该类的状态查询方法均内联，给六种 `rah_air::state` 调用增加拒绝包装后，必须重建 CEulerSolver，不能只重建旧 v165 的两个对象。

| 必须重建的对象 | 必要原因 |
| --- | --- |
| `SU2_CFD/src/libSU2core.a.p/variables_CNSVariable.cpp.obj` | 以准确 spatial-wall CNSVariable 为底，接 v165 同向限步、原始/旧态/回溯阶段作用域；新增 CSU2IdealAir 类型识别和上下文头依赖。 |
| `SU2_CFD/src/libSU2core.a.p/solvers_CNSSolver.cpp.obj` | 安装按线程作用的现有几何/迭代/线性量上下文提供器。 |
| `SU2_CFD/src/libSU2core.a.p/solvers_CEulerSolver.cpp.obj` | 重新实例化 CSU2IdealAir 内联查询包装，并在 BC_Riemann 等明确选定位置安装阶段/节点作用域。 |

这里的三个对象结论以新上下文头仅被 CSU2IdealAir、以上三个实现文件或它们的局部上下文头包含为前提。不得为方便把新头加入 CFluidModel、CEulerSolver.hpp、CNSVariable.hpp 等公共头；若实现改变了这些公共头或增加其他调用点，必须重新核对实际依赖，不能沿用三对象声明。

当前仅保留 7 个 `.d` 文件，其中 CEulerSolver 和 CNSVariable 的旧依赖记录包含 CSU2IdealAir。它们提供交叉核对，**不是完整当前父包的依赖图**：CNSVariable 记录来自共享数值缓存，父 spatial-wall 的 CNSVariable 实际尚未包含 CSU2IdealAir。本评审的对象范围以准确源扫描、既有对象身份及新包含关系为依据。

不直接复制 v165 的已编对象。v165 的父包是 `rah-ideal-air-admissible-v1`，未含此处查询包装和空间 Tw 组合；应复用其有限算子思想及受测源码，针对 v171 准确源生成新的明确差异。CEuler 重建仅增加已声明上下文与查询包装，不能混入 v159 边界数值改动。CTurbSSTSolver 使用父包准确对象，不移入 v167 的 SST 能量约束。

## 可直接恢复的三个父对象

来自 `outputs/v171_validation/build01/produced_objects`：

| 对象 | 准确 SHA256 |
| --- | --- |
| `Common/src/libSU2Common.a.p/geometry_CGeometry.cpp.obj` | `f802fce250cf18b0430588f738403aba5a47edd3d35867a4d9bafbb10ba2a8b6` |
| `SU2_CFD/src/libSU2core.a.p/python_wrapper_structure.cpp.obj` | `01b996a23413e848184c7c23846af982845a9dce23049fdd05fafaeca47ea84e` |
| `SU2_CFD/src/libSU2core.a.p/solvers_CTurbSSTSolver.cpp.obj` | `1e0a7d570de5e4371da312cd47333d8285b6d9d5e0139cb703f9a70b7bf72105` |

第四份缓存 CNSVariable 为 `355fdf49bf2dc4b89643cb94ed09be44126f359cde97d3ba4d9c736673230291`，身份有效，但此方案已重建该对象，无需再将其当成纯恢复对象重编或复制。

两个静态库均重新读取确认 magic 为 `!<thin>\n`：Common 库 113 成员、SU2core 库 184 成员，六个将触及的对象均在对应库引用中。因此复制三个准确缓存、重建三个对象后，原链接命令会读取实际替换对象；无需重建两个库或全量引擎。余下一对象为既有主程序对象，其异常观察不在此次重建范围。

新组合应有 295 个对象逐字继承 spatial-wall 父包（298 减去三个重建对象）。相对开始时共享缓存，六对象可能被临时改动，其他 292 对象必须逐字不变；这两个计数不得混称。

## 建议事前预算与前置检查

建议沿用已受监督的直接编译路线，固定现有 `compile_commands.json` 和父链接命令，不调用全量 Ninja / Meson：

- 仅三个对象，每个最多 60 秒，累计编译最多 120 秒；串行执行，并从累计剩余额度扣减。
- 仅一次链接，最多 30 秒；父监督最多 210 秒；4 GiB 内存上限，自动重试/追加 0。
- 失败、超时或依赖身份变化即停止并恢复；不在失败后自动扩大对象数或预算。
- 构建只是独立运行包准备。任何实际小任务或新真实候选须另有事前合同、次数与总资源预算，本文件不授权数值重跑。

历史时间仅支持这个预算可行，不是新实测承诺：v165 CNSVariable 2.703 秒、CNSSolver 9.906 秒；v159 单 CEulerSolver 16.000 秒；v171 链接 2.172 秒，包含准确源审计和恢复的父任务 71.047 秒。当前编译器为固定 w64devkit g++，现有 C++17 / OpenMP / x86-64 及浮点选项沿用，不改变优化参数。

正式构建前，最小局部检查应覆盖新 header 的多翻译单元链接、嵌套作用域恢复和 OpenMP 线程隔离，以及合法查询逐值保持、domain_error 的查询模式/两输入/阶段保留、其他异常透传。v165 fixture01 曾出现 `TLS init function for rah_step::provider` multiple definition，fixture02 才通过；该原失败保留。现有成功头采用函数内 `static thread_local`，新上下文实现应先在同一工具链验证实际多翻译单元行为，不能只靠单文件语法检查。

查询包装不能扩大数值保护的声明：它能在实际抛出时保留模式和输入，不自动使重构态或边界态满足可接受域。没有安装作用域的位置必须明确标为阶段未知，不得沿用上一次节点/阶段。边界仍严格拒绝；合法路径保持现有状态求值与一次性 assign 顺序。

## 完整恢复与发布合同

v165 `overlay()` 的 `finally` 仅恢复源码，旧 builder 会在共享对象/可执行文件留下构建产物，**不能直接作为此次共享状态恢复合同**。应复用 v171 已执行的源和二进制双恢复结构，并扩展到三重建对象。

1. 开始前核对父 runtime 全部文件、19 份准确覆盖源、完整父 `source_tree.json`、298 个当前对象、三个准确缓存、编译器/链接器与编译链接命令。固定新计划与实现输入 SHA；记录当前对象、原 exe、三个 `.d` 和每个临时源的存在/不存在状态。
2. 对源使用 v171 的 19 路径覆盖，加新 `ConservativeStepTrace.hpp` 和 `PropertyQueryContext.hpp` 两个路径；按此实现预期共 21 路径。四个既有文本发生明确修改：CSU2IdealAir.hpp、CNSVariable.cpp、CNSSolver.cpp、CEulerSolver.cpp；其他父准确源保持逐字。具体落点最终以实施清单为准，不允许无记录扩张。
3. 二进制恢复集合为六个触及对象、原共享 `SU2_CFD.exe`、三条编译命令 `-MF` 指定的依赖文件，共 10 路径；三份 `.d` 本次只读观察均已存在，仍应在开工时重新捕获。保持两个 thin archive、build.ninja、.ninja_deps、.ninja_log 不变。
4. 在 `try/finally` 内完成三个准确对象复制、三个编译、295 父对象继承/292 共享未触及对象审核、一次链接、执行源及新产物保存。新包输出到独立具名路径，不覆盖 spatial-wall / observed / 主默认；复制新 exe 必须在恢复共享 exe 之前完成。
5. 所有终态恢复源字节及原不存在状态，再恢复上述 10 二进制路径。最后重新枚举且逐 SHA 核对全部 298 对象等于本次开工前，源恢复逐项为真，输入 SHA 保持。恢复失败必须显式失败，不能发布为可用。
6. 为后续最小构建保存本次六个必要对象（三个新编、三个准确继承），不复制全部库。独立 manifest 绑定新 source tree、父 runtime SHA、对象来源、实现/命令/预算/恢复回执。明确 `samePhysicsAndNumerics=false`；只给合法查询路径或弱热公式的窄等价描述，不把限步算子称为原数值完全相同。

此路线不修改物性温区/能量参考、不裁温、不清 k、不补能量、不更换空间 Tw；保持 v165 原收益失败及 v172 未知具体调用者/异常失败。构建成功也不构成收敛、网格、真实热实验或区域资格。该评审期间只读源/对象/运行包，无编译、CFD、共享源/运行包修改、历史清理或 Git 操作；唯一新增产物为本文件。
