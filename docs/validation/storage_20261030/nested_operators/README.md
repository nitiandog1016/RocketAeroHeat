# 旧结果中的预条件矩阵清理

2026-10-30按用户“outputs占空间，核实后直接清理”的要求执行。本轮只维护存储，没有运行MATLAB、改变计算方法或新增物理解面积，未将v102草稿和OpenCode并发修改混入交付。

## 已执行的范围与空间

| 项目 | 实际结果 |
| --- | ---: |
| 退役原MAT容器 | 10份，9,338,390,756字节 |
| 保留的数值结果MAT | 10份，6,903,568,796字节 |
| MAT容器净减少 | **2,434,821,960字节，约2.43GB** |
| 迁存前至独立审计结束的磁盘可用空间增加 | 2,431,717,376字节，约2.43GB |
| 清理后outputs去重并计入现有透明压缩的数据长度 | 63,546,302,758字节，约63.55GB |

磁盘空闲变化还受其他进程影响。`execution.json`中的删除窗口包含已提前另存的6.90GB结果，不能把原容器合计9.34GB当作净回收。数据长度也不等同于磁盘簇统计，资源管理器按路径累计还会重复计算硬链接。

仅处理`outputs/v83_validation/krylov01`—`krylov08`、`krylov18`、`krylov19`。十项实际计算均已终态、exit0、源码稳定，且待删除原文件只有一条硬链接。冻结目录、当前Q41、当前方向/上下文、原几何和v101内核/最终工作台不在删除范围。

## 删了什么，保留了什么

删除的是`result.firstOperator`中44项预条件矩阵/分解字段：`diagonal`、`lower`、`upper`，以及对应旧研究中的`minusCoupling`、`plusCoupling`、`secondLower`、`secondUpper`和`wlsNormalFactor`。后者核实仅含`inverse`、`upper1`、`upper2`、`effectiveLower1`四项矩阵，不夹带其他结果。移除的原始未压缩MAT数据元素共2,867,993,856字节，该值不是磁盘回收量。

初始/最终Q、残差及完整info、面元压力和热流、边界与气体参数、迭代/线搜索记录、方向、JVP对照和失败证据全部保留。预条件器的算法、伪时间步、尺度、面/层身份、颜色分组、误差和WLS构造报告也保留。没有按文件名含“cache”就删除数值证据。

保留位置为`outputs/retained_operator_results_20261030/v83_validation/krylovXX/wall_velocity_krylov_evidence.mat`。每个原研究目录新增`OPERATOR_CACHE_PRUNED_20261030.json`，记录原容器身份、精确保留位置和移除字段。

迁存采用流式MAT5处理，不把GB级结构整体加载进内存。仅更新两级结构长度/字段表并移除白名单元素；**662个其他字段的原MAT矩阵数据元素按类型、形状、长度和全部字节SHA核对一致**。其他顶层变量的原压缩记录直接复制，没有重新编码科学数组或降低精度。新容器需要重新压缩result，因此不能套用原容器SHA，也不承诺重建原MAT容器SHA。

## 当前使用与历史复查

当前Q41仍在原路径，v101的完整科学MAT和最终view05保持。用户日常运行不必附带全部旧输出；本轮保留的独有历史数值尚未上传GitHub，删除整个outputs仍会丢失这些结果。

历史研究若只查看数值结果，按各目录清理标记加载迁存MAT。它是明确缺少预条件矩阵的结果文件，不能当成原完整缓存、直接运行旧预条件器审计或复用旧文件SHA绑定。

若需旧矩阵，应从迁存结果的`initialQ`、气体/边界、原几何及保留的尺度/伪时间步/模式出发，在新目录按准确历史构造源码重建法向线、横向图和所需WLS修正/分解，再建立新的来源绑定。相关入口包括`buildRocketShellLineOperator`、`buildRocketShellGraphPreconditioner`、`buildCutCellSelectedDirectionCache`和`buildRocketWlsNormalPreconditioner`。不能把后续接受状态的方向或旧已清理粗基当成原构造仍可直接加载；本轮没有执行缓存重建。

164项实际历史源码版本已核对，其中原来仅在outputs中保留的20份源码原字节另存至本目录`historical_sources/`，共201,453字节，随小型回执同步GitHub。`source_archive.json`记录准确原位置、恢复路径与SHA；原冻结副本未删。

## 核验与入口

- `tests.json`：四项实际单元检查通过，覆盖NaN/正负无穷/负零/中文/整数类型和原字节保留、损坏压缩流拒绝、缺少白名单拒绝、WLS夹带额外结果拒绝；小型迁存MAT经SciPy实际读取。
- `plan.json`、`inventory_before.json.gz`：固定十份原容器、44项缓存、662个保留字段及全部原路径身份。
- `execution.json`、`operations.jsonl`：原文件SHA及保留副本核实后，以固定路径逐文件删除。
- `postcheck.json`：另一个实际进程重新读取十份保留MAT，通过字段与原记录检查；其余28,166个原文件身份不变、79项保留依赖和164项历史源码SHA通过。
- `outputs_README_before.md.raw`、`readme_update.json`、`final_postcheck.json`：上述审计后仅更新outputs用途说明，并记录最终核验。

入口为`scripts/mat5NestedCache.py`、`scripts/compactNestedOutputOperators.py`和`scripts/Invoke-NestedOperatorCleanup.ps1`。规划/迁存阶段不删除原文件；执行器先复核固定白名单与全部SHA，再执行单文件删除，无递归目录删除。维护回执不能当作新增物理或实验验证。

`audit`对应更新outputs用途说明之前的完整数值审计；当前文档更新后的状态用`python scripts/compactNestedOutputOperators.py final-audit`复核，允许的原文件变化仅为已保存原文与SHA的outputs/README.md，不重复运行MATLAB。
