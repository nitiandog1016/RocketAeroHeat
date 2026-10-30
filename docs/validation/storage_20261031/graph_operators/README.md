# 八份旧图预条件研究的存储回收

2026-10-31按用户“outputs占空间，能删就删”的要求执行。本轮仅维护存储，没有运行MATLAB、修改物理算法或重算流场。

| 项目 | 已执行结果 |
| --- | ---: |
| 原MAT容器 | 8份，8,710,739,258字节 |
| 保留科学结果MAT | 8份，5,228,726,855字节 |
| 净减少的数据长度 | **3,482,012,403字节，约3.48GB** |
| 迁存前至独立核验结束磁盘空闲增加 | 3,477,745,664字节，约3.48GB |
| 清理后去重并计入已有NTFS压缩的数据长度 | 60,080,118,223字节，约60.08GB |

磁盘空闲也受其他进程影响；数据长度不等同于磁盘簇精确占用。execution中的删除窗口包含已提前另存的结果，不能当作净回收量。

## 已清理内容与完整保留内容

范围固定为`outputs/v83_validation/krylov09`—`krylov15`和`krylov17`。移除57项`result.firstOperator`构造字段：`diagonal/lower/upper`、`minusCoupling/plusCoupling`、`lineInverse/lineUpperFactor`及k17的`secondLower/secondUpper/wlsNormalFactor`。每项普通缓存须是5×5×N的double块；WLS结构须且只能包含四项对应分解矩阵，夹带其他结果会拒绝。

初始与最终Q、完整残差和info、墙面压力/热流、气体/边界、迭代/线搜索、数值方向、JVP和失败记录等756个字段的原MAT矩阵数据元素逐项核对类型、形状、长度和SHA。图着色/分组/面号、构造选项、误差和WLS报告亦保留。其他顶层变量直接复制原压缩记录，没有降低精度、重新编码科学数组或把未知值填零。

七项原进程实际exit0/源稳定；krylov15对应recycledtrial01实际exit1，原因是数值已保存后的旧UI导出失败，不能称宿主成功。其数值结果和失败记录全部保留；recycledexport01另行重导出实际exit0/源稳定。本轮没有重算数值。首稿按编号寻找不存在的krylov09_process.json，在迁存和删除之前拒绝；现显式绑定真实进程名与MAT中的引擎ID，首稿盘点与失败说明保留在preflight。

## 保留位置与继续使用

保留MAT位于`outputs/retained_graph_results_20261031/v83_validation/krylovXX/wall_velocity_krylov_evidence.mat`。八个原目录新增`OPERATOR_CACHE_PRUNED_20261031.json`，提供原身份、准确新路径和缺少字段。这些是缺少旧预条件矩阵的科学结果文件，不能当成完整缓存或直接套用原MAT容器SHA。

当前Q41状态/方向/上下文、v101最终结果、v102全部MAT、原几何以及全目录其他独有结果、失败和检查点保持。用户日常运行程序不需要携带全部历史输出；大型结果未上传GitHub，源码同步不是它们的备份。

历史预条件器复核应从保留结果的initialQ、gas、boundary、residualScale、pseudoTimeStepSeconds、fixedSymmetricSweeps及原几何出发，按plan/source_archive绑定的准确历史源码，在新目录重建`buildRocketShellGraphPreconditioner`；k17还需要相应选中方向缓存与`buildRocketWlsNormalPreconditioner`。上游若已按此前回执清理，也须在新目录重建全套构造与新来源，不能将旧审计入口当作原缓存仍在。本轮未执行重建，不保证历史MAT的SHA复现。

## 核验资料与同步范围

- `plan.json`、`inventory_before.json.gz`：固定八份容器、57项缓存、756个保留字段及原路径身份。
- `tests.json`：6项测试通过；覆盖特殊浮点值、负零、整数类型、中文、原字节保留、非法块形状及损坏/夹带字段拒绝。
- `execution.json`、`operations.jsonl`：固定白名单、绝对路径及全部SHA复核后逐文件删除，没有递归删除目录。
- `postcheck.json`：另一个进程重新读取全部迁存MAT，核对28,254个其余原文件身份、93项保留依赖、100项准确历史源码。
- `source_archive.json`：准确历史源码绑定；新增26份原字节小型归档，其余复用现有归档或固定Git字节。
- `readme_update.json`、`outputs_README_before.md.raw`、`final_postcheck.json`：核验后仅更新outputs用途说明，保留其原文并复核。

源码、测试、小型回执和目标记录按约定同步公开GitHub。v102未收口内容和OpenCode并发修改不纳入本次提交；不上传大型outputs，不生成工程ZIP。正式热流覆盖仍66.5626327424%、本轮新增求解面积0；Q41仍未稳态，完整气动热目标继续有效。

入口为`scripts/compactRetiredGraphOperators.py`的plan/archive-sources/audit/final-audit，以及固定执行器`scripts/Invoke-RetiredGraphOperatorCleanup.ps1`。已有回执不可覆盖，不应重复执行本轮删除。
