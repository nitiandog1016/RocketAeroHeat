# 再次清理 outputs 的实际记录

用户要求解释输出用途并直接清理可删除文件。本次仅做存储维护，没有执行气流求解或修改物理结果。

- 删除：四组旧工况的24个派生矩阵，文件长度8,329,067,703字节；磁盘空闲实增8,329,019,392字节，约8.33GB。
- 保留：最新Q39、旧输入与独有结果、原几何、方向与检查点、公开参考和失败证据、v92小域试验、原有三维界面。
- 剩余：唯一文件长度约79.68GB；路径累计约157.75GB，差额主要来自既有硬链接。唯一文件长度不是磁盘分配空间的精确测量；实收以执行前后空闲差为准。

`plan.json`记录24项原SHA、输入、绑定源码与保护文件；`inventory_before.json.gz`保存完整文件目录，压缩以减少维护资料占用。`execution.json`与`operations.jsonl`记录实际删除。`postcheck.json`是独立Python核对结果：缺失项严格等于清单24项，其余原数据文件的大小、修改时间和文件ID不变，保护文件及源码完整SHA再核对通过。`readme_update.json`单独记录输出目录说明的更新，不将说明文字的修改混入计算数据保持不变的结论。

首次尝试因本机`gzip.compress`不支持`mtime`参数而失败。执行器在预检查阶段拒绝，删除数为0；修为`gzip.GzipFile`并重新生成完整计划后执行。旧计划和失败说明保留，不作为实际执行计划。

## 历史复查

v84三组矩阵使用的历史`cutCellExactEulerBranchJacobian.m`与当前不同，完全匹配的原字节保存在`historical_sources/`，SHA为`1018f50acad9bde439d48e1615f59e1da09d1353ac43e87630bed34bc86ab507`。其余绑定模块与Git提交`267af60ab27fadd56417dcacf8237481eb78f982`逐项核对。取该版本的独立源码目录，按各组`sourceModules.restorePath`放置历史原文件，再使用保留输入，在新目录调用`prepareRocketCurrentPeakSystem`。不得覆盖当前开发源码或复用已删产物的旧清单。必须重新通过原NS回放和来源检查；本次未运行重建，也不保证历史MAT文件SHA逐字节复现。

v91准备系统属于Q38，已经清理；最新Q39本来就需要重新构造自己的导数和方向。原Q39保持未稳态，正式合格热流覆盖66.5626327424%，新增合格面积0。

本次使用`scripts/planHistoricalOutputCleanup.py`生成固定范围清单，`scripts/Invoke-HistoricalOutputCleanup.ps1`逐项检查后用PowerShell删除明确文件，`scripts/auditHistoricalOutputCleanup.py`独立复核；没有递归删除目录。这些是本次固定范围工具，不可通过重复执行恢复旧缓存。
