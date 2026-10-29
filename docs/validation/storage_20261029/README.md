# outputs 二次盘点与旧构造缓存清理

2026-10-29按用户“占空间、可以删就删”的要求完成存储维护。本轮没有运行气流求解，没有改变算法、当前接受状态或正式覆盖结论。

- 删除：`outputs/v83_validation/currentparent01/`中的operator、packed、factor、image_scale四项构造缓存；`peakweightreal01/`中的cache、operator、packed、factor、scale五项构造缓存。
- 文件长度：3,219,102,007字节；磁盘空闲实际增加3,219,111,936字节，约3.22GB。
- 全目录核对：删除项恰为清单9项；其余27,430个原文件的大小、修改时间、文件ID不变，仅新增两份`CACHE_PRUNED.json`。
- 保留：Q40当前端点、Q39前驱、Q20/Q22历史输入、原几何、原候选/模型与失败证据、当前方向及现有三维界面。16项关键保留文件和776项历史源码版本逐项SHA核对。
- 剩余：按文件ID去重后的唯一文件长度78,036,741,856字节，约78.04GB。路径累计包含硬链接重复；不能把删除一个共享路径的文件长度当作实际释放空间。

## 文件用途与恢复边界

`outputs`是本地运行数据，不是源码。矩阵缓存用于减少重复构造；状态、几何和唯一验证数据用于续算与追溯；HTML/JS/CSV用于查看和导出结果。GitHub排除了大型输出，不能把源码同步当作这些状态已备份，也不能整体删除`outputs`后承诺直接接续原计算。

两组原实际任务的`*_process.json`均已terminal、exit0且源稳定。构造缓存分别属于旧Q20和Q22，当前已经是Q40，实际推进必须生成Q40自己的系统。旧候选和完整模型结果均保留；删除缓存后，依赖原矩阵的完整历史方向审计需要重建。

若需重建，先在独立源码目录恢复`plan.json`的`sourceCommit`，逐项核对`sourceFiles`中的历史SHA，并从`archivedFile`覆盖不在该提交中的原字节，放回`restorePath`。不要覆盖当前开发源码。按原`studyRocketCurrentReservoirParent.m`及`studyRocketCurrentPeakDirection.m`的**构造阶段**在新输出目录重建，使用记录中的原输入：前者保留`reservoirjacreal01`的原缓存，后者使用Q22、原几何和`updateradiusreal01`种子。本次未运行重建；不保证新MAT的历史SHA相同，不能将旧清单或旧方向资格直接套到新准备系统。

三个不参与这些缓存构造的旧模块（显示导出、主线接续、峰值方向求解）明确不纳入构造源码恢复承诺，列表见`excludedNonConstructionModules`。本记录不承诺重放它们的全部历史行为。历史Jacobian的匹配版本已保留；其他未入Git的小型MATLAB源快照保存为原始`.raw`文件，维护资料合计约3.4MB，没有生成大型工程副本或ZIP。

## 核验资料

`plan.json`记录固定九项、输入与源码绑定；`inventory_before.json.gz`为完整原目录；`execution.json`及`operations.jsonl`记录实际删除；`postcheck.json`记录独立全目录和关键SHA复核。清理器只删除绝对路径核对后的九个明确文件，不递归删除目录。

第一次只读规划发现原运行锁包含未入Git的`outputs`源快照，已在任何删除前拒绝；修复为保留这些文件的原字节后完成规划。`preflight_note.json`记录该情况，没有删除计算数据。

本轮维护验收通过。此前v93试验代码仍未测试、未交付，本次不会将其随维护同步。总目标仍active；Q40仍未稳态，正式合格热流覆盖仍为66.5626327424%，新增0。
