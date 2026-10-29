# outputs 后续清理：旧分解因子与Q40准备矩阵

2026-10-29按用户“outputs占空间，可以删就删”的要求完成第二轮核对与清理。此次维护已验收，没有运行新的气流求解或改变物理门槛。

- 删除固定11个MAT缓存，文件长度4,267,156,908字节，磁盘空闲实际增加4,267,143,168字节，约4.27GB。
- 全目录独立核对：只有清单11项消失，其余27,617个原文件的大小、修改时间与文件ID均不变；新增6份`CACHE_PRUNED.json`。
- 49项保留依赖及734项历史MATLAB源码版本逐项SHA核对。Git中没有对应原字节的必要小源码另存raw，没有工程压缩包或大型源码副本。
- 去重后的唯一文件长度由81,722,261,682降为77,455,110,739字节，约77.46GB。此数不声称精确等于磁盘簇占用，也未把硬链接路径长度重复当成可释放空间。

## 具体范围

| 目录 | 已删除缓存 | 保留数据及复查入口 |
| --- | --- | --- |
| v83 `patchhierarchyreal01` | 局部层次分解因子 | 同Q原矩阵、基础因子、右端、数值结果；`studyRocketPatchHierarchy`构造段 |
| v83 `blockfillreal01` | 一级填充分解因子 | 原矩阵、模式、输入、Arnoldi失败证据；`studyRocketBlockFillLevel1`构造段 |
| v83 `overlapFreal01` | 受限重叠分解因子 | 原F图、冻结拓扑、旧强度因子、右端与结果；`studyRocketRestrictedOverlap`构造段 |
| v83 `floworderfactor01` | 流向排序分解因子 | 几何、原矩阵、排序模式和作用结果；`studyRocketFlowOrderingFactor`构造段 |
| v83 `spatialreductionreal01` | 粗层局部分解因子 | 原消元装配、右端、原分解输入和结果；`studyRocketSpatialReduction`构造段 |
| v93 `continuous01/step01/current` | Q40的cache、operator、packed、factor、scale、metric六项 | Q40/Q41原状态、几何、原准备清单、方向、自由前缀和恢复检查点；`prepareRocketCurrentPeakSystem` |

五组v83原任务均terminal、exit0、源码稳定；最新v93原失败任务terminal/exit1，恢复及显示导出任务terminal/exit0，源码稳定。实际Q41路径和SHA已核对，恢复绑定及原约束/单元/界面独立审计通过后才清理Q40的旧准备系统。没有把原失败改记为成功，也没有把Q41的数值接受写成物理稳态或新增合格热流面积。

## 保留与重建边界

Q41原状态SHA为`dc9fea4e454f137fa4e4ead736bed37ffb5407eed2cbc6b5f5012ca9e4ffd1ce`，实际文件位于`outputs/v93_validation/continuous01/step01/krylov01/`，`outputs/v83_validation/krylov41/`保留硬链接别名。Q40前驱、原几何、各轮线性检查点、方向、模型、原失败与历史源码、现有工作台均保留。

下一次实际推进必须从Q41新建系统。需复查历史矩阵时，在独立源码目录恢复计划中的`sourceCommit`及对应raw原字节，按清单SHA确认输入；使用各原研究入口的**构造段**，显式读取保留输入并写入新输出目录。不要覆盖当前项目源码、原结果或硬链接历史文件。此次没有执行重建，不保证新MAT的历史SHA相同；旧清单和方向资格不能直接套在新系统上。涉及被清理因子的历史审核命令须先完成该重建与新的来源绑定。

源码同步到GitHub不会备份被忽略的大型MAT输出。本次只同步维护脚本、用途说明、最新目标记录与小型清理证据；v93数值已测，UI交付收口尚待后续，本次没有将其算法修改一并发布。全面积/全周期目标仍active，正式热流合格面积66.5626327424%、新增0保持。

## 核验资料

- [计划](plan.json)：固定十一项、输入、原任务、来源与源码绑定。
- [执行回执](execution.json)及[逐项日志](operations.jsonl)：实际删除项和磁盘空闲变化。
- [独立复核](postcheck.json)：全目录及关键SHA检查。
- [原目录清单](inventory_before.json.gz)：删除前全部路径的大小、时间及文件ID。
- [规划摘要兼容修复记录](preflight_note.json)：首次清单已完整写入，摘要打印因Python3.8不支持字典并集合并而失败；修正摘要后独立核对清单，任何删除发生前已解决。没有覆盖原清单或重算数值。

`scripts/Invoke-OutputCacheFollowup.ps1`逐项验证绝对路径及全部保留依赖后，仅用原生PowerShell删除固定文件，不递归删除目录。已有执行回执时拒绝重跑；后续清理需要新的独立范围与清单。
