# outputs 目录用途与保留规则

`outputs/` 是本地运行数据区，不是源码目录。它包含历次计算结果、续算检查点、矩阵缓存、验证证据和三维工作台导出。整个目录受 `.gitignore` 排除，GitHub 源码同步没有备份其中的大型计算结果。

| 内容 | 用途 | 保留规则 |
| --- | --- | --- |
| 原几何、已接受状态、边界与气体参数、数值方向/检查点 | 继续当前计算，追溯已有结果 | 保留唯一副本；不能只靠 GitHub 源码恢复已算出的状态 |
| `current_peak_*_evidence.mat` | 缓存某一状态的导数、矩阵、分组因子与度量，避免重复构造 | 保留当前正在使用的系统；旧工况在输入与构造源码核对后可清理 |
| 工作台 HTML/JS/CSV、旧界面导出、重复隔离数据 | 查看三维分布，复核历史界面与交付 | 内容一致时共享磁盘数据，保留原路径；不再默认生成大型重复工程包 |
| 日志、小型 JSON、来源 SHA 与验证回执 | 记录成功、失败、方法和适用范围 | 保留，不能因为清理缓存就改写历史计算结论 |

## 当前续算需要的数据

- 最新接受状态：`outputs/v83_validation/krylov38/wall_velocity_krylov_evidence.mat`。
- 原几何：`outputs/v80_validation/domain02/shell_domain_evidence.mat`。
- 最近方向：`outputs/v90_validation/recovery01/active_direction.mat`。
- 最近准备系统：`outputs/v90_validation/continuous01/step01/current/`。该系统属于当时的输入状态，不得当作 Q38 的导数直接复用；下一步应从 Q38 重建。
- 保留历史输入、实际数值结果、自由/受限检查点、公开参考审计及当前工作台。

## 旧缓存清理后的复算

本轮只纳入构造源码仍与清单 SHA 完全一致的 v86—v89 十组旧准备系统。每组的 `peak_system.json` 保留为原始历史记录；清理后另加 `CACHE_PRUNED.json`，列明输入和被移除的六项派生缓存。旧清单不能直接恢复，需在新目录调用：

```matlab
addpath('internalFlow'); setupCutCellFlow;
prepareRocketCurrentPeakSystem(newFolder, sourceFile, geometryFile, 'domain');
```

`sourceFile` 和 `geometryFile` 来自清理计划对应条目。后续使用新生成的准备清单，不能绕过来源检查，也不能冒称重建 MAT 文件的历史 SHA 会完全相同。重建需要额外计算；本轮没有为了清理重新执行数值求解。

v84 三组旧缓存的一项构造模块已经变化，本轮没有清理。原始数值结果、检查点和失败记录没有因为“看起来旧”就删除。

## 后续约定

每次阶段结束检查缓存是否仍在使用；能由已保留输入和对应源码重建的旧缓存及时回收。保留当前工作所需数据、独有数值结果和小型验证记录，不再默认累积全套旧矩阵、工程 ZIP 或隔离工程副本。具有相同内容的历史大文件可共享一份磁盘数据；资源管理器按路径累计的大小会高于实际占用。

本轮逐项计划和执行回执位于 `docs/validation/storage_20261028/outputs_cleanup_*.json`。实际释放空间以执行回执中的磁盘可用空间变化为准，不能用所有路径长度之和计算。
