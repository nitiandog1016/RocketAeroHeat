# 2026-10-30 outputs维护

`outputs/`保存本地计算结果、续算状态、验证证据、派生矩阵和三维界面导出。日常运行不需要携带全部历史输出，但目前独有的大型MAT没有随源码上传GitHub；整体删除会失去已有状态续算和部分历史复核能力。

本轮按用户“占空间，能删就删”的授权执行，没有修改计算算法或重跑MATLAB。删除4个文件，内容长度642,891,500字节；执行窗口磁盘空闲实际增加 **642,883,584字节，约643MB**。

## 删除范围及保留依据

| 删除文件 | 文件长度 | 保留依据 |
| --- | ---: | --- |
| v83 `localrocket02/local_coarse_construction_evidence.mat` | 287,883,394字节 | 9个原变量在同目录`local_coarse_evidence.mat`中逐值完全相等 |
| v83 `localunit02/local_coarse_construction_evidence.mat` | 352,379字节 | 同样逐变量完全相等，保留最终结果及追加验证变量 |
| v91 `continuous01/step01/direction_context.mat` | 177,327,834字节 | 旧QR因子可重建；初始温度行、准备SHA与原`active_direction.mat`完全相同 |
| v92同名旧上下文 | 177,327,893字节 | 同上；实际方向、模型、原失败和线性检查点保留 |

四项均只有一条硬链接路径；本轮没有把重复路径长度当作空间实收，也没有原位改写共享历史MAT。

另核对16份旧`*_raw.mat`，因最终文件没有包含全部原变量或存在差异而未删除。各项差异列于`plan.json`的`rejectedCandidates`；不能仅按文件名中的raw、cache或旧版本号判断可删。

## 空间及核验结果

清理前共有27,982条文件路径，按路径累计150,910,352,121字节；按文件ID去除硬链接重复后为73,503,376,227字节。进一步按Windows压缩文件长度接口计入现有NTFS压缩后约70.20GB，这仍不是精确的磁盘簇统计。

清理后、更新目录说明前，独立遍历确认只有计划4条路径消失，新增4份清理标记，其余27,978个原文件的大小、修改时间和文件ID不变。396项保留文件和411项相关历史MATLAB源码版本SHA通过。当前Q41状态/方向/上下文、原几何、独有结果/失败/检查点和v94—v98全部保留；本轮不交付仍待收口的v98。

去重后的内容长度约72.86GB；计入原有透明压缩后的数据长度约69.55GB。保留数据中v83历史研究占主要部分。数据长度与执行窗口磁盘空闲变化分别记录，不能混用。单独更新`outputs/README.md`后又核对整个目录，只有该说明文件有额外文档变更，详见`final_postcheck.json`。

## 历史复查

- 两份重复中间MAT的全部变量可直接从各自保留最终MAT读取；如需恢复旧文件名，在新目录只导出`plan.json`列出的变量。数值完全保留，但不承诺MAT容器头或文件SHA相同。
- v91/v92旧准备矩阵此前已清理，本次进一步回收QR因子。完整旧审计不能直接运行；先按计划中原状态、几何和准确历史源码在新目录重建准备系统，再由`studyRocketTemperatureLookahead`构造段重建QR，生成新来源绑定。不把旧结果套用新准备清单，不承诺历史MAT的SHA自动复现。
- 原进程的完整来源快照保留。重建所需计算模块及实际驱动依赖逐项核验，未声称无关测试的全部历史版本均已存档。与当前不同的相关源码原字节保存在`historical_sources/`或计划指向的既有档案。

## 记录与执行入口

- `plan.json`、`inventory_before.json.gz`：固定范围、精确比较结果、保留依赖、源码SHA和全目录原身份。
- `execution.json`、`operations.jsonl`：实际删除记录和磁盘空闲变化。
- `postcheck.json`：清理后全目录及保留SHA核验。
- `outputs_README_before.md.raw`、`readme_update.json`、`final_postcheck.json`：目录说明原文、独立文档更新和最终全目录核验。
- `preflight_notes.json`：三次只读规划拒绝及一次PowerShell编码解析失败；这些尝试没有删除文件，修正后才通过全部删除前检查。

规划/审计入口为`scripts/planRetiredOutputCleanup20261030.py`；实际删除入口为`scripts/Invoke-RetiredOutputCleanup20261030.ps1`。已有回执时拒绝重复执行。GitHub同步维护脚本、说明和小型回执，不上传大型输出，不生成ZIP或工程隔离副本。
