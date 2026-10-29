# outputs 旧 QR 构造缓存清理

2026-10-29按用户“占空间，能删就删”的授权，清理已结束研究中的可重建构造缓存。本轮没有重算气动热，没有交付仍在收口的v95，也没有删除当前续算状态。

- 删除11个独立数据对象的25条硬链接路径。按路径累计为4,425,819,194字节，独立对象文件长度为1,947,178,983字节；磁盘空闲实际增加 **1,946,992,640字节，约1.95GB**。
- 删除后、单独更新目录说明前，全目录27,816个保留原文件的大小、修改时间、文件ID不变。目录说明原字节另存，最终仅允许这一文件有文档变更。
- 独立复核354项保留文件和836项历史MATLAB源码版本SHA；九份上下文的温度约束行，与保留的`active_direction.mat`逐值完全相等。
- 最终去除硬链接重复后的唯一文件长度约 **74.23GB**。这不是磁盘簇精确占用；不能用所有硬链接路径长度之和判断实际占用。

## 删除范围

| 原文件 | 数量 | 保留内容 |
| --- | ---: | --- |
| v85 `qr_real01/augmented_qr_factor.mat` | 1对象、2路径 | 准备清单、输入、QR方向、完整模型、实际原NS结果与验证记录 |
| v86 `comparison01/temperature_projection_system.mat` | 1对象、2路径 | 准备清单、自由/受限方向、约束行、分支结果及验证记录 |
| v87 `continuous01/step01`—`step03`的`direction_context.mat` | 3对象、6路径 | 原方向、同值约束行、实际模型、线性检查点与接受状态 |
| v88 `continuous01/step01`、`step02`和`recovery01`的同名上下文 | 3对象、6路径 | 原失败与恢复研究、方向、检查点及实际状态；失败的step03上下文未删除 |
| v89 `continuous01/step01`—`step03`的同名上下文 | 3对象、9路径 | 方向、同值约束行、模型、检查点与接受状态 |

硬链接路径还包括旧隔离验证目录，全部别名均位于本工程`outputs/`内，删除前按文件ID、链接数与SHA核对。没有原位改写历史共享MAT文件。

## 为什么不删除整个 outputs

源码和日常运行不需要附带全部历史输出。不过，当前Q41状态、原几何、最近方向及构造上下文是续算输入；旧研究的独有数值结果、失败与线性检查点用于复核历史结论。这些大型文件没有上传GitHub。整体删除将失去已有状态续算和部分历史复核能力。

本轮保留当前Q41上下文、全部独有结果、失败、线性检查点、v94及v95本地验证数据和工作台；未启动任何MATLAB计算。正式合格热流面积仍为66.5626327424%，本次新增0，总目标保持active。

## 旧研究如何恢复

各旧准备系统的六项大矩阵此前已回收，本次再回收下游QR构造。原`peak_system.json`保留为历史记录，各被清理目录新增`QR_CACHE_PRUNED.json`；旧完整审计和续接入口不能直接使用。

1. 按`plan.json`的源码SHA恢复对应版本，使用所记录的原状态、几何及参考方向；历史差异源码原字节在`historical_sources/`及此前清理档案中。不要覆盖当前开发源码。
2. 在新目录通过对应历史版本的`prepareRocketCurrentPeakSystem`重建准备系统，得到新清单及来源绑定。
3. 按历史`studyRocketAugmentedPatchQR`、`studyRocketTemperatureProjectedQR`、`continueRocketTemperatureActive`或相应恢复入口的构造段，重建QR因子、温度约束行及投影。原参考方向、活跃循环和检查点继续保留。

本轮没有执行重建。重建需要计算，且不能承诺MAT容器的历史SHA相同；不能把旧方向或旧资格直接套到新清单。保留的历史方向审计若读取已删除的上下文，需先按上述来源重建整个匹配研究及新绑定。

## 核验资料

- `plan.json`、`inventory_before.json.gz`：固定范围、全部硬链接别名、保留输入和源码身份，以及清理前全目录快照。
- `execution.json`、`operations.jsonl`：实际逐文件删除和磁盘空闲变化。
- `postcheck.json`：删除后的独立全目录与保留SHA核验。
- `outputs_README_before.md.raw`、`readme_update.json`、`final_postcheck.json`：单独更新目录说明的原字节、前后身份和最终全目录核验。
- `preflight_note.json`：初次预检查因一项隔离验证源码版本不同而拒绝；没有删除文件，找到SHA完全匹配的历史原文后才重新通过。

规划与复核脚本：`scripts/planRetiredQRCacheCleanup.py`；实际删除脚本：`scripts/Invoke-RetiredQRCacheCleanup.ps1`。使用原生PowerShell按固定清单逐文件删除，不递归删除目录；已有回执时拒绝重复执行。GitHub只同步这些维护源码、说明和小型证据，大型`outputs/`仍不纳入版本库。
