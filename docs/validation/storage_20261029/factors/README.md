# outputs旧Q17矩阵与分解因子清理

2026-10-29按用户授权直接清理已核实可删除的本地缓存。本次只维护存储，没有重新运行MATLAB求解，也没有交付尚未提交的v94算法修改。

- 实际删除4个文件，文件长度1,288,130,041字节；磁盘空闲实际增加1,288,122,368字节，约1.29GB。
- 原目录27,711个文件中，计划四项消失，其他27,707项在首次独立核验时大小、修改时间和文件ID不变；新增三份`CACHE_PRUNED.json`。
- 随后单独修正过时的`outputs/README.md`：Q40更新为实际Q41，补充最新本地数据与保留规则。原字节、前后身份见`readme_update.json`；最终复核只允许这一文档变更，其余27,706项原文件不变。
- 核对99项保留文件和710项历史MATLAB源码版本；当前Q41、全部v94输出、几何、原始数值结果、失败/中断证据和检查点保留。去重后的唯一文件长度约76.18GB。

## 删除范围及恢复入口

| 原目录 | 删除缓存 | 保留输入与恢复方法 |
| --- | --- | --- |
| `outputs/v83_validation/blockilureal01` | `steady_graph_packed_evidence.mat`、`block_ilu_factor_evidence.mat` | 保留Q17原状态（SHA `71dc35fe56cce987bbe807aa52899a1a17383acd60b667304750f55e99691746`）。使用历史`studyRocketBlockILU0`的构造段，从`result.firstOperator`调用`packRocketSteadyGraph`，再调用`factorRocketBlockILU0`；不执行后续NS试算。 |
| `outputs/v83_validation/coupledpatchreal01` | `coupled_patch_factor_evidence.mat` | 保留原右端`innergraphreal01/inner_graph_real_evidence.mat`和实际作用/求解结果。先从Q17重建packed，再按历史`studyRocketCoupledPatches`的分组与局部因子构造段恢复。 |
| `outputs/v83_validation/fullFstrengthreal02` | `full_F_strength_factor_evidence.mat` | 保留`fullFreal01`的原图、输入、数值结果及`spatialreductionpreflight01`分区。按历史`studyRocketFullFStrength`构造段重建F内分组、因子和覆盖元数据。 |

各原进程均已终态且源码稳定。`blockilureal01`原进程为exit1，其余两组为exit0；这里只回收可构造缓存，原失败与后续诊断证据保留，未把原失败改记成功。未清理`momentumlinereal03`的中断研究或任何线性检查点。

重建应在独立源码目录恢复`plan.json`记录的原Git提交和历史raw原字节，然后按原构造段显式读取保留输入，输出到**新目录**，生成新的来源绑定。原研究入口中包含历史绝对路径或SHA护栏，不能不加适配直接重跑全研究，也不能将新MAT套用旧清单/方向资格或承诺MAT容器的历史SHA相同。本轮没有执行重建。

此前已清理的patch层次、重叠、流向排序等历史研究，若依赖本轮已移除的旧packed/因子，也须先恢复其上游构造。以前的保留清单描述的是当时状态，不代表这些缓存永久存在；以本轮计划和各目录清理标记为准。原结果、原状态与历史来源记录继续保留。

## 核验资料

- `plan.json`：固定四项、输入、源码与目录快照身份。
- `execution.json`、`operations.jsonl`：实际删除与磁盘空闲变化。
- `inventory_before.json.gz`：清理前全部文件的大小、修改时间和文件ID。
- `postcheck.json`：删除后的首次完整独立检查。
- `readme_update.json`、`final_postcheck.json`：单独更新目录说明后的最终独立检查。
- `preflight_note.json`：第一次只读预检查因v94工作文件与旧源码身份不同而拒绝；补充精确Git字节回退并保存历史源码后重新通过，失败时尚未删除任何文件。
- `historical_sources/`：需要按原字节保留的小型历史源码，没有复制大型工程或计算输出。

执行脚本为`scripts/Invoke-FurtherOutputCleanup.ps1`；只用原生PowerShell删除固定四个文件，不递归删除目录。规划与独立核验见`scripts/planFurtherOutputCleanup.py`，有历史回执时拒绝重复执行或覆盖核验记录。

GitHub同步只包含维护源码、说明和小型证据，`outputs/`大型数据仍仅在本地。当前Q41仍未稳态，正式合格热流面积66.5626327424%、本次新增0，总目标保持active。
