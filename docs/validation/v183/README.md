# v183验证索引

已实现原生实数TVD RK3/局部伪时间入口，唯一6步小例及MATLAB只读调用链通过；未计算整箭，物理资格0。

- [事前计划](PLAN.md)、[结果与限制](../../EILMER_EXPLICIT_ENTRY_V183.md)、[下一项](NEXT_DECISION.md)。
- [证据清单](evidence_manifest.json)、[预算与存储](budget_receipt.json)、[暂存字节与保护核对](delivery_receipt.json)。
- `evidence/small01_prepare/report.json`：原网格/物理初态保持，实际配置差异清单。
- `evidence/small01_run/report.json`：原始误拒仍失败；`process00/report.json`记录唯一原生任务正常退出。
- `evidence/small01_recovery/report.json`：只读恢复，明确输出绑定时刻，原失败不改写。
- `evidence/verify01/report.json`：原面有符号q独立复核、20拒绝与旧NK只读兼容；残差上升如实记录。
- `evidence/matlab_parent.json`及`matlab01/report.json`：实际MATLAB inspect/read/旧失败保持，无新CFD。
- `evidence/source_at_prepare`、`source_at_run`与`source_contract`：执行版本与准确原生源码；最终适配源另列SHA，不能混称同一版本。

保留独有失败及结果；只收集本轮小证据，保护29项并发文件。不给定Tw/有限值/exit0授物理资格，不将伪时间用于轨迹。

本轮183项证据共740076字节，Git暂存原始blob SHA逐项通过；MATLAB执行期间1897项源码/显式输入稳定。没有大场新增或独有删除。
