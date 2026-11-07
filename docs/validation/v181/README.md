# v181验证索引

[结果与限制](../../EILMER_REAL_NK_V181.md)、[总计划](PLAN.md)、[实数算子合同](REAL_OPERATOR_CONTRACT.md)、[原小例事前计划](SMALL_RUN_PLAN.md)、[下一步](NEXT_DECISION.md)。

| 内容 | 证据 |
|---|---|
| 独立解析尺度 / 6拒绝 | [profile01](evidence/profile01/report.json) |
| 实际准备 / 运行 | [prepare](evidence/small01_prepare/report.json)、[run](evidence/small01_run/report.json)、[实际日志](evidence/small01_run/process00/stdout.txt) |
| **未通过的数值比较** | [small01_check](evidence/small01_check/report.json)，不能以其wallExport子项通过替代总体失败 |
| 首次读取失败源 | [run_small.py](evidence/run_small.py)，错误字段`tke`；准确字段为`tq-tke` |
| 同状态残差与准确源码有限复核 | [difference01](evidence/difference01/report.json)，未认定唯一根因 |
| 失败方案大域/延长拒绝 | [admission01](evidence/admission01/report.json)、[含通用准备入口的最终4拒绝](evidence/admission02/report.json) |
| 具体方法/失败资格导出 | [identity01](evidence/identity01/report.json) |
| 运行时准确源 | [source_at_small_run](evidence/source_at_small_run/report.json) |
| 时间 / 字节与存储 | [budget](budget_receipt.json)、[manifest](evidence_manifest.json) |

v181只有一次实数小域CFD，对照复用v180；无整箭/MATLAB/新主UI/引擎构建或新增物理资格。准入/身份变化在原运行之后，分别有定向检查；不冒称最终源码都在原生小例当时执行。完整大网格仍在既有E盘位置，本轮没有重新打包或上传它们。
