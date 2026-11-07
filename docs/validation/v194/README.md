# v194证据索引

本轮接受完整原壁Euler网格并完成真实M5/20°原生准备，没有流场发展或新增热资格。见[实施报告](../../HYBRID_ROCKET_EULER_PREPARATION_V194.md)和[下一决定](NEXT_DECISION.md)。

- 事前[唯一优化计划](PLAN.md)、[只读恢复计划](RECOVERY_PLAN.md)及[原生准备计划](PREPARATION_PLAN.md)。
- [原优化失败](evidence/netgen01/report.json)与[原准确源码](evidence/netgen01/executed_optimize.py)；[首次恢复失败](evidence/netgen_recovered01/report.json)；[恢复通过](evidence/netgen_recovered02/report.json)。原始134159单元输出仍只在E盘，未重跑优化。
- [独立原壁/体积/质量审核](evidence/geometry01/report.json)；[最终批次选择回放和5拒绝](evidence/replay01/report.json)。原官方源码及许可证在[evidence/sources](evidence/sources)。
- 两次原生准备前的代码失败：[Linux专属导入](evidence/design01.stdout.txt)、[GBK读取](evidence/design02.stdout.txt)，相应准确源保留；[成功设计](evidence/rocket_design03/manifest.json)。
- [唯一原生准备](evidence/preparation01/report.json)、[父监督](evidence/prepare01_process.json)、[完整原生配置/初场审核](evidence/native_audit01/report.json)。均匀初场不是求出的压力/热流。
- [13项合同拒绝](evidence/contract01/report.json)，明确manifest替换的测试范围；[最终源版本](evidence/final_source)。
- [证据清单](evidence_manifest.json)：逐文件SHA、准确进程时间、仅本地大场路径及29保护项核对。小证据不等于大网格/完整历史的GitHub备份。

11父监督23.252秒，包含4次准确失败，未超各自预算。原生优化仅1次、原生准备仅1组，CFD推进/构建/MATLAB/主UI/computer-use均0。最终优化器新增批次选择逻辑仅回放已存原始数组测试，没有第二次优化。无工程包、无独有删除。
