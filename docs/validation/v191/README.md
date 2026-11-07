# v191证据索引

结论及限制见[实现报告](../../HYBRID_EULER_CONTROL_V191.md)；下一项见[NEXT_DECISION](NEXT_DECISION.md)。本轮新热流面积与物理资格均为0。

- 事前预算：[准备/接口](PLAN.md)、[唯一400更新](RUN_PLAN.md)。
- [原生准备](evidence/preparation01/report.json)、[实际配置/拓扑](evidence/audit01/prepared_audit.json)、[同EOS理论参照](evidence/audit01/reference.json)。
- [原生19状态查询](evidence/control01/gas-query-v191.csv)、[物性核对及6拒绝](evidence/thermo01/report.json)。
- [真实运行监督](evidence/run01/native/report.json)、[原生输出](evidence/run01/native/stdout.txt)、[准入SHA](evidence/admission.json)。
- [最终审核](evidence/final02/report.json)、[240原控制面CSV](evidence/final02/wall_comparison.csv)、[科研图](evidence/final02/control_comparison.png)。
- [独立二进制/CSV、面积/总焓和12拒绝](evidence/verify01/report.json)。final01与final02数值及图片相同，仅最终导出绑定更严格的帧序号读取器。
- [清单](evidence_manifest.json)：118份小证据2,392,047字节，所有复制项逐字SHA核对；29保护项保持。局部输出4,287,493字节，不作工程包，无独有删除。

大场政策沿用：本次5个rawbinary路径仍仅本地（含逐字网格别名），路径/大小/SHA在清单`localOnly`中。准确运行库、WSL及原历史大场未备份到GitHub。其余小输入/配置/metadata/状态交换结果/源与回执公开。需要重建时使用固定Eilmer版本与v187包，严格依准备器先生成一个**新的**目录；不能覆盖历史路径或直接重跑已用admission。原唯一400步结果作为历史保留，不自动续算。

本机只读复核可以使用：

```powershell
D:\anacondapython\python.exe scripts/auditHybridEulerControl.py --case outputs/v191_validation/control01 --preparation outputs/v191_validation/preparation01/report.json --run-receipt outputs/v191_validation/run01/report.json --reference outputs/v191_validation/audit01/reference.json --out outputs/v191_validation/new_read_only_audit
```

上面只读取已保存场并创建新审核目录，不运行CFD。必须先确认`--out`不存在。当前场必须保持原SHA，无法从GitHub小证据冒充恢复独有原生场。
