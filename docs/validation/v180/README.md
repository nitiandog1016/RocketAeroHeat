# v180证据索引

本轮真实候选初始化失败，未产生新的非线性结果或物理资格。[交付说明](../../EILMER_FULL_DOMAIN_DEVELOPMENT_V180.md)、[初始计划](PLAN.md)、[唯一候选事前计划](REAL_RUN_PLAN.md)、[下一步](NEXT_DECISION.md)。

| 内容 | 准确证据 |
|---|---|
| 旧配置兼容 / 小域阶段 | [profile01](evidence/profile01.json)、[small01_phase](evidence/small01_phase.json) |
| 首次资源失败 / 修复 | [原失败](evidence/memory01/limit/report.json)、[修复检查](evidence/memory02/report.json)、[取消/超时](evidence/cancel01/report.json) |
| 真实准备 / 原运行回执 | [prepare](evidence/real01_prepare/report.json)、[旧误报回执](evidence/real01_run/report.json)、[原生日志](evidence/real01_run/process00/stdout.txt) |
| 当前真实失败结论 | [inspection](evidence/real01_inspection.json)；优先于旧`passed=true`解释 |
| 完成审核 / 旧7面接口 | [completion02](evidence/completion02/report.json)、[legacy_interface01](evidence/legacy_interface01/report.json) |
| MATLAB | [结果](evidence/matlab01/report.json)、[父监督](evidence/matlab_supervision.json)、[1885源锁](evidence/eilmer_development01_process.json) |
| 下一实数路线依据 | [准确源码有限审核](evidence/memory_route/report.json)，没有原生实数运行 |
| 新门槛兼容 / 源码追溯 | [Run14只读核对](evidence/closeout01.json)、[MATLAB时点源恢复](evidence/matlab_source_recovery.json) |
| 编码启动失败 | [bootstrap01](evidence/bootstrap01_failure.json)，原生任务0 |
| 时间 / 存储 / 字节绑定 | [budget](budget_receipt.json)、[evidence_manifest](evidence_manifest.json) |

原回执、原失败、运行时准确脚本和追溯SHA保持；最终源码在本提交中，不复制整个工程。`.log.raw`是原日志逐字副本，扩展名改变不改变内容。大于1MiB的大场/网格只列本地路径与SHA，不上传；依赖它们的测试不是仅克隆仓库即可离线完整重算。小输入、生成代码与固定引擎安装信息可用于重新准备，不能称大场已备份。

程序分项/嵌套时间分开，MATLAB锁检查只证明当时执行窗口内源稳定；后续完成审核名称修正等另有说明和测试。旧控制/v167研究/本轮小域均不提供新的真实全面积物理解。
