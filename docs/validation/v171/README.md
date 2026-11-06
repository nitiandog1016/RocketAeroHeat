# v171 验证证据导航

- [事前构建、小组合与唯一实际任务计划](PLAN.md)
- [准确observed基线与有限构建](RUNTIME_REVIEW.md)
- [入口、稀疏保存与旧合同审查](CONTRACT_REVIEW.md)
- [空间壁热与原实验段功率](HEAT_REVIEW.md)
- [实际发展、边界诊断及下一步依据](DEVELOPMENT_REVIEW.md)
- [原输运表的有限补证与超时记录](TRANSPORT_REVIEW.md)
- [实现与实际结果](../../NACA_SPATIAL_WALL_V171.md)
- [下一执行起点](NEXT_DECISION.md)
- [逐字证据清单](evidence_manifest.json)
- [仅重复别名去重回执](storage_receipt.json)

小证据保留准确源、实际输入、构建对象/源码SHA、监督回执、失败、热审核、发展统计、科学图。清单按SHA复用相同存档，`archive`字段是实际存放位置；文件数量/实际存档字节见清单，不把路径数量当独有文件数量。

主要报告为 `audit_real02/report.json`、`development_real01/report.json`、`saved_frame_study02/report.json`，末科学图为 `plot_real03/naca_spatial_fields.png`。它们的源路径位于 `outputs/v171_validation`，公开副本按清单的archive查找。原生任务超时，正常成功为false，4完整保存帧的热接口通过；这些状态不能合并成“本次求解成功”。

全部大网格、初场/体场/原生表面、exe/DLL/对象库和原PDF留本地，仅记录身份，未备份到GitHub。独有结果未删；两末帧别名去重后原路径和SHA保持。无工程包、无新整箭/MATLAB/原生UI，实验与新增物理资格0，总目标继续。

后处理墙钟322.922秒超过原180秒计划，6进程实际累计92.156秒；资料补证334秒超过300秒。两项偏差如实保留，没有延长原1800秒数值预算或重复CFD。
