# v173 证据索引

本轮节点状态保护/空间壁温组合已构建并测试，实际4次试探恢复；唯一实际候选异常退出、残差收益失败，不采纳/不续算/不改默认，新增物理资格0。

- [事前计划](PLAN.md)、[准确覆盖](COVERAGE_REVIEW.md)、[组合评审](COMBINATION_REVIEW.md)、[构建审查](BUILD_REVIEW.md)。
- [终止证据](FAILURE_REVIEW.md)：退出原因未知，4条恢复事件不是最终异常定位。
- [原生推进路径限制](NEXT_ROUTE_REVIEW.md)、[官方低焓参考路线](REFERENCE_ROUTE_REVIEW.md)、[下一决定](NEXT_DECISION.md)。
- [完整说明](../../NACA_SPATIAL_STEP_V173.md)、[存储回执](storage_receipt.json)。

小证据由`evidence_manifest.json`列出原路径、归档路径与SHA256。重复内容引用既有v171/v172证据，不复制旧证据树或运行包。`finalSources`绑定最终实现，实际失败时的准确源码单独保存；大场、引擎exe/对象只记身份，不在GitHub备份。

主要本地报告位于`outputs/v173_validation/`：`coverage_review01`、`fixture01`、`build01`、`contract01`、`small01`、`audit_small01`、`verify_small01/02`、`compare_small01/02`、`real01`、`audit_real01`、`development_real01`、`verify_real01`、`termination_review01`、`screen_real01/02`及`history_plot01`。父目录保留监督回执与日志。第一次覆盖启动在创建进程前失败，首次小场清点缺元数据、逐字比较失败及筛选序列化失败均保留。

关键身份：原生exit134；正常任务/正常完整保存未通过；仅完整local0/系列1502的严格保存帧接口通过，最后history56/系列1558无场。验证器passed=true只表示复核完成，不是求解或实验通过。最终Tw/y+及两完整帧p/q变化待验，不能用首帧替代。

`history_plot01/spatial_step_history.png`和SVG已实际查看，显示57行原生历史、四分之一域有符号功率和4条事件，明确最终原因未知及局部/系列身份。没有新MATLAB/主UI或整箭实算。
