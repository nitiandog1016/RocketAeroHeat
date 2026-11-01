# v146 验证索引

本轮原19752面完整域SST二阶控制/一阶启动候选各100更新；没有新几何或小型格式求解。候选恢复筛选失败，不运行恢复二阶，不更改默认，新增物理资格0。

- [事前计划](PLAN.md)：固定共同源149、唯一重构开关、100更新预算和恢复条件。
- [固定源码依据](source_basis.json)：官方提交7048fd35的归档核验、原文件及COPYING许可。
- [配置与恢复条件测试](contract_tests.json)：2项测试通过；包括错误源、CFL、门槛、湍流阶数和恢复否决反例。
- `evidence/control/`、`evidence/first_order/`：实际配置、输入清单、真实终态、七方程历史、热接口审核。
- `evidence/startup_comparison.json`：整域节点、原面面积、数值阶数和预设筛选；`restoreAllowed=false`。
- `evidence/workbench01_process.json`、`evidence/workbench01/report.json`、`evidence/export_review.json`：主工作台只读和独立原面/截面导出核对，实际完成后归档。
- `source_versions/`和`delivery_manifest.json`：准确执行源码、文件SHA、预算口径和本地大场索引，归档后才构成交付证据。

主要运行入口为`prepareSU2SSTStartupStudy.py control/first_order`，接既有`runEmbeddedSU2.py`、`auditSU2SSTState.py --actual-surface`及`compareSU2SSTStartupStudy.py`。准备器拒绝覆盖已有任务；这些命令依赖本地v145源场，仓库没有上传大型重启动文件，不能声称仅克隆即可复现本次真实整箭结果。

MATLAB完整替换文件在项目根目录及本轮源码快照中。现有OpenCode的WebGL改动单独保留在工作树；不混入提交，其准确版本仅本地保存并列明散列。无新原生鼠标/截图验收。
