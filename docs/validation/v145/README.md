# v145 验证索引

一般混合网格、独立壁距初场、真实完整域SST及主五步骤已接通；物理收敛/热流筛选未过，新增物理资格0。

- [事前计划](PLAN.md)、[小例修正前补充](FIXTURE_ADDENDUM.md)、[解析/穷举/置换检查](unit/report.json)。
- [主运行](evidence/phase/main_supervisor/process.json)、[真实原生终态](evidence/actual/engine_process.json)、[同阶段热审核](evidence/actual/sst_state_audit.json)。
- [独立原面/壁距/偏移列复核](evidence/phase/independent_review.json)、[k下限传播](evidence/phase/startup_diagnostic.json)。
- [只读主工作台复验](evidence/phase/main02/report.json)、[小例独立能量/距离](fixture02_independent.json)。
- [存储回执](evidence/phase/storage_receipt.json)、[完整来源与文件清单](delivery_manifest.json)、[结果与下一步](../../GENERAL_HYBRID_INITIALIZATION_V145.md)。

`fixture01_sources/`、`main01_sources/`、`reviewer01_sources/`、`diagnostic01_sources/`保留修改前准确版本；失败身份不覆盖。只读复核修正、重开等待与存储合并没有数值重跑。`source_versions/`保留本轮最终源码及相关共享依赖。

公开证据包含小型混合网格及小型实际输入输出。真实大MAT/SU2重启动/完整字段仅本地；本轮没有为GitHub生成工程包。当前工作树既有的OpenCode UI改动未混入提交；清单记录其哈希，准确副本保留在本地`outputs/v145_validation/uncommitted_ui_sources/`，不能声称本次公开提交独自包含全部已测依赖。
