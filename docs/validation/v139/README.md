# v139 验证索引

唯一固定边界小例通过；唯一真实箭体候选质量失败。新物理解面积0，无SU2/MATLAB/新UI测试。

- `plan.md`：执行前预算、选择规则、门槛与停止条件。
- `runtime/`：独立保点入口构建与准确库/源码身份，未重新编译上游MMG。
- `seed_test/`、`seed_test_process/`：已知体积小例输入、完整细分谱系、原生输入输出、日志与通过结果。
- `actual/`：真实几何、质量拒绝、细分事件、独立质量/谱系复核及原生日志；大MAT留在`outputs/v139_validation/actual/`。
- `actual_process/`、`review_process/`：完整候选及只读复核的监督回执；原生exit0与候选质量exit1分开。
- `source_versions/`：本轮准确执行源码，`storage_and_sources.json`：哈希与本轮文件长度，不冒称磁盘实占或完整历史备份。

总体解释与下一步见[阶段说明](../../MMG_SEED_OPTIMIZATION_V139.md)。未通过质量的候选未接入流场；不得清除质量失败标志或使用本轮网格自动续算。
