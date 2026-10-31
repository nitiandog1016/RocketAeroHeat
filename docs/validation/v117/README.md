# v117实际证据

状态：冻结方向诊断及工作台身份已测试；恢复候选2/4未过启动门槛，没有600步执行；新增合格面积0。

- `PLAN.md`：运行前机制和停止条件。
- `build/`：固定受控编译/单位检查进程、准确源码清单、补丁、安装身份和合同测试。
- `cases/diagnostic/`：两步真实SU2终态、原输入合同、逐节点实际守恒增量、计数与具名输出的SHA。
- `checks/diagnostic_gate.json`：原节点完整状态匹配、独立终点内能及2/4失败决定。
- `checks/launch_gate_rejection.json`：600步准备实际被拒绝，未建候选目录、未启动求解器。
- `checks/review01/`：生产WebGL组件结果、原面CSV、截面CSV和运行时来源。
- `checks/independent_audit.json`：原面、EOS、CSV、具名第1步及四个冻结节点独立核对。
- `processes/`：MATLAB一次终态exit0与全源码稳定记录。
- `source_versions/`、`delivery_manifest.json`：本轮准确版本与逐项SHA，不因未来源码演进改写历史。
- `prepare_failed01.*`：UTF-8准备失败原文；修复后实际运行不覆盖此失败。

原51,547节点体网格、初场/重启动及运行时仅在本地`outputs/v117_validation/`和`tools/runtime/`；GitHub不代表大数据已备份。网格、初态、投影复用硬链接，不得原位写回历史输入。构建可再生依赖和后续回收见`storage_receipt.json`。本轮不改OpenCode或共享历史对象库，无工程包。

详见[v117方法与结论](../../SU2_FROZEN_ENERGY_DECISION_V117.md)。
