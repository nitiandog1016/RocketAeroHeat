# v102 已测试研究与实际UI验证

状态：数值内核及科研查看流程已测试；球鼻候选未通过真实区域资格，新增正式热流面积0。完整目标继续有效。方法与结果见[阶段说明](../../SHOCK_FITTED_DYNAMICS_V102.md)。

| 证据 | 覆盖的实际内容 |
| --- | --- |
| dynamics_final_audit.json | 非均匀Sod四面体体平均、Riemann、ALE/CFL、末步和全步五量账 |
| motion_audit.json | 全局共享节点投影、独立SVD/QR及可表示速度/GCL |
| cap_coarse_audit.json、cap_fine_audit.json、cap_resumed_audit.json | 三组真实曲面瞬态；明确列出未稳态、相容差与几何误差 |
| export_audit.json | 8算例16帧、56412场量及42CSV的独立SVD交线/KD树轮廓核对 |
| sections01/node_audit.json | 实际导出JS的流程、全部探针、CSV、消息隔离和保相机叠加 |
| browser/independent_download_audit.json | 实际浏览器下载的36段周向CSV数值与几何 |
| browser/ | 原生五步骤、3D拾取、滑块、初末帧、未知状态、方法/SHA、截图和实际下载 |
| provenance.json、raw/processes/、raw/logs/ | 9项实际终态/锁源记录，883精确来源版本，12显式输入 |
| legacy_ui_node.json、legacy_export_audit.json | 旧v101模式兼容与旧实际MAT映射 |
| figures/ | 从既有数值证据生成的定量图，没有新增数值求解 |

raw/dynamics01保留报告结构问题前的源码/实际日志；raw/cap01保留逐节点平均的原始问题版本。raw/export_audit_attempt01/02为独立读取器字段和状态标签修正前原文。raw/provenance_attempt01保留全快照中未调用模块缺版本的原始报告；实际motion调用闭包不含该求解器。raw/provenance_attempt02保留Git换行原字节归档前的来源索引，最终以provenance.json及raw/source_versions为准。以上均不能改写成区域资格失败或覆盖掉原数值。

完整MAT仍在`outputs/v102_validation`。当前cap04、cap03续算输入、早期独有结果与检查点保留；本次没有再扫描/清理60GB历史研究区。源码、必要的小型复现输入及验证资料按GitHub约定交付，没有生成工程包。现有MAT与来源身份绑定，不可原位改写；复现或继续研究使用新目录和相匹配的新来源。
