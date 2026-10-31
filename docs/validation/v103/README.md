# v103 验证索引

已测试，纳入本次GitHub交付；所有候选生产资格false，新增正式热流面积0。方法、实际数值、界面、限制和复现见[阶段说明](../../SHOCK_FITTED_INTERFACE_V103.md)。

| 记录 | 内容 |
| --- | --- |
| `baseline_audit.json`、`cap01_audit.json`至`cap05_audit.json` | 独立三维通量、QR/Brent节点、GCL与全步五量账 |
| `comparison.json` | 同壁面/气体/上游、前四组初态逐字节相同及六组实际指标 |
| `boundary_flux_summary.json` | 9个解析移动面/坐标变换、旧六步及两种运动续算 |
| `trace_reconstruction_summary.json`、`nodal_motion_summary.json`、`nodal_dynamics_summary.json` | 原数学/真实状态试探与未合格身份 |
| `export_audit.json` | 68148场量/诊断、当前帧壁压、时间及独立真实截线核对 |
| `sections01/node_audit.json` | 当前页11668探针、36份CSV、42列字段及类型5外壳 |
| `legacy_compatibility/node_audit.json` | 原v102数据与当前界面资源兼容 |
| `compatibility_source_binding.json`、`compatibility_runner.js` | 仅替换两个JS资源路径的原测试绑定，无数据副本 |
| `browser/interaction.json`、`browser/native_download_audit.json` | 实际浏览器与真正下载文件的独立数值核对 |
| `provenance.json`、`raw/source_versions/` | 8实际进程、17输入、49精确来源版本的保守调用闭包 |
| `raw/processes/`、`raw/logs/` | 原运行回执和日志，不改写退出码 |
| `raw/before_riemann_boundary/` | 新通量以前的精确源码原文 |
| `raw/compatibility_preflight/` | 测试包装器执行前路径错误及原文，无数值重算 |

完整MAT只存本地`outputs/v103_validation/`。旧六步小型输入随源码保存，按原SHA核对。未知热流、Te/Tr及壁温没有补零；线性拟合、守恒账和通量缺陷零值均不能替代区域资格或实验验证。此阶段没有建立新的实时全面积主程序，完整目标继续有效。
