# v94验证与交付索引

本版压力驱动表面流线/流管度量和条件研究工作台已完成实际计算、独立复核及浏览器核验。源码、测试、小型冻结输入和本目录证据随GitHub交付；不生成工程ZIP或重复大型副本。它尚未求解新的真实箭体压力、热流或壁温，合格热流面积保持66.5626327424%、本轮新增0；Q41仍为19次接受/82次主NS且未稳态，总目标active。

## 真实计算与来源

四个原MATLAB进程均terminal、exit0且运行期源码稳定，未为本轮交付收口重算：

| 回执 | MATLAB引擎PID | 本地输出 | 身份 |
| --- | ---: | --- | --- |
| `raw/surfaceflowunit94_process.json` | 50120 | `outputs/v94_validation/unit01` | 最初三射线版本，未覆盖 |
| `raw/surfaceflowdual94_process.json` | 10276 | `outputs/v94_validation/unit02` | 最终五射线/双间距核心；小型MAT随仓库保留 |
| `raw/surfaceflowview94_process.json` | 37284 | `outputs/v94_validation/study01` | 资源修复前研究，保留原数值 |
| `raw/surfaceflowview94final_process.json` | 31788 | `outputs/v94_validation/study02` | 修复后最终研究/工作台 |

最终研究MAT的SHA256为`193d4d0ceb6b8915744d76ba3ba3424138dca4a4ad4e97043c3f18874956a176`，只保存在本地。核心冻结MAT位于`tests/data/surface_pressure_flow/`；其实际进程和源身份可由当前独立审核器复核。来源PDF、实际输出及旧MAT没有因GitHub同步而自动备份或删除。

## 数值和界面核对

- `frozen_input_replay.json`：当前审核器直接读取仓库冻结输入，16组结果/80射线通过；RK误差比约16、16，标签误差比约4、4，解析聚焦s=0.8 m后保持未知。
- `final_study_replay.json`：当前审核器直接读取最终研究MAT/JSON，三组闭合箭体和31组流管/155射线/1271站点/1240线段通过；给定场量与数值方向/度量分列。
- `final_node_replay.json`：20个五步骤事件、20条无关/过时消息拒绝、164行CSV和未知输入不变零通过。Node不提供浏览器视觉或GPU证明。
- `browser_verification.json`：实际1280×720四模型、固定五步骤侧栏/独立滚动、三维点选、滑块End键与真实CSV下载记录。锥段管5第21站s=0.173185 m、h=0.0551727 m/rad、h/h₀=5.32963，输入p=71.5357 kPa、T=820.01 K；q与求解壁温未知。
- `cone_tube5_actual_download.csv`、`actual_download_audit.json`：真实下载逐行与41个保存站点核对，SHA为`1e0dbd2d245f4495b4cd56e28d1317fa78a25df707aaf2a49fd2ebc18cb8bf05`。

四张实际界面图为`cone_flow_probe.png`、`cylinder_helical_flow.png`、`spherical_nose_flow.png`和`vortex_pressure_turn.png`。`workbench_before_asset_fix.png`保留初次观察，不作为修复后的截图。

## 原错误与修订保留

- `raw/dual_spacing_revision.json`及`raw/source_epochs/`：三射线在聚焦点留下差分残宽，修订后用五射线双间距判定；单间距旧MAT和源版本仍在本地。
- `raw/asset_fix_binding.json`：资源正则误匹配JSON下载链接、隐藏选择器被CSS覆盖，均已修复并回归。
- `raw/frontend_numeric_equivalence.json`：前端修复前后115,295个保存数值/报告对象逐值一致；只排除MAT容器头。
- `raw/comparator_first_failure.json`、相应旧比较器：原比较器未递归MATLAB对象数组中的空字符数组，修正对象递归后通过，未放宽数值容差。
- `raw/compact_probe_patch.json`、`raw/frontend_missing_input_guard.json`：紧凑探针只改前端；未知输入换算不变零。
- `raw/*before_report_path*`、`audit_report_path_revision.json`：为冻结输入复核增加可选新报告路径；原审核器及原报告保留，数值检查条件保持，最终审核已再次通过。
- `browser_control_review.json`、`raw/browser_control.html`：同类MutationObserver错误在无项目脚本控制页复现；`browser_verification.json`中的原两条错误原样保留。没有逐栈断言每一原事件来源或声称全部控制台干净。

`manifest.json`逐项绑定本次源码、冻结输入、原回执、修订证据、截图和下载；`delivery_check.json`核对当前源字节、原运行终态、物理对象身份和本次交付范围。此类SHA核对不能替代物理或实验验证。

## 边界与后续

本版属于解析校核与条件数值研究，尚无新增独立实验验证。完整底图不代表全区域已解，解析压力/热状态不是本模块新求外缘，h也不能直接换算热流。下一步接内部合格外流状态与NASA7/变比热的一致采样合同，处理真实附着区域及边界层入口；分离/再附/底面仍独立推进。方法、复现命令与明确限制见[方法说明](../../PRESSURE_DRIVEN_SURFACE_FLOW_V94.md)。
