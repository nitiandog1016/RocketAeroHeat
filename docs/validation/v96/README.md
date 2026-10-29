# v96验证索引

本版交付逐入口标签H/S、物理横向梯度和入口分辨率检查；完整方法及限制见[正式说明](../../LABELLED_PRESSURE_SURFACE_FLOW_V96.md)。它是已测试的条件流线模块，不是整箭气动热解，新增合格面积0。

- `frozen_input_independent_audit.json`：对仓库内最终MAT重新执行的八组/40射线独立核验，包括物理梯度、Crocco关系、独立斜激波根、入口分辨标记与NaN前缀。
- `frozen_resolution_exact_equivalence.json`：门槛前后两份冻结MAT的18378个公共对象逐值相同，排除且只排除两项新增元数据。
- `raw/final_numerical_audit.json`：原最终`entry02`数值核验。
- `raw/final_view_audit.json`、`raw/final_view_export.json`：最终view03三模型/四流管/84站点/80线段的MAT→JSON→三维映射及几何/资源检查。
- `raw/final_node_workbench.json`：63行、48列模式的物理梯度/来源/分辨率字段、9次新增曲线选择、五步骤和无关消息拒绝。符号零归一化明确计数，非零值没有容差。
- `raw/legacy_v95_current_program.json`、`raw/legacy_v94_current_program.json`：当前JS和测试器复核旧数据，分别93/164CSV行通过。
- `raw/*_process.json`、`raw/*_stdout.log.txt`：两次数学核验及三次仅重导出UI的真实终态、实际引擎与运行期源码身份，全部exit0/源稳定；没有新增原NS。
- `raw/source_epochs/`、`raw/view01_check_notes.json`、`raw/csv_signed_zero_test_note.json`：文案/比较器及入口门槛前源码与阶段说明。早期UI断言失败不是MATLAB数值失败。

## 实际浏览器与下载

- `browser_verification_final.json`：最终view03五步骤、三种新增曲线、三维点选第2管第11站（s=0.1m）及明确的未分辨总焓变化提示。
- `labelled_cylinder_final.png`、`labelled_cylinder_probe_final.png`：最终三维工作台及总焓/熵/梯度、分辨状态、Te/Tr/Tw与方法探针。
- `labelled_cylinder_tube2_final.csv`、`downloaded_csv_final_audit.json`：实际下载21行48列；32数值列和16文本/布尔列逐项对应保存数据，热流/求解壁温unknown。
- `raw/final_surface_flow_payload.json`：真实下载对照的完整小型数据；复核命令见正式说明。

未带final后缀的三张截图、CSV及`browser_verification.json`来自入口分辨率提示加入前的view02，属于阶段证据，不能当作最终版本截图。最终物理数值仍与该阶段逐值一致。

最终浏览器记录保留一条MutationObserver参数类型诊断。其签名与[v94控制页记录](../v94/browser_control_review.json)相同，但本次没有项目URL/完整栈，因此不宣称全部归因或控制台干净。实际交互和数值校核结果分列。

`SOURCE_INDEX.json`绑定最终计算叶函数、验证程序及本目录小型证据SHA；MAT冻结输入位于`tests/data/labelled_surface_flow/`。大型本地状态未入库，没有生成ZIP或完整隔离副本。当前正式热流66.5626327424%、新增0，目标active。
