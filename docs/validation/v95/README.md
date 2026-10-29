# v95验证索引

本版交付同物性压力恢复、总状态锚点转换、连续cp积分及对应三维工作台。详细方法和限制见[正式说明](../../THERMAL_PRESSURE_SURFACE_FLOW_V95.md)。实际主计算、独立复核及最终浏览器检查通过；新增合格热流面积0，未完成单时刻全面积目标。

## 数值与来源

- `frozen_input_independent_audit.json`：对仓库内337730字节核心MAT和原失败MAT重新执行的独立复核，六组通过。
- `raw/original_thermal_audit.json`：原`thermal03`独立复核；包括独立热力积分、横向动量、法向激波根和Taylor–Maccoll。明确不是实验验证。
- `raw/caloric_regression_audit.json`、`raw/caloric_exact_equivalence.json`：旧16组/80射线回归和32672个保存对象完全相等；仅新增物性模式元数据。
- `raw/historical_v94_frozen_replay.json`、`raw/legacy_v94_data_with_current_program.json`：旧冻结数据按原Git来源核对及当前JS兼容检查，不冒称当前MATLAB重算。
- `raw/final_view_audit.json`、`raw/final_view_export.json`：最终view05的3模型/4流管/124站点、真实锥面、几何拓扑、相机和MAT→JSON逐值映射。
- `raw/final_workbench_node.json`：15步骤、15无关/过时消息拒绝、93CSV行与逐变量方法/状态核对。
- `raw/*_process.json`、`raw/*_stdout.log.txt`：两次真实失败、最终主计算和各次仅重导出界面的原进程证据。均已终态，按各自实际exit0/exit1保留。
- `raw/source_epochs/`：连续物性、标签、相机、阶段视角及CSV方法列修订前的小型源码原字节。

## 实际浏览器与下载

- `browser_verification.json`：最终view05同页五步骤、低温末站和真实三维锥面站点点选；探针Te/Tr/Tw及方法身份分列。
- `thermal_surface_hot_plane.png`、`thermal_surface_cold_plane.png`：最终程序下跨1000K和150K的Te曲线；均为条件反算。
- `thermal_surface_cone.png`：实际点击第2流管、第18站，s=0.17m，h/h0=5.25。
- `thermal_surface_cone_probe_methods.png`：对应探针完整数值、各变量方法和未知热量状态。
- `thermal_surface_cone_tube2.csv`、`downloaded_csv_audit.json`：实际浏览器下载的31行CSV；所有数值逐值相同，方法与状态列逐项相同，热流和求解壁温为unknown。

旧图`before_nested_role_fix.png`、`conditional_cone_temperature.png`、`conditional_cone_probe_final.png`分别属于标签/相机/CSV方法列修订前的中间版本，只作为修订证据，不是最终交付截图。历史文件名中的final不能替代此处版本说明。

当前浏览器新标签页未记录警告/错误；此前v94浏览器控制页的原记录仍在原验证目录，没有改写。Node和独立数据核对不能代替实际浏览器检查，两者证据分列。

`SOURCE_INDEX.json`记录本次交付相关源码与小型证据SHA。大型输出仍位于被忽略的`outputs/v95_validation/`，没有复制完整工程或生成ZIP。维护提交`1cf7cbd007f6062d333d0fcf63c9554ec02c7ddd`另行清理了约1.95GB旧QR缓存，保留全部v95数据。
