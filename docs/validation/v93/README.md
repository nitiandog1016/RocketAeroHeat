# v93 验证与交付索引

状态：数值与实际浏览器检查已完成；源码、测试和本目录小型证据随本版GitHub提交交付。最新Q41仍未稳态，正式合格热流覆盖66.5626327424%，新增合格面积0；全面积及全周期目标未完成。

| 范围 | 实际证据 | 结论及边界 |
| --- | --- | --- |
| 候选路径、小域稠密参考 | [路径审计](coordinate_unit01/independent_trial_coordinates.json) | 两排列通过；人工面和制造矩阵不等于闭合NS几何或物理Jacobian |
| 原默认兼容 | [原路径](default_compatibility.json)、[行预处理后](conditioned_default_compatibility.json) | 5,658项已有非计时对象逐值不变；只排除耗时/进程号 |
| 精确行预处理 | [小域](row_unit02/independent_row_conditioning.json)、[真实首投影](row_real02/independent_row_conditioning.json) | 全部条件保留；80项精确有理关系、完整P转置跨度通过；不证明近重复差异的物理意义 |
| 原失败与恢复 | [秩诊断](projector_rank_diagnostic.json)、[恢复绑定](recovery_binding.json)、[逐值身份审计](independent_recovery_identity.json) | 原exit1保留；自由方向/坐标完全相同，复用128步、新受限256步，准备/探针/主NS为1/2/4 |
| 真实当前方向 | [原H、目标和投影审计](continuous01/independent_temperature_active.json) | C未达1e−3，完整模型/原NS缺陷约0.05432通过0.1；仍须实际线搜索 |
| 单元/内部界面 | [独立复核](independent_cell_interface_continuation.json) | 全部CV、输入分区和并集通过；公开界面抽样289/607/607，不声称全部独立Riemann重解 |
| 墙面、五量与更新 | [Q41独立审计](krylov41_independent_audit.json) | 原墙构成、Fourier、Jacobian、散度及实际候选算术通过；不是独立完整NS求解器 |
| 全部开放边界 | [储库审计](krylov41_independent_riemann_reservoir.json) | 19,752开放面及墙/散度/端点/历史通过；不代表稳态或实验验证 |
| 最终工作台数据与事件 | [Node](workbench01/node_wall_velocity_summary.json)、[原面字段](workbench01/independent_wall_velocity_view.json) | 592,560项CSV值、6,624截线样本、15同窗事件及滚动重置通过 |
| 实际浏览器 | [检查记录](browser_verification.json)、[压力](q41_workbench.png)、[热流](q41_heat.png)、[同窗诊断](q41_diagnostics.png) | 1280×720视口、固定五步骤/独立滚动、面元8901、算法/温度身份、任意截面与色标已实际检查 |
| 源码运行身份 | [796项MATLAB源码](recovery_source_stability_final.json) | 当前MATLAB源与真实恢复任务一致；计算结束后另改UI，无新物理重算 |
| 工程/实验资料 | [来源记录](research/source_triage.json)、[选择说明](../../ENGINEERING_SOURCE_TRIAGE_V93.md) | 已读原PDF，未增加区域闭合或实验误差验证；M6主分布不能冒充20°测点 |

[汇总与文件SHA](progress_coordinate_preview.json)记录实际Q41指标、任务终态、源码与证据身份。[RMS/峰值图](temperature_active_continuation.png)及[逐轮温域图](temperature_repreview_rounds.png)对应真实记录，亦提供同名PDF。原未修复侧栏截图和旧色标/布局前回执保留，其资源SHA不得用于宣称最终页面通过。

原失败投影源码保存在`source_epochs/`。核验器的两次字段/数组处理失败及修正前源码另存，未修改MAT数据或放宽物理/数值门槛。原完整MAT、自由/受限检查点、状态和大矩阵不在GitHub；Q40旧准备矩阵已按用户授权清理，完整历史矩阵审计需要新目录重建并建立新来源绑定，不承诺历史MAT SHA相同。独立恢复身份与最终工作台审计在清理后仍通过，未重复准备或求解。

本版不是全表面物理解验收。后续优先推进压力驱动的表面流线/流管接口与工程区域衔接，同时完善科研UI；内部NS若继续则从Q41建立新系统。累计热量、轨迹与跨流态扩展仍后置。
