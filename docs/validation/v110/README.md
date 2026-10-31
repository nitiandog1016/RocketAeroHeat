# v110验证证据

范围和实际结论见[专题说明](../../SU2_CHECKPOINT_DIAGNOSTICS_V110.md)，运行前计划见[PLAN.md](PLAN.md)。2项SU2实际任务终态：1配置失败、1诊断迭代成功；5项MATLAB实际终态：4失败、1最终成功。失败不能从统计中删去。新增物理资格面积0。

`processes/`保留真实MATLAB终态与stdout；`engine/probe01`和`probe02`保留实际配置、回执和短历史。`checks/`含每轮已有审核及最终主UI/导出，`source_versions/`绑定当前及失败源码。`independent_audit.json`独立词法读取59,122节点并核对RMS、时刻及CSV元数据。`delivery_manifest.json`列全部SHA、来源、终态及本地存储大小。

原生体场、初场和重启动仍在`outputs/v110_validation/`及所引用v108/v109目录，未随GitHub备份；不能仅凭公开小型报告重做全部数值审计。共享对象库未改动。没有工程包或重复工程副本。

测试反例的进程回执明示`syntheticTestFixture`，不是新求解任务；只上传反例结论和构造源码，不把测试夹具冒称真实数值输出。截图来自v107已保存主窗口快照，不代表v109/M5三维场。完整研究输出未重新求解。
