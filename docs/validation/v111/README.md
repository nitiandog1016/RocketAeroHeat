# v111 验证资料

结论与下一起点见[专题说明](../../SU2_NORMALIZATION_CHECKPOINT_V111.md)。总目标未完成，新增合格面积0。

- `PLAN.md`：运行前计划及有记录的调整；初版按原字节保存在`source_versions`。
- `normalization.json`：旧3750检查点经v110一次新评价的残差，按重构控制体积区分三种范数。归一化源码初版SHA绑定其原始版本。
- `comparison.json`：解析几何校核、三个实际终态、两项否决、全原面映射、正状态和双面积口径。`farfieldInitialCharacteristicCounts`是几何法向Mach分类，未经原生边界算子核验。
- `engine/`：三项真实SU2进程回执、配置、完整历史及小型表面文件；只上传基础表面文件，具名副本的原SHA在回执中保留。
- `checks/`及`processes/`：MATLAB两项实际任务终态、映射结果、writer配置及拒绝反例。反例使用内存回执／缺基础文件夹具，没有启动新求解。
- `failures.json`：失败筛选与两次Python代码问题；不把SU2退出0称为物理通过。
- `source_versions/`：所用相关源码准确字节，包含初始计划／准备器和两次Python失败版本；并非整个工程副本。
- `upstream/`：官方定义摘录及许可证来源。
- `duplicate_merge.json`、`storage_audit.json`：8组重复文件合并和独立SHA／文件身份核对，所有路径保留。测量含历史输入硬链接，独立文件长度不是磁盘分配量。
- `delivery_manifest.json`：文件SHA、使用版本、原输入输出绑定、任务终态及交付范围。体场、完整重启动、大几何不在GitHub；需要保留本地依赖才能重新审计，不能直接运行历史准备器覆写原目录。

本轮没有主UI布局修改或浏览器截图。具名检查点调用选项和审核逻辑已实测，主UI开关／超时检查点受控重开仍待实现。保存结果均未通过原残差门槛；不能把完整云图和逐面映射当作全面积物理资格。
