# v106证据入口

结论见[阶段说明](../../SU2_ROCKET_INTERFACE_V106.md)。`PLAN.md`是运行前预算和判据，`independent_audit.json`是独立复核，`delivery_manifest.json`绑定真实终态／准确字节。

- `processes/`：8项MATLAB任务；5成功、3失败，不抹去超时、尾轴节点错误或日志正则错误。
- `cases/`：6项SU2进程。`rocket_volume02`两次是虚拟网格dryrun；真正加载三维网格的是`rocket_startup01/angle01`和`rocket_startup02/angle02`，只有一步，未收敛。
- `mapping/report.json`：全部2,448原面映射、独立的标量面积热流与引擎节点热负荷口径。
- `inputs/`：一份三维网格及外域对照网格；相同几何的不同攻角只改变配置，不必重复网格副本。
- `source_versions/`：当前与发生失败时的准确源版本；未变的v105源版本可由清单引用上一目录。
- `plate_outer_audit.*`：外域数值检查；加密试验没有壁面结果，未伪造同类误差文件。

本地依赖还包括v105 plate03、v106生成几何及投影MAT、restart与完整历史，它们没有全部上传。按源版本和小型输入可重算，不保证逐位相同浮点结果。不重建已清理的旧矩阵，也不删除当前失败历史或共享对象库。

此目录不是完整物理验证或主UI验收。面积映射完整不等于全面积气动热完成；面积口径差、未收敛、三维近壁质量、公开实验与突起支持仍须按目标推进。
