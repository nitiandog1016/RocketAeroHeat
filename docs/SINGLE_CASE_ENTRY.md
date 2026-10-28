# 显式单工况计算入口

`runAeroHeatSingleCase(mesh,cfg,progressFcn)` 供统一科研工作台的“计算”步骤调用。它不创建窗口、不改变网格、不触发导出；返回的结果可直接交给 WebGL 查看器，切换结果页无需重算。

```matlab
[current,freestream,report] = runAeroHeatSingleCase(mesh,cfg,@onProgress);
```

输入网格须先完成 `calcFaceGeometry → fixNormalDirection → buildSurfaceCoordinates → analyzeLocalGeometry`。当前边界条件使用 `cfg.flight`；长度和高度为 m，温度为 K，攻角/侧滑为 deg。与专门用于飞行历史每站初始壁温的 `calcFlightStationState` 不同，本入口保留 `cfg.heating.wallTemperatureByFace`，未知面元壁温仍为 NaN。

入口依次执行自由来流、区域分类、修正牛顿/所选底压、压力载荷、驻点、内部工程外流、轴向激波/膨胀、准流线、热边界和面积审计。模型、阈值、未解处理和算法来源沿用既有后端。目前完整链要求 M>1；工况不满足时拒绝执行，不返回部分压力作为已完成解。

进度回调接收一个结构体，包含 `phase`、`status`、`fraction`、`message`、`elapsedSeconds`、`errorIdentifier` 和 `errorMessage`。11 个阶段依次为 `inputs/freestream/regions/pressure/loads/stagnation/externalFlow/axisymmetricMarch/streamlines/heating/coverage`。`status` 为 `started`、`completed` 或 `failed`；`fraction` 表示已完成阶段数/11，不是物理解面积、预计耗时或阶段内部进度。失败时发出真实阶段通知并重抛原错误；通知自身失败不能遮盖后端错误。调用方只在成功返回后接受快照。

`report` 保存实际接受的 `cfg/condition`、面元数、总面积、压力求解/给定/未解面积比例、热流求解面积比例、载荷状态、耗时和阶段事件。给定 Cp 不计入已求解压力面积；“链执行完成”也不等于全面积求解完成。每个区域和变量的适用范围继续读取结果中的真实算法、状态与来源。

验证由 `tests/testRunAeroHeatSingleCase.m` 使用真实经典锥柱工程链检查：当前非均匀/未知壁温、压力与载荷同步刷新、旧历史清除、亚声速/旧几何/晚期未知模型的真实失败阶段与原错误保护。状态以实际测试日志为准；该测试不构成实验精度验证。
