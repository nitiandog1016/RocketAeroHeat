# v95小型冻结输入

`thermal_pressure_surface_flow_evidence.mat`为实际MATLAB任务`thermalContinuous95`保存的六组结果，337730字节，SHA256：

`adf156e52224a72e5abe3d0a19bca99e383c4fe3d78f12f22d62d1740b8e4aa8`

`failure/thermal_case1_before_invariant.mat`为修正前1000K分段问题的真实诊断输入，49094字节。它用于保留失败与修正后守恒误差的对照，不能标记为通过结果。

来源回执、匹配源码、独立复核和实际浏览器证据位于`docs/validation/v95/`；方法、复现命令与限制见`docs/THERMAL_PRESSURE_SURFACE_FLOW_V95.md`。可直接进行Python独立复核，或用MATLAB从此MAT在新目录导出工作台，不需要重跑原六组积分。

本数据是条件流线核验，不是整箭气动热解，不含新求得热流或壁温；不是实验测量。
