# 压力驱动表面流线的小型冻结输入

`pressure_surface_flow_evidence.mat`来自2026-10-29最终五射线核心实跑`surfaceflowdual94`，长度461,095字节，SHA256：

```text
45b8166707d566c072c706a719ed14365dcc13895831ecbc4bc9708d671da326
```

这里保存MATLAB实际积分结果、独立解析参考、加密结果、拒绝及有效前缀记录，用于无需MATLAB的独立复核。它是解析验证输入，压力/外缘状态为给定量，没有实验热流或真实20°整箭求解资格。对应实际进程回执在`docs/validation/v94/raw/surfaceflowdual94_process.json`。

在项目根目录运行，报告须使用尚不存在的新路径：

```text
python tests/auditPressureDrivenSurfaceFlow.py tests/data/surface_pressure_flow docs/validation/v94/raw/surfaceflowdual94_process.json --report outputs/surface_flow_core_recheck.json
```

依赖Python、NumPy与SciPy。复核检查16组结果/80条射线、曲面轨迹、双间距横向尺度、四阶/二阶加密、聚焦与NaN前缀，以及实际进程和当前MATLAB源身份。它读取冻结输入，不调用MATLAB，也不覆盖旧证据。完整方法与限制见`docs/PRESSURE_DRIVEN_SURFACE_FLOW_V94.md`。
