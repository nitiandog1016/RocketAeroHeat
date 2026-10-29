# v97 共同激波入口原型：验证索引

共同三维激波上的守恒入口和原面元几何关联已测试。完整区域求解器仍为草稿：形状未反算、跨层质量/熵输运和壁面匹配未完成，未提供新壁面压力、热流或壁温。原科研UI没有改动，本轮没有新GUI验收或实验验证。正式合格热流面积仍66.5626327424%，新增0；Q41仍未稳态。

## 实际任务

| 任务 | 实际输入及结果 |
| --- | --- |
| commonShockEntry97 | 原v56真实箭体快照，仅几何和自由流；MATLAB引擎18132，terminal/exit0/全部来源稳定 |
| shockBodyGeometry97 | 上项已保存MAT，未重算跃迁；引擎39032，terminal/exit0/全部来源稳定 |
| commonShockFrozen97 | 334167字节小型几何/自由流输入，独立重算；terminal/exit0/全部来源稳定 |

任务原回执与原stdout在`raw/`，数值原文保存在`tests/data/shock_layer_entry/`。原18MB输入与原值的绑定见其中`input_binding.json`；小输入不含旧外流或热流，不能当作原MAT字节副本。

## 独立校核

- `raw/initial_shock_audit.json`、`raw/frozen_shock_audit.json`：完全气体闭式、独立Brent压缩根、独立五点差分及QR物理梯度。18共同曲面点、4平面激波；最大相对状态差1.0962438042e-14、导数差1.6419540994e-10。分别核对1126/1130项实际来源SHA。
- `raw/body_mapping_audit.json`：另一种平面交点与重心方程，9条法线的最近原面号一致，最大位置差2.8609792491e-17m。未知/内部起点保持NaN，几何映射不赋物理量。
- `raw/frozen_input_equivalence.json`：原输入和小输入重算的1745个保存对象逐值完全相等，仅排除MAT容器元信息；没有使用数值容差。两个MAT文件SHA不同不影响载荷相等，不宣称容器SHA复现。
- `raw/comparison_note.json`：两次比较器类型处理失败及部分证据复制中止的说明。问题在验证工具的列表/嵌套MATLAB结构处理，修正后独立通过；没有改MAT、重算物理场或放宽容差。
- `raw/copied_evidence_manifest.json`：小型冻结数据和原回执复制后逐字节SHA核对。

## 从GitHub源码复现

在工程根目录，输出路径必须尚不存在。只读重放已有证据，无需加载旧大型MAT或运行MATLAB：

```text
python tests/auditShockSurfaceInlet.py tests/data/shock_layer_entry/shock_surface_inlet_evidence.mat docs/validation/v97/raw/commonShockEntry97_process.json --report outputs/shock_entry_audit_new.json
python tests/auditShockNormalBodyIntersections.py tests/data/shock_layer_entry/shock_surface_inlet_evidence.mat tests/data/shock_layer_entry/shock_body_geometry_evidence.mat docs/validation/v97/raw/shockBodyGeometry97_process.json --report outputs/shock_body_audit_new.json
```

重算需要MATLAB，不需要外部CFD数据。完整函数文件带中文模块用途注释：

```matlab
setupEngineeringShockLayer; addpath('tests');
testShockSurfaceInlet('outputs/shock_entry_new', ...
    'tests/data/shock_layer_entry/actual_rocket_geometry_and_freestream.mat');
testShockNormalBodyIntersections('outputs/shock_body_new', ...
    'outputs/shock_entry_new/shock_surface_inlet_evidence.mat');
```

对新结果调用相应Python审核器，并传对应实际进程回执。若比较两个重算输入的全部MAT载荷，可用：

```text
python docs/validation/v97/compare_frozen_input.py ORIGINAL.mat RECOMPUTED.mat --report outputs/shock_input_comparison_new.json
```

审核器的来源绑定按实际MATLAB快照检查；后续若修改绑定源码，应建立新任务和新回执，不能将旧回执当作新源码的测试证明。

范围、方程、原始资料及下一步见[方法草稿](../../SHOCK_LAYER_ENTRY_V97_DRAFT.md)。本轮未保存全套矩阵或工程ZIP；源论文与读页图仅在本地研究目录，源码、必要小输入和小型证据随GitHub提交同步。
