# v141 验证索引

连续包络小例、真实壳几何和CGAL无自交通过；唯一完整外域有31个低质量四面体，拒绝。无SU2/MATLAB/新UI测试，新资格0。

- `plan.md`：运行前固定层高规则、保护范围、解析位移界及预算。
- `fixture/`与`fixture_process/`：变形八面体解析体积、源锁及通过结果。
- `actual/`与`actual_process/`：真实准备、CGAL与Gmsh、完整几何/质量拒绝和准确日志。
- `review_process/`：独立层高/坐标、原生连接/边界质量复核。
- `tetgen_research/`：官方手册与固定版本元数据，只是下一路线依据，尚未集成。
- `source_versions/`与`storage_and_sources.json`：准确执行源码、输入/输出身份和本轮文件长度。

大网格在本地`outputs/v141_validation/actual/`，尚未取得网格接受或物理资格。总体结论和下一起点见[阶段说明](../../CONTINUOUS_PRISM_ENVELOPE_V141.md)。
