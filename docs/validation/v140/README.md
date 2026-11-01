# v140 验证索引

小例通过，真实渐变终止壳因138个低质量四面体拒绝。无新完整外域、CGAL自交、MATLAB、SU2或UI验证，新资格0。

- `plan.md`：运行前规则、预算和停止门槛。
- `fixture/`与`fixture_process/`：已知体积5的小例输入/三类单元及通过结果。
- `actual/`与`actual_process/`：真实准备、壳体积分、源锁及质量拒绝原文。
- `review_process/`：独立体积/逆矩阵质量、原列身份与失败层位置的只读复核。
- `source_versions/`：准确执行源；`storage_and_sources.json`：本轮哈希、文件长度及未执行分支说明。

大网格保留于本地`outputs/v140_validation/actual/prepared_transition.mat`，不是完整外域网格。不要将本轮shell体积或节点数解释为全流场域结果，也不要绕过质量门槛继续启动外域或流场。详见[阶段说明](../../GRADED_LAYER_TERMINATION_V140.md)。
