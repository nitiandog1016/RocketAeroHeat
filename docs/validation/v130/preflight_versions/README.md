# v130启动前失败与准确版本

这些记录不是失败的SU2计算。最终只有control与column_final获得原生启动；其他准备保持未运行身份。

1. `prepare_unbounded.py.raw`产生column：37个正k低于SU2原下限，未启动。
2. `prepare_global_column_norm.py.raw`产生column_bounded：原整列误差通过；`verify_initial_strict.py.raw`独立逐节点检查失败，omega最大5.827785019638249e−7、k最大6.956172542712249e−12，检查线1e−10。旧场保留。
3. `prepare_local_global_linesearch.py.raw`加严局部停止条件，但仍使用全局绝对残差下降的线搜索；column_local准备在Newton下降检查中失败，未写初态，仅生成配置/网格；该目录NOT_RUN记录身份。
4. `prepare_local_merit.py.raw`用局部相对范数下降生成column_final。真实全部列准备成功，但两端k均取下限的独立合成测试触发线搜索停滞，原输出见`local_merit_unit_failure.log`。最终增加当前迭代固定权重的备用下降范数，仍保持同一严格停止条件和80次上限。4项测试通过；最终运行准备器逐字段重建确认已生成column_final完全一致，故无需重写原场。

最终独立有量纲斜率核对见`evidence/initial_state_independent.json`：omega局部误差3.154701e−13、k互补缺陷8.961101e−12；全部端点、未选节点、其他字段及热内能保持。最终目录名称为column_final，LOCAL_RESIDUAL事前文件中的column_local是先前拟用名称，文件保留原字节以维持绑定。

一次临时诊断包装通过stdin执行时未设置`__file__`，在解析根目录时失败；明确传入文件路径后只读打印上述缺陷，没有产生新原生任务或改变数值。原始失败堆栈当时由终端返回；未保存的终端输出不冒称完整原始日志。这里仅保存实际已有源版本、日志和目录标记。
