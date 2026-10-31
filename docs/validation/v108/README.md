# v108验证资料与复现边界

本目录对应[版本说明](../../SU2_INITIALIZATION_MESH_V108.md)。`PLAN.md`在数值运行前确定预算和筛选尺度；`delivery_manifest.json`绑定实际终态、准确源码及证据SHA。源码、小型来源和局部结果随GitHub同步，完整大型输出仍在本地。

- `processes/`：4项MATLAB终态及原日志，3成功/1失败。`failures.json`说明初场配置格式拒绝和含NaN结构的测试误报，没有把补写的说明伪装成原stderr。
- `engines/warm/`：实际收敛重启动的配置、回执、壁面VTK、日志和**仅末两行**的历史摘录；完整23478行历史的SHA在manifest，摘录不是全历史。
- `engines/normal/`：实际58560单元两步运行，同样保持未收敛身份；与原输入方向完全匹配。
- `checks/`：独立初场/平板数值比较、几何审计、原面映射及UI定向验证；`independent_mapping_audit.json`独立从VTK复算三种原面场与两种面积负荷。
- `initialization_seed/`：762111字节的内部粗解CSV与其原始迭代历史的无损gzip（解压后SHA必须与manifest一致）。这是小型数值来源，不是整个工程包。其余粗解配置、原生网格、回执、壁面和日志复用v105已交付文件，路径/SHA见manifest的`sharedCoarseInputs`。
- `source_versions/`：按准确字节保留各阶段源码，包括失败版本。不是新工程副本；不要把当前源码套到旧日志。
- `upstream/`：SU2固定8.6.0格式/方向依据，附官方链接和原版权/许可证。

本轮原始结果在`outputs/v108_validation/`。`plate_warm01/restart_flow.csv`是已收敛细解，`initial_flow.csv`是其真实初场，不能混用。`normal01/normal_volume.mat`和`normal01/native_input/`是后续三维工作起点；两步重启动只可作研究检查点。它们未全部上传，不得因GitHub同步而删除。

无须重算流场即可复核主UI网格行为：在保留本地几何结果的工作区，MATLAB执行（输出目录须尚未存在）：

```matlab
addpath('tests');
testSU2MeshModeUI('outputs/v108_ui_recheck', ...
    'outputs/v108_validation/normal01/normal_volume.mat');
```

从干净检出复现几何时，先调用`studySU2NormalLayers('outputs/v108_geometry_reproduction')`；它生成完整原生几何与配置而不启动求解器，然后以该目录的`normal_volume.mat`作为上述UI测试输入。`auditSU2NormalLayers.py <目录>`独立核对原生方向和首层距离。历史UI/取消试验不需再次运行。

细网格初始化复现：在新`outputs`目录建立粗解来源子目录，按manifest的共享输入表复制v105原字节，将`coarse_restart_flow.csv`命名为`restart_flow.csv`，用Python `gzip.decompress`恢复`coarse_history.csv.gz`为`history.csv`并核对SHA。MATLAB生成新细网格：

```matlab
buildSU2FlatPlateCase('outputs/v108_reproduction/fine', ...
    struct('nx',120,'ny',96,'firstCellHeight',1e-6,'iterations',30000));
```

随后运行`python scripts/prepareSU2FlatPlateWarmStart.py <新粗解来源目录> outputs/v108_reproduction/fine`及`python tests/auditSU2WarmStart.py outputs/v108_reproduction/fine`。这不启动SU2；来源配置与物性必须一致。实际重启动另用`runEmbeddedSU2Case`明确2线程/720秒预算，保持输出目录唯一；研究计划未变时不重复跑已通过的算例。数值算法可复现不意味着不同机器墙钟时间或文本每一末位都必然相同。

`auditSU2WarmConvergence.py`、`auditSU2V108.py`与证据汇总脚本用于核对本次原路径和已留存输入；不是任意目录的通用复现入口。修改路径做新研究须保留新的来源绑定，不覆盖本次历史证据。
