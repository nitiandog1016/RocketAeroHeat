# 当前交接

2026-11-21。民用光滑 B 型热敏珠已通过独立安装的单工况全表面数值复算。在 6 m/s、距构造初场 7.809360865 µs 时，全周 56 层网覆盖 **5,216 / 5,216 面、3.08983026 mm²**，整珠放热 **10.69725502 mW**。与 v295 相比新增 **0 面、0 m²**。本版从程序内部恢复初场，在另一安装中重新推进，逐面结果与参考表逐字节一致。按论文面积换算的启动时刻平均换热系数比实验读数高 18.28%，局部壁压与热流的联合实验验证面积为 **0 m²**。[v296 报告](validation/v296/RESULTS.md) · [原面和实测量对照](validation/v296/evidence/installation_result.json)

<!-- records:last_completed=v296 -->
<!-- records:active=none -->
<!-- records:base_commit=c7f0bcd6d0ec2ab120bb7fa1dc4a6d13c70e2036 -->

## 热敏珠当前结果

从“选用 B 型热敏珠”进入网格页，选择“6 m/s · 360° 全周 · 56 层”，准备体网、确认边界后可新算、读取已完成计算或接续。第二套安装的任务 `d457118c-1e76-4618-ab24-5bcab5bfce09` 已完成全部 11 阶段，MATLAB 耗时 45.63 min。末步收支差 −0.9432 nW，热峰位于工作台面元 5211。初场冷态恢复与新算在同机两个安装目录验收。另一台电脑复算需安装固定依赖。[本机三维页](../outputs/v296_validation/desktop_final/outputs/thermistor_view_fb14f0d9-dab2-4f6f-84e0-e505ac1861cb/workbench/index.html) · [离线逐面核对](validation/v296/evidence/offline.json) · [安装方法](validation/v296/INSTALL.md)

v295 曾在一次独立壁反力读取中断，原读数终态未知，现场留存。同一工作台已验证从成功阶段接续这一任务，再读完整两半步。[中断与接续](validation/v295/RESULTS.md)

在主侧栏的 B 型热敏珠菜单选择“查看 6 m/s · 近壁网格对照”。主值为 56 层，可核对两网同刻压力、热流及相同位置的温度。实际截面 X = 0.172 mm 导出 121 行，逐行对应同一批原面。6 m/s 运行、热量核对及界面验收均有终态回执。[界面回执](validation/v294/evidence/ui.json) · [本机三维页](../outputs/v294_validation/view01/export/workbench/index.html)

v292 在 6 m/s、293.15 K 空气、给定壁温 303.15 K 下，以一个整步及两个半步到达相同时刻。48 层网的局部最大步长差分别为 0.002772 Pa、0.006122 W/m²，放热峰均为工作台面元 5211。发展后尾迹和空间收敛阶待验证。[同刻差异](validation/v292/evidence/comparison.json)

两个半步各在第 260 次内迭代保存终点，末尾分别有 38、72 次六项残差连续达标。初次回放把不同更新阶段的压力与标量读到一起。改用固定端点装配后，两步全域压力与接受场的无量纲最大差均低于 1.35 × 10⁻¹⁵。第二步放热 10.726572344 mW，收支差 −2.498 pW。主侧栏“查看 6 m/s · 同刻时间步对照”支持同原面差场、统一色标、探针和截面导出。[首步热量账](validation/v292/evidence/first_energy.json) · [第二步热量账](validation/v292/evidence/second_energy.json) · [本机对照三维页](../outputs/v292_validation/view02/export/workbench/index.html)

6 m/s 全周快照以原生强温差壁反力保守映回原面。绝对壁压 101309.77–101351.96 Pa，入壁热流 −6606.75 至 −1152.12 W/m²，首层流体温度 302.665–303.068 K。工作台热峰面元为 5211。整珠放热 10.726572556 mW，远场流出 10.726560354 mW，空气储热率 12.192505 nW，收支差 9.764 pW。在主侧栏的 B 型热敏珠菜单选择“查看 6 m/s · 全周启动瞬态”。MATLAB 工作台、离线三维页、探针和下载保留相同物理时刻与原面值。[同阶段热量账](validation/v290/evidence/surface.json) · [界面核对](validation/v290/evidence/ui.json) · [本机三维页](../outputs/v290_validation/ui03/export/workbench/index.html)

3 m/s、293.15 K 空气和 4 mW 净功率下，15° 与 7.5° 网分别求得 298.26613 K 与 298.19793 K 的均匀集总壁温。全部 5,216 面均有压力、热流、首层流体温度与壁温，放热积分为 4 mW。[3 m/s 结果](validation/v286/RESULTS.md) · [本机细网三维页](../outputs/v286_validation/fine01/workbench/index.html)

主侧栏可选择给定壁温或净加热功率。净功率按扣除引线、辐射等损失后的对流份额输入。零功率时各面热流为零，壁温为 293.15 K。各自网格的支持上限为 9.38210、9.50885 mW。15° 网用独立半温差原生功率反求 299.15 K，误差为 5.72 × 10⁻⁷ K。7.5° 网面对介于两上限之间的功率时可提供主结果，单列 15° 同功率对照超域。[15° 验收](validation/v286/evidence/grid15.json) · [细网验收](validation/v286/evidence/grid7p5.json) · [直接原生对照](validation/v285/evidence/native.json)

给定壁温支持 293.15–305.15 K。在 305.15 K 参照点，两网平均 h 相对图 7 读数 218.36 W/(m²·K) 的偏差为 +15.73% 和 +17.29%。真实珠体含引线；名义模型采用光滑旋转长椭球、沿长轴来流与等温假设。逐面压力和热流的联合独立验证面积为 0 m²。给定壁温模式保留原有变量身份；功率模式提供集总热平衡壁温，内部固体温度梯度未知。历史弱焓超温、二维重建差额及给定温差直接计算均留档。[图 7 读数](validation/v282/evidence/comparison.json) · [细网独立热量账](validation/v283/evidence/energy.json) · [给定壁温与物理零](validation/v285/RESULTS.md)

## 下一项

[当前下一决定](validation/v296/NEXT_DECISION.md) 是在 3 m/s 基准工况下加入论文给出线径的引线，分别验收珠体、引线和连接区的原面、壁压、热流及功率。引线露出长度、接入位置、端面和热边界需在开算前声明。当前结果属于光滑珠体。[冷态恢复](validation/v296/evidence/cold.json) · [空缓存网格](validation/v296/evidence/empty_mesh.json)

6 m/s 名义横径雷诺数 360.20，沿轴长椭球首次对称破缺的临界图读约 241。30° 镜面扇区只提供基流初猜。原论文图 7 的 6 m/s 测点测得整珠平均 h 约 292.32 W/(m²·K)，温差读数约 10.09 K。[图读与适用范围](validation/v287/evidence/measurements.json)

## 自绘模型保留的状态

自绘模型有 19,752 个原面，总面积 1.9646292565715553 m²。冻结工程结果中 9,725 面同时有壁压和热流，共 **0.846754 m²（43.10%）**，联合缺口 1.117875 m²，全表面热功率未知，联合独立验证面积为 0 m²。936 个已有工程壁压的原面共 0.0957078182 m²，三种存档直接补值路线已排除。v263 的侧壁工作保留未完成状态。[原面覆盖](validation/v263/RESULTS.md) · [上次核对](validation/v281/RESULTS.md) · [原在研计划](validation/v263/PLAN.md)

保留原近壁 32 层和 2 µm 首层、历史 19 层及尾缘资产，原物理门槛保持。v252 旧三链父监督终态仍未知。v280 风洞圆柱的两组首个物理步均未通过残差门槛。历史失败和原始回执可通过 [执行索引](EXECUTION_INDEX.md) 查阅。29 个保护文件、54 项工程文件和两对候选均已核对。大型场继续保存在 E 盘，发布仅包含本任务的小型数据、证据和代码。
