# 当前交接

## 简化运载火箭飞行工况与低速风洞

风洞 9,216 面已算。M0.30 的 6.25 µs 研究场入壁功率 −1,045.057 W。原生边界与储能之和 −13.53 kW；扣除参考焓的质量输运后为 +140.18 W。下一项在同刻检验最外两层加密。五点最终功率待验。[原生守恒账](validation/v370/RESULTS.md) · [下一项](validation/v370/NEXT_DECISION.md)

<!-- records:last_completed=v370 -->
<!-- records:active=none -->
<!-- records:base_commit=31fd2d9d7cfb89be223923c39178333d6d717c18 -->

## RFZ 独立计算与公开对照

deg190 目录标称 Mach 5.35，工作簿列 Mach 5.38。论文表 4 编号 3 对应工作簿第 4 行 Case 2。72,706 个作者样点可连同来源导出。壁压零点、热流正向待核实。v348 实跑核对了 130,044 面、740.729003 m² 和 300 K 壁温，入口随后以 `RAH:RFZ:SurfaceInput` 退出。独立壁压与热流同面覆盖为 0 面、0 m²。第 13 区有 2,409 面的壁邻四面体四点同温，其余 20,443 面的内部对顶点温度待求。近壁分辨率及外域精度待验。[运行回执](validation/rfz_entry_20261126/receipt.json) · [标准工况](../cases/rfz_rec/full_cad_input.json) · [节点温界](validation/v346/RESULTS.md) · [v348 对照条件](validation/v348/NEXT_DECISION.md)

## 公开试件的既有结果

公开实验试件 Short Test 65 单独核对了 32,480 / 32,480 面计算壁压与热流。球面量热输入 0.827123193 W，球温比实验均值低 11.157 K，七站 Nu 的均方根误差为 4.8882。[试件独立验收](validation/v345/TEST65.md)

20D 的两半步放热为 2.005187691 W，较 14D 减少 0.01221%。[远场对照](validation/v335/RESULTS.md)

双支管 128 扇区的稳态场是本次扰动基线。法向和周向加密对照均已完成。119.8° 测站顺向、166.7° 反向。[稳态与周向对照](validation/v332/RESULTS.md) · [近壁诊断](validation/v331/RESULTS.md) · [法向对照](validation/v330/RESULTS.md)

独立光滑球对照在扰动后 82.10 ms 已有 4,464 / 4,464 面同时刻壁压与热流，细步总放热 0.977093924 W，比实验折算高 16.68%。该段尚无周期统计，联合实验验证面积 0 m²。[瞬时结果](validation/v322/RESULTS.md) · [本机三维页](../outputs/v322_validation/read01/wakeFine/workbench/index.html) · [测量定义](validation/v321/SOURCE.md)

B 珠已有 10,276 面、3.196650633 mm² 全表面结果。第三档珠部换热比图读高 22.80%，两轮网格排列候选已关闭。带线模型 4 mW 集总壁温为 297.694199 K，热峰在迎风引线 207 号面。[重剖子面验收](validation/v320/RESULTS.md) · [网格试验](validation/v307/RESULTS.md) · [分带热量账](validation/v305/evidence/workbench.json)

## 带线试件

v298 配对中，带线珠体含接头放热 11.94750 mW，相对同方法光滑珠体减少 0.36120 mW。原带线试件覆盖 10,272 / 10,272 面，6,144 个共同珠面可逐面比较。[带线三维页](../outputs/v298_validation/wired_ui03/workbench/index.html) · [光滑三维页](../outputs/v298_validation/smooth_ui03/workbench/index.html) · [三项配对](validation/v298/evidence/comparison.json)

从“选用带引线 B 型热敏珠”进入五步侧栏。工作台已通过网格、边界、冷态新算、归档、三维云图与七区探针验收。最大放热在迎风侧引线面元 207。后引线、后接头和珠体背流侧共有 1,536 面出现壁剪切与首层流速同时反向，两套近壁网辨认出同一批原面。新算任务 `thermistor_leads_run_d5d2003c-4d02-45b2-8fd1-0dd2641a07e0` 已有正常终态和逐面核对。引线长度、接头圆角、全体等温均为名义假设。[本机三维页](../outputs/v297_validation/ui03/workbench/index.html) · [科研出口](validation/v297/evidence/workbench.json) · [近壁数据](validation/v297/evidence/nearwall.json)

## 光滑热敏珠

v299 光滑控制完成了 6,624 原面的中段加密对照，放热由 12.30870 降至 10.62365 mW；中段与原结构网的额外放热缩小 71.98%。[原面分带](validation/v299/evidence/comparison.json) · [三维页](../outputs/v299_validation/midbody_ui01/workbench/index.html)

两套安装已验收光滑珠体的读档和新算；另一台电脑复算需要固定依赖。[离线逐面核对](validation/v296/evidence/offline.json) · [安装方法](validation/v296/INSTALL.md)

v295 曾在一次独立壁反力读取中断，原读数终态未知，现场留存。同一工作台已验证从成功阶段接续这一任务，再读完整两半步。[中断与接续](validation/v295/RESULTS.md)

侧栏“查看 6 m/s · 近壁网格对照”用 56 层为主值，X = 0.172 mm 截面已核对全部 121 行。[界面回执](validation/v294/evidence/ui.json) · [本机三维页](../outputs/v294_validation/view01/export/workbench/index.html)

v292 在 6 m/s、293.15 K 空气、给定壁温 303.15 K 下，以一个整步及两个半步到达相同时刻。48 层网的局部最大步长差分别为 0.002772 Pa、0.006122 W/m²，放热峰均为工作台面元 5211。发展后尾迹和空间收敛阶待验证。[同刻差异](validation/v292/evidence/comparison.json)

v292 更正跨阶段读取后，48 层网第二半步放热 10.726572344 mW，收支差 −2.498 pW。[原步热量账](validation/v290/evidence/surface.json) · [半步热量账](validation/v292/evidence/second_energy.json)

3 m/s、293.15 K 空气和 4 mW 净功率下，15° 与 7.5° 网分别求得 298.26613 K 与 298.19793 K 的均匀集总壁温。全部 5,216 面均有压力、热流、首层流体温度与壁温，放热积分为 4 mW。[3 m/s 结果](validation/v286/RESULTS.md) · [本机细网三维页](../outputs/v286_validation/fine01/workbench/index.html)

主侧栏可选择给定壁温或净加热功率。净功率按扣除引线、辐射等损失后的对流份额输入。零功率时各面热流为零，壁温为 293.15 K。各自网格的支持上限为 9.38210、9.50885 mW。15° 网用独立半温差原生功率反求 299.15 K，误差为 5.72 × 10⁻⁷ K。7.5° 网面对介于两上限之间的功率时可提供主结果，单列 15° 同功率对照超域。[15° 验收](validation/v286/evidence/grid15.json) · [细网验收](validation/v286/evidence/grid7p5.json) · [直接原生对照](validation/v285/evidence/native.json)

给定壁温支持 293.15–305.15 K。在 305.15 K 参照点，两网平均 h 相对图 7 读数 218.36 W/(m²·K) 的偏差为 +15.73% 和 +17.29%。真实珠体含引线；名义模型采用光滑旋转长椭球、沿长轴来流与等温假设。逐面压力和热流的联合独立验证面积为 0 m²。给定壁温模式保留原有变量身份；功率模式提供集总热平衡壁温，内部固体温度梯度未知。历史弱焓超温、二维重建差额及给定温差直接计算均留档。[图 7 读数](validation/v282/evidence/comparison.json) · [细网独立热量账](validation/v283/evidence/energy.json) · [给定壁温与物理零](validation/v285/RESULTS.md)

## 自绘模型保留的状态

自绘模型有 19,752 个原面，总面积 1.9646292565715553 m²。冻结工程结果中 9,725 面同时有壁压和热流，共 **0.846754 m²（43.10%）**，联合缺口 1.117875 m²，全表面热功率未知，联合独立验证面积为 0 m²。936 个已有工程壁压的原面共 0.0957078182 m²，三种存档直接补值路线已排除。v263 的侧壁工作保留未完成状态。[原面覆盖](validation/v263/RESULTS.md) · [上次核对](validation/v281/RESULTS.md) · [原在研计划](validation/v263/PLAN.md)

保留原近壁 32 层和 2 µm 首层、历史 19 层及尾缘资产，原物理门槛保持。v252 旧三链父监督终态仍未知。v280 风洞圆柱的两组首个物理步均未通过残差门槛。历史失败和原始回执可通过 [执行索引](EXECUTION_INDEX.md) 查阅。29 个保护文件、54 项工程文件和两对候选均已核对。大型场继续保存在 E 盘，发布仅包含本任务的小型数据、证据和代码。
