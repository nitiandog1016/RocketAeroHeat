# 当前交接

2026-11-23。Short Test 65 名义双支管银球（Re = 882）两套网均有 **32,480 / 32,480 面、0.000700423853595 m²** 的稳态壁压与热流，联合实验验证面积 0 m²。v331 已核对近壁流向，新增壁压热流 0 面、0 m²。两网的 119.8° 测站均顺向、166.7° 均反向，沿 +Z 测线的变号区间为 143.4375°–149.0625°。反向面积五层 42.118、十层 45.110 mm²，热峰面 25,055 均顺向。细网外露银球折合量热比实验高 38.61%，七站 RMSE 为 4.659 Nu。[近壁诊断](validation/v331/RESULTS.md) · [本机流向云图](../outputs/v331_validation/ui02/workbench/index.html)

<!-- records:last_completed=v331 -->
<!-- records:active=none -->
<!-- records:base_commit=923acc184349702c31da1ff0da548ba1014d18ad -->

五区共享 6,016 个源壁片。法向首层从 25.4 减至 12.7 µm，整件放热由 2.011094 变为 2.009415 W。两网近峰区交并面积比 99.08%。三维流向入口按所选网格载入，七站可并列核对。[法向网格对照](validation/v330/RESULTS.md) · [独立冷算](validation/v329/RESULTS.md)

独立光滑球对照在扰动后 82.10 ms 已有 4,464 / 4,464 面同时刻壁压与热流，细步总放热 0.977093924 W，比实验折算高 16.68%。该段尚无周期统计，联合实验验证面积 0 m²。[瞬时结果](validation/v322/RESULTS.md) · [本机三维页](../outputs/v322_validation/read01/wakeFine/workbench/index.html) · [测量定义](validation/v321/SOURCE.md)

B 珠已有 10,276 面、3.196650633 mm² 全表面结果。第三档珠部换热比图读高 22.80%，两轮网格排列候选已关闭。带线模型 4 mW 集总壁温为 297.694199 K，热峰在迎风引线 207 号面。[重剖子面验收](validation/v320/RESULTS.md) · [网格试验](validation/v307/RESULTS.md) · [分带热量账](validation/v305/evidence/workbench.json)

## 带线试件

v298 配对中，带线珠体含接头放热 11.94750 mW，相对同方法光滑珠体减少 0.36120 mW。原带线试件覆盖 10,272 / 10,272 面，6,144 个共同珠面可逐面比较。[带线三维页](../outputs/v298_validation/wired_ui03/workbench/index.html) · [光滑三维页](../outputs/v298_validation/smooth_ui03/workbench/index.html) · [三项配对](validation/v298/evidence/comparison.json)

从“选用带引线 B 型热敏珠”进入五步侧栏。工作台已通过网格、边界、冷态新算、归档、三维云图与七区探针验收。最大放热在迎风侧引线面元 207。后引线、后接头和珠体背流侧共有 1,536 面出现壁剪切与首层流速同时反向，两套近壁网辨认出同一批原面。新算任务 `thermistor_leads_run_d5d2003c-4d02-45b2-8fd1-0dd2641a07e0` 已有正常终态和逐面核对。引线长度、接头圆角、全体等温均为名义假设。[本机三维页](../outputs/v297_validation/ui03/workbench/index.html) · [科研出口](validation/v297/evidence/workbench.json) · [近壁数据](validation/v297/evidence/nearwall.json)

## 光滑热敏珠

v299 光滑控制完成了 6,624 原面的中段加密对照，放热由 12.30870 降至 10.62365 mW；中段与原结构网的额外放热缩小 71.98%。[原面分带](validation/v299/evidence/comparison.json) · [三维页](../outputs/v299_validation/midbody_ui01/workbench/index.html)

从“选用 B 型热敏珠”进入网格页，选择“6 m/s · 360° 全周 · 56 层”，准备体网、确认边界后可新算、读取已完成计算或接续。第二套安装的任务 `d457118c-1e76-4618-ab24-5bcab5bfce09` 已完成全部 11 阶段，MATLAB 耗时 45.63 min。末步收支差 −0.9432 nW，热峰位于工作台面元 5211。初场冷态恢复与新算在同机两个安装目录验收。另一台电脑复算需安装固定依赖。[本机三维页](../outputs/v296_validation/desktop_final/outputs/thermistor_view_fb14f0d9-dab2-4f6f-84e0-e505ac1861cb/workbench/index.html) · [离线逐面核对](validation/v296/evidence/offline.json) · [安装方法](validation/v296/INSTALL.md)

v295 曾在一次独立壁反力读取中断，原读数终态未知，现场留存。同一工作台已验证从成功阶段接续这一任务，再读完整两半步。[中断与接续](validation/v295/RESULTS.md)

侧栏“查看 6 m/s · 近壁网格对照”用 56 层为主值，X = 0.172 mm 截面已核对全部 121 行。[界面回执](validation/v294/evidence/ui.json) · [本机三维页](../outputs/v294_validation/view01/export/workbench/index.html)

v292 在 6 m/s、293.15 K 空气、给定壁温 303.15 K 下，以一个整步及两个半步到达相同时刻。48 层网的局部最大步长差分别为 0.002772 Pa、0.006122 W/m²，放热峰均为工作台面元 5211。发展后尾迹和空间收敛阶待验证。[同刻差异](validation/v292/evidence/comparison.json)

v292 更正跨阶段读取后，48 层网第二半步放热 10.726572344 mW，收支差 −2.498 pW。[原步热量账](validation/v290/evidence/surface.json) · [半步热量账](validation/v292/evidence/second_energy.json)

3 m/s、293.15 K 空气和 4 mW 净功率下，15° 与 7.5° 网分别求得 298.26613 K 与 298.19793 K 的均匀集总壁温。全部 5,216 面均有压力、热流、首层流体温度与壁温，放热积分为 4 mW。[3 m/s 结果](validation/v286/RESULTS.md) · [本机细网三维页](../outputs/v286_validation/fine01/workbench/index.html)

主侧栏可选择给定壁温或净加热功率。净功率按扣除引线、辐射等损失后的对流份额输入。零功率时各面热流为零，壁温为 293.15 K。各自网格的支持上限为 9.38210、9.50885 mW。15° 网用独立半温差原生功率反求 299.15 K，误差为 5.72 × 10⁻⁷ K。7.5° 网面对介于两上限之间的功率时可提供主结果，单列 15° 同功率对照超域。[15° 验收](validation/v286/evidence/grid15.json) · [细网验收](validation/v286/evidence/grid7p5.json) · [直接原生对照](validation/v285/evidence/native.json)

给定壁温支持 293.15–305.15 K。在 305.15 K 参照点，两网平均 h 相对图 7 读数 218.36 W/(m²·K) 的偏差为 +15.73% 和 +17.29%。真实珠体含引线；名义模型采用光滑旋转长椭球、沿长轴来流与等温假设。逐面压力和热流的联合独立验证面积为 0 m²。给定壁温模式保留原有变量身份；功率模式提供集总热平衡壁温，内部固体温度梯度未知。历史弱焓超温、二维重建差额及给定温差直接计算均留档。[图 7 读数](validation/v282/evidence/comparison.json) · [细网独立热量账](validation/v283/evidence/energy.json) · [给定壁温与物理零](validation/v285/RESULTS.md)

## 下一项

在十层网基础上进行一次壁面切向分辨率对照，检查两处背风测站的换热及变号区间。[下一决定](validation/v331/NEXT_DECISION.md)

## 自绘模型保留的状态

自绘模型有 19,752 个原面，总面积 1.9646292565715553 m²。冻结工程结果中 9,725 面同时有壁压和热流，共 **0.846754 m²（43.10%）**，联合缺口 1.117875 m²，全表面热功率未知，联合独立验证面积为 0 m²。936 个已有工程壁压的原面共 0.0957078182 m²，三种存档直接补值路线已排除。v263 的侧壁工作保留未完成状态。[原面覆盖](validation/v263/RESULTS.md) · [上次核对](validation/v281/RESULTS.md) · [原在研计划](validation/v263/PLAN.md)

保留原近壁 32 层和 2 µm 首层、历史 19 层及尾缘资产，原物理门槛保持。v252 旧三链父监督终态仍未知。v280 风洞圆柱的两组首个物理步均未通过残差门槛。历史失败和原始回执可通过 [执行索引](EXECUTION_INDEX.md) 查阅。29 个保护文件、54 项工程文件和两对候选均已核对。大型场继续保存在 E 盘，发布仅包含本任务的小型数据、证据和代码。
