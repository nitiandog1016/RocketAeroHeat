# 执行历史索引

这里用于检索历史证据，不确定当前执行顺序。当前状态唯一入口是 [CURRENT_STATE.md](CURRENT_STATE.md)。历史报告的下一步不覆盖当前页。

## 当前与独立工作入口

- [v329 双支管银球单工况全表面结果](validation/v329/RESULTS.md) · [热流云图与七站](validation/v329/evidence/surface.png) · [下一项](validation/v329/NEXT_DECISION.md)。统一远场通量后 32,480 面 p/q 全覆盖；两次独立冷算面表一致，外露银球折合量热高 38.43%。

- [v328 双支管银球的有界稳态及远场热量差](validation/v328/RESULTS.md) · [逐次历史、边界与原始回执](validation/v328/evidence/receipts.zip) · [下一项](validation/v328/NEXT_DECISION.md)。600 次更新有界且残差收敛，远场与壁面差 0.640012%。

- [v327 双支管银球扩散精度检验](validation/v327/RESULTS.md) · [三档网格结果](validation/v327/evidence/cuts.csv) · [下一项](validation/v327/NEXT_DECISION.md)。沿边幅值限制在斜网解析场上完全丢失横向截面通量，该候选已关闭。

- [v326 双支管银球的标量通量对照](validation/v326/RESULTS.md) · [首次越界诊断图](validation/v326/evidence/flux.png) · [下一项](validation/v326/NEXT_DECISION.md)。质量通量对照仍越界，首次最低值节点已核实八条边和隐式更新。

- [v325 双支管流体体网与原生启动](validation/v325/RESULTS.md) · [网格与越界位置](validation/v325/evidence/mesh.png) · [下一项](validation/v325/NEXT_DECISION.md)。6,016 个壁片完成五区守恒映射；两套外区的首段均因主标量越界未获接受。

- [v324 双支管银球几何与逐原面外露面积](validation/v324/RESULTS.md) · [接头图](validation/v324/evidence/geometry.png) · [下一项](validation/v324/NEXT_DECISION.md)。88 个原球面受支管影响，两处接缝连续，七站均位于外露球面。

- [v323 银球球后回流与六截面空间分辨率](validation/v323/RESULTS.md) · [原节点速度剖面](validation/v323/evidence/profiles.png) · [下一项](validation/v323/NEXT_DECISION.md)。轴线回流长 1.474 mm，两处背风测站的首层流向相反。

- [v322 Short 银球 82.10 ms 全表面瞬时解](validation/v322/RESULTS.md) · [尾迹与背风换热](validation/v322/evidence/wake.png) · [下一项](validation/v322/NEXT_DECISION.md)。两种步长、六个同刻、4,464 面的结果通过数值核对，背风两站与实验的趋势差异仍在。

- [v321 Short 银球全周基流与七站对照](validation/v321/RESULTS.md) · [局部 Nu 曲线](validation/v321/evidence/stations.png) · [下一项](validation/v321/NEXT_DECISION.md)。4,464 面均有基流压力和热流，119.8° 处高估 139.3%。

- [v320 带线 B 珠的局部重剖与细壁积分](validation/v320/RESULTS.md) · [逐子面数值](validation/v320/evidence/verification.json) · [下一项](validation/v320/NEXT_DECISION.md)。两个邻区共六个子面，五条路径核对完整的 10,276 个导入面。

- [v319 普通 STL 的 B 珠全表面新算与导出](validation/v319/RESULTS.md) · [示例 STL](../cases/thermistor_b/probe-mm.stl) · [下一项](validation/v319/NEXT_DECISION.md)。3 m/s 下按导入面号交付全部 10,272 面，长引线新算放热 11.14834 mW。

- [v318 电阻等效温度与 B 珠表面温度](validation/v318/RESULTS.md) · [测温与功率对照](validation/v318/evidence/thermometry.png) · [下一项](validation/v318/NEXT_DECISION.md)。主测点内外温差 0.11058 K，珠部功率减少 0.09442 mW，剩余偏差 +25.22%。

- [v317 B 珠长引线热边界与同步均温解](validation/v317/RESULTS.md) · [沿线温度与分段放热](validation/v317/evidence/wire.png) · [下一项](validation/v317/NEXT_DECISION.md)。引线降温使整件放热减少，珠部放热增加 0.01800 mW。

- [v316 带线 B 珠独立冷启动](validation/v316/RESULTS.md) · [超差原面](validation/v316/evidence/outside_tolerance.csv) · [下一项](validation/v316/NEXT_DECISION.md)。全表面 10,272 面均已求解。与 v315 的局部差异集中在两根引线和端面。

- [v315 B 珠 3 m/s、12 K 风洞测点](validation/v315/RESULTS.md) · [平均换热对照图](validation/v315/evidence/mean_heat.png) · [下一项](validation/v315/NEXT_DECISION.md)。96 档名义侧吹珠部平均 h 比原图读数高 26.16%，换算与直接热边界逐面等步核对通过。

- [v314 带线热敏珠 48／96 档周向对照](validation/v314/RESULTS.md) · [迎背风曲线](validation/v314/evidence/section.png) · [下一项](validation/v314/NEXT_DECISION.md)。整件放热下降 0.352%，同截面背风近邻放热密度增加 12.77%。

- [v313 带线热敏珠侧吹前端面近壁细化](validation/v313/RESULTS.md) · [端面放热图](validation/v313/evidence/cap_map.png) · [下一项](validation/v313/NEXT_DECISION.md)。端面首层距离减半，最强放热密度下降 0.00228%，两共边镜像面交换排序。

- [v312 带线 B 珠侧吹远场半径对照](validation/v312/RESULTS.md) · [周向差值图](validation/v312/evidence/section_domains.png) · [下一项](validation/v312/NEXT_DECISION.md)。20D 扩至 30D 后，整件放热减少 0.00473%，两域热峰均在第 10267 面。

- [v311 带线 B 珠全周侧吹](validation/v311/RESULTS.md) · [实际周向分布](validation/v311/evidence/circumference.png) · [下一项](validation/v311/NEXT_DECISION.md)。3 m/s 沿 −Y 侧吹的 10,272 面均已计算，整件向空气放热 5.74320 mW。

- [v310 带线热敏珠珠部平均强迫对流预测](validation/v310/RESULTS.md) · [独立复算](validation/v310/evidence/verification.json) · [下一项](validation/v310/NEXT_DECISION.md)。3 m/s、299.15 K 工况已可在侧栏比较式（28）和三维珠部功率，并保留表 1 极限 Nu 的差异。

- [v309 B 型热敏珠 6 K 实验测点核对](validation/v309/RESULTS.md) · [工作台与逐面核对](validation/v309/evidence/verification.json) · [下一项](validation/v309/NEXT_DECISION.md)。两幅原图的相邻测点均与 6 K 分开，计算系数与未知实验偏差已接回同一工况。

- [v308 带线热敏珠可输入均匀壁温](validation/v308/RESULTS.md) · [逐面数值核对](validation/v308/evidence/verification.json) · [下一项](validation/v308/NEXT_DECISION.md)。两档网在同一 299.15 K 工况下导出全表面结果，第二档逐面比较独立冷态。

- [v307 带线热敏珠前中段固定节点换边](validation/v307/RESULTS.md) · [成品几何](validation/v307/evidence/geometry.json) · [下一项](validation/v307/NEXT_DECISION.md)。169 次换边保住了原 16 条测线，22 条交错测线的最大间距超限，4 个近退化棱柱未通过面锥检查。

- [v306 热敏珠前中段切向重排](validation/v306/RESULTS.md) · [十六条测线](validation/v306/evidence/rays.json) · [下一项](validation/v306/NEXT_DECISION.md)。最靠近壁法向的边的平均夹角降低 8.30%，62 号条带最宽间距增大 26.61%，该候选在几何阶段结束。

- [v305 两档热敏珠网的分带热量收支](validation/v305/RESULTS.md) · [逐带热量与交换矩阵](validation/v305/evidence/budget.json) · [下一项](validation/v305/NEXT_DECISION.md)。珠部放热减少 318.303 µW，前中段占 58.8%，两档十个带的同阶段收支均闭合。

- [v304 带线热敏珠的第三档中段网距](validation/v304/RESULTS.md) · [三档同面数据](validation/v304/evidence/face_comparison.zip) · [下一项](validation/v304/NEXT_DECISION.md)。珠部放热再降 3.10%，4 mW 壁温升至 297.694199 K，热峰位于迎风引线 207 号面。

- [v303 带线热敏珠的 20D 与 30D 外边界对照](validation/v303/RESULTS.md) · [两种热边界的逐面数据](validation/v303/evidence/face_comparison.zip) · [下一项](validation/v303/NEXT_DECISION.md)。外边界移至 30 倍珠径后，给定 305.15 K 时整件放热降低 0.00427%，同一 4 mW 下壁温升高 0.00018847 K。

- [v302 带线热敏珠整件功率边界](validation/v302/RESULTS.md) · [半温差独立核验](validation/v302/evidence/half.json) · [下一项](validation/v302/NEXT_DECISION.md)。4 mW 时 10,272 面共用 297.562541 K 求解壁温，七区放热总和复现输入。

- [v301 带线热敏珠温差输运对照](validation/v301/RESULTS.md) · [全域温界与同迭代流场](validation/v301/evidence/frames.json) · [下一项](validation/v301/NEXT_DECISION.md)。一次有限斜率 MUSCL 试验收敛，412 个节点持续越过温度下限，关闭这一配置。

- [v300 带引线热敏珠的中段网距和引线效应](validation/v300/RESULTS.md) · [七区与四场](validation/v300/evidence/comparison.json) · [下一项](validation/v300/NEXT_DECISION.md)。带线珠部放热从 11.94750 降至 10.27423 mW，整件放热 10.87809 mW。

- [v299 光滑热敏珠中段过渡网加密](validation/v299/RESULTS.md) · [功率与热峰](validation/v299/evidence/comparison.json) · [下一项](validation/v299/NEXT_DECISION.md)。同一光滑原面上的放热降低 1.68505 mW，整珠平均 h 与实验相差 +31.04%。

- [v298 光滑与带线热敏珠的配对控制](validation/v298/RESULTS.md) · [三项功率](validation/v298/evidence/comparison.json) · [网距与额外放热](validation/v298/evidence/spatial.json) · [下一项](validation/v298/NEXT_DECISION.md)。同一光滑几何换用混合网的放热差为 +2.79986 mW，匹配方法下加入引线的珠部差为 −0.36120 mW。

- [v297 带两根引线的 B 型热敏珠全表面换热](validation/v297/RESULTS.md) · [七区壁压与热流](validation/v297/evidence/fine_surface.json) · [反向流动](validation/v297/evidence/nearwall.json) · [工作台](validation/v297/evidence/workbench.json) · [下一项](validation/v297/NEXT_DECISION.md)。3 m/s 对称稳态覆盖 10,272 个原面，整件放热 12.53996 mW，珠部平均换热比实验图读高 47.37%。

- [v296 民用热敏珠独立安装与全表面复算](validation/v296/RESULTS.md) · [安装方法](validation/v296/INSTALL.md) · [5,216 面同刻结果与整珠实验](validation/v296/evidence/installation_result.json) · [下一项](validation/v296/NEXT_DECISION.md)。程序冷态重建初场，独立工作台新算复现原面表，整珠放热 10.69725502 mW。

- [v295 热敏珠 6 m/s 五步工作台计算与中断接续](validation/v295/RESULTS.md) · [新算逐面对照](validation/v295/evidence/recovered.json) · [结果页](validation/v295/evidence/resumed_ui.json) · [下一项](validation/v295/NEXT_DECISION.md)。从侧栏接续真实中断，全部 5,216 面与参考表一致，整珠放热 10.69725502 mW。

- [v294 民用热敏珠全表面近壁网格对照](validation/v294/RESULTS.md) · [两网逐面统计](validation/v294/evidence/comparison.json) · [科研界面](validation/v294/evidence/ui.json) · [下一项](validation/v294/NEXT_DECISION.md)。56 层与 48 层在同一启动时刻覆盖 5,216 面，整珠放热相差 −0.2733%，热峰均在工作台面元 5211。

- [v293 热敏珠近壁加密与第一半步残差定位](validation/v293/RESULTS.md) · [网格与初场](validation/v293/evidence/grid.json) · [终场](validation/v293/evidence/first_state.json) · [下一项](validation/v293/NEXT_DECISION.md)。最内 8 个区间加密，首次装配差异定位到加密带；末窗焓未达标。

- [v292 热敏珠 6 m/s 同刻全表面时间步长对照](validation/v292/RESULTS.md) · [同面差值](validation/v292/evidence/comparison.json) · [科研界面](validation/v292/evidence/ui.json) · [下一项](validation/v292/NEXT_DECISION.md)。两个半步的固定端点热量账通过，5,216 面最大热流差 0.006122 W/m²。

- [v291 热敏珠第一半步的连续残差核对](validation/v291/RESULTS.md) · [实际终场](validation/v291/evidence/first_half.json) · [下一项](validation/v291/NEXT_DECISION.md)。240 次内迭代末尾有 18 行六项全部通过，第一半步未达 20 行要求。

- [v290 热敏珠 6 m/s 启动时刻的全表面压力与热流](validation/v290/RESULTS.md) · [同阶段热量账](validation/v290/evidence/surface.json) · [界面核对](validation/v290/evidence/ui.json) · [下一项](validation/v290/NEXT_DECISION.md)。5,216 个原面数值覆盖完整，适用时刻为距构造初场 7.809 µs。

- [v289 全周松弛与 6 m/s 热敏珠首次接受的物理步](validation/v289/RESULTS.md) · [200 次复核](validation/v289/evidence/first_step_200.json) · [下一项](validation/v289/NEXT_DECISION.md)。同一 7.809 µs 启动步在第 188 次内迭代形成首个连续窗，终点六项通过。

- [v288 同速基流与热敏珠全周首步](validation/v288/RESULTS.md) · [同节点连续性与焓右端](validation/v288/evidence/rhs.json) · [下一项](validation/v288/NEXT_DECISION.md)。末窗五项通过、焓未达标。首次装配的主要误差位于下游轴线。

- [v287 热敏珠 6 m/s 的实验点和全周首步](validation/v287/RESULTS.md) · [原生残差与完整场核对](validation/v287/evidence/first_step.json) · [下一项](validation/v287/NEXT_DECISION.md)。已排除 30° 扇区，147,456 单元的完整周向首步能在本机运行，三项残差尚未达标。

- [v286 热敏珠净功率边界与稳态壁温](validation/v286/RESULTS.md) · [细网核对](validation/v286/evidence/grid7p5.json) · [下一项](validation/v286/NEXT_DECISION.md)。4 mW 时两网均覆盖 5,216 面，稳态壁温相差 0.06820 K。

- [v285 B 型热敏珠可调给定壁温](validation/v285/RESULTS.md) · [中间温差与零温差原生核对](validation/v285/evidence/native.json) · [下一项](validation/v285/NEXT_DECISION.md)。两档网格支持 293.15–305.15 K，零温差的 5,216 个面元均有物理零热流。

- [v284 B 型热敏珠的五步主工作台](validation/v284/RESULTS.md) · [15° 冷启动原面核对](validation/v284/evidence/fresh_audit.json) · [下一项](validation/v284/NEXT_DECISION.md)。两档网格可确认边界、回读或计算；15° 新运行逐字节复现原面表。

- [v283 B 型热敏珠的周向细化与同面差值](validation/v283/RESULTS.md) · [定量对照](validation/v283/evidence/comparison.json) · [下一项](validation/v283/NEXT_DECISION.md)。从 15° 到 7.5°，平均换热系数增加 1.35%，迎流端放热峰值增加 6.18%。

- [v282 公开 B 型热敏珠全表面壁压与热流](validation/v282/RESULTS.md) · [论文图 7 与三网对照](validation/v282/evidence/comparison.json) · [下一项](validation/v282/NEXT_DECISION.md)。新增 5,216 个珠体原面，整珠放热 9.38210 mW，平均换热系数比实验读数高 15.7%。

- [v281 缺热原面的三份存档与同面导出](validation/v281/RESULTS.md) · [逐量核对](validation/v281/evidence/comparison.json) · [下一项](validation/v281/NEXT_DECISION.md)。目标 936 面上两份工程热流存档均缺项，研究热流有值 936 面，新增工程面积 0 m²。

- [v280 同一物理步的风洞圆柱对照](validation/v280/RESULTS.md) · [九分量与原面诊断](validation/v280/evidence/trial.json) · [回到单工况原面缺项](validation/v280/NEXT_DECISION.md)。固定 120 次、自适应 720 次均为 0/9 达标，工作台可切换归档，自绘模型本轮新增 0 m²。

- [v279 给定热边界的三维风洞圆柱](validation/v279/RESULTS.md) · [首个物理步](validation/v279/evidence/trial.json) · [下一项](validation/v279/NEXT_DECISION.md)。220 原面有工程壁压和热流，2160 面试算均保留未收敛身份，安装壁交线热通量已按圆柱标记复核。
- [v278 NASA Lewis 光滑圆柱换热实测对照](validation/v278/RESULTS.md) · [六站数据](validation/v278/evidence/stations.csv) · [下一项](validation/v278/NEXT_DECISION.md)。铜带平均 Nu 的六站相对均方根差 3.26%；有量纲热流、壁压和背面换热单列未知。
- [v277 CM-011 周向网格敏感性](validation/v277/RESULTS.md) · [锥面压力与热流曲线](validation/v277/evidence/cone_section.png) · [下一项](validation/v277/NEXT_DECISION.md)。10,080 面总放热变化 0.1604%，肩部 202 面的表压变化达到或超过 0.0002 Pa。
- [v276 同一低速试件的两套 STL 剖分](validation/v276/RESULTS.md) · [四个截面与原面积分](validation/v276/evidence/comparison.json) · [下一项](validation/v276/NEXT_DECISION.md)。10,080 与 14,140 面均有同工况全壁数值，原面积分功率相差 0.2075%。
- [v275 从 STL 在五步工作台完成低速 CM-011 计算](validation/v275/RESULTS.md) · [原面与截面核对](validation/v275/evidence/verification.json) · [下一项](validation/v275/NEXT_DECISION.md)。14,140 导入面和同一输入版本的云图、探针、两条截面 CSV 已核对。
- [v274 CM-011 低速 +10° 三维全壁计算](validation/v274/RESULTS.md) · [周向曲线](validation/v274/evidence/cone_ring.png) · [下一项](validation/v274/NEXT_DECISION.md)。14,140 原面和同一截面迎背流差均已复算，MATLAB 从空目录计算到导出通过。
- [v273 CM-011 三维等温试件全表面压力与热流](validation/v273/RESULTS.md) · [三维能量账](validation/v273/evidence/energy.json) · [下一项](validation/v273/NEXT_DECISION.md)。14,140 面映射核对通过，12D 三维扇区壁面与远场差额为壁面放热的 0.01142%。
- [v272 CM-011 全壁放热与远场流出的同阶段核对](validation/v272/RESULTS.md) · [工作台核对](validation/v272/evidence/workbench.json) · [下一项](validation/v272/NEXT_DECISION.md)。12D 细网原生边界差额为 8.900 µW，归档旧节点重建差额的四项贡献。
- [v271 CM-011 等温缩比试件的全表面计算](validation/v271/RESULTS.md) · [热流与压力图](validation/v271/evidence/wall_distribution.png) · [下一项](validation/v271/NEXT_DECISION.md)。低速细网 14,140 面均有压力与热流，远场热收支差额仍待解释。
- [v270 Apollo AS-202 原量热测点与单工况适用范围](validation/v270/RESULTS.md) · [图读数据](validation/v270/evidence/station.json) · [下一项](validation/v270/NEXT_DECISION.md)。区分 CA1253R 飞行散点、相邻压力孔和给定壁温下的工程预测。
- [v269 单工况同面覆盖与热功率](validation/v269/RESULTS.md) · [下一项区域交付](validation/v269/NEXT_DECISION.md)。常用入口、混合入口与三维结果侧栏共用原面统计，并核对部分、完整、零及负功率。
- [v268 等温微球的全面积压力与热流](validation/v268/RESULTS.md) · [下一项交付](validation/v268/NEXT_DECISION.md)。2,608 个原面及双向穿流接口通过解析和独立数值核对，自绘算例新增 0 m²。
- [v267 混合求解入口和连通域壁面结果](validation/v267/RESULTS.md)。工程方法与局部黏性求解器实际交换切面状态，低速两壁的双变量计算覆盖 2/2 m²；原箭体覆盖另列。
- [v266 单时刻热功率](validation/v263/SIDEWALL_TEMPERATURE.md#单时刻热功率v266) · [导出复算](validation/v263/evidence/snapshot_power/export_check.json)。主工况已算面功率与整表面未知量分列；公开风洞算例使用原壁曲面积。
- [v265 原面超域分类与编号查询](validation/v263/SIDEWALL_TEMPERATURE.md#区域工作台按原面编号核对v265) · [原面核对](validation/v263/evidence/edge_domain/original_check.json)。主工况 336 个直接超域面按源记录分列，936 个下游缺热面的数值身份保持可查。
- [v264 自绘模型的侧壁温度与上游来源](validation/v263/SIDEWALL_TEMPERATURE.md) · [NASA 风洞独立测站结果](validation/v263/TMR_RESOLVED.md) · [原面下一项](validation/v263/NEXT_DECISION.md)。在研 v263 阶段的本地代码交付。936 个缺热原面的保存路径和探针通过对照，风洞终态与实验按全部原始测站复核。
- [v263 主工况原面基线与侧壁缺口](validation/v263/RESULTS.md)。已核对联合覆盖和候选重叠，旧研究值的工程归类已更正，区域求解待完成。
- [v263 NASA 侧壁同站实验对照](validation/v263/TMR_JOINT.md)。80 个截断域固壁环带，12 个压力与热流共同位置，含项目两套网格及 WIND 数值参照。
- [v262 静止空气节点探针](STATIC_AIR_PROBE_V262.md)、[离线节点页](validation/v262/evidence/workbench/index.html)、[下一决定](validation/v262/NEXT_DECISION.md)。主窗口查询和下载按实际节点、帧、量的身份与反算方法核对。
- [v261 静止空气初始化](STATIC_AIR_FACTORY_V261.md)、[全节点核对](validation/v261/evidence/audit/report.json)、[下一决定](validation/v261/NEXT_DECISION.md)。两种配置的物性对象、温度回读与越界恢复通过本机单线程对照。
- [v260 量热反算](STABLE_CALORIC_INVERSE_V260.md)、[独立审核](validation/v260/evidence/audit02/report.json)、[下一决定](validation/v260/NEXT_DECISION.md)。195 个固定状态的温度精度达标，修复 60 K 附近的域误拒并完成端点对照。
- [v259报告](EILMER_MARCH_REVIEW_V259.md)、[最终回执](validation/v259/evidence/final_report.json)、[当时下一决定](validation/v259/NEXT_DECISION.md)。固定lmr块推进入口未实现、旧分支反馈/驻留语义不满足直接采用，候选拒绝；无新求解，本地交付。
- [v258报告](DIRECTION_CONSISTENCY_V258.md)、[最终回执](validation/v258/evidence/final_report.json)、[当时下一决定](validation/v258/NEXT_DECISION.md)。固定方向诊断完成，半步一致性未过、0更新；原面方向出口与本地交付，回H1具名区域候选复审。
- [v257报告](KINSOL_COUPLED_CONTROL_V257.md)、[最终回执](validation/v257/evidence/final_report.json)、[当时下一决定](validation/v257/NEXT_DECISION.md)。固定KINSOL完整耦合接口/出口通过，首轮线性精度拒绝，接受更新0；本地交付，组合关闭。
- [v256报告](FULL_DOMAIN_RESIDUAL_V256.md)、[最终回执](validation/v256/evidence/final_report.json)、[当时下一决定](validation/v256/NEXT_DECISION.md)。来源拒绝后实现完整成熟F(x)，零时间接口/原面诊断通过；尚无新非线性发展，本地交付。
- [v255报告](REGULARIZED_ENTROPY_INTEGRAL_V255.md)、[最终回执](validation/v255/evidence/final_report.json)、[当时下一决定](validation/v255/NEXT_DECISION.md)。驻点正则化及全原球帽目标数值对照通过；实际前体入口仍待闭合，本地交付。
- [v254报告](COUPLED_ENTROPY_INTEGRAL_V254.md)、[最终回执](validation/v254/evidence/final_report.json)、[当时下一决定](validation/v254/NEXT_DECISION.md)。变熵积分/46交界通过，全球帽起步未收敛；本地交付。
- [v253报告](ENTROPY_INLET_V253.md)、[最终回执](validation/v253/evidence/final_report.json)、[当时下一决定](validation/v253/NEXT_DECISION.md)。真实共同环入口与一次耦合拒绝；本地交付。
- [v252报告](PHASE_NONSIMILAR_PROFILES_V252.md)、[最终回执](validation/v252/evidence/final_report.json)、[当时下一决定](validation/v252/NEXT_DECISION.md)。是否仍为最新请读当前页。
- [公开实验模型资料](validation/public_benchmarks_20261116/README.md)。
- [工程算法报告交付](reports/20261112_algorithm_inventory/README.md)。
- [Harness第九批验收](../harness_tasks/acceptance_20261116_batch9/REPORT.md)。

## 原记录无损快照

- [原AGENTS完整原文（含v252）](validation/record_layout_20261116/AGENTS.pre_apply.md)：全部用户决定、阶段摘要、失败和存储记录。
- [原目标完整原文（含v252）](validation/record_layout_20261116/PROJECT_GOALS.pre_apply.md)：原能力表与执行记录。
- [开始整理时的v251快照及并发接续](validation/record_layout_20261116/README.md)。
- [更早目标历史](PROJECT_GOALS_HISTORY_20261031.md)、[存储/恢复说明](OUTPUTS_STORAGE.md)。
- [迁移SHA与来源](validation/record_layout_20261116/migration_manifest.json)。快照相对路径按原文件位置解释，不作现行指令。

## 原AGENTS版本入口

下表行号指含v252的原AGENTS快照，仅索引有明确版本标题的段落；未匹配的零散摘要和用户决定仍完整保存在快照。没有现存链接就不猜造。

| 版本 | 专题报告 | 当时下一决定 | 原AGENTS行 |
|---|---|---|---:|
| v252 | [报告](PHASE_NONSIMILAR_PROFILES_V252.md) | [决定](validation/v252/NEXT_DECISION.md) | 5 |
| v251 | [报告](NATIVE_ILU_ACCEPTANCE_V251.md) | [决定](validation/v251/NEXT_DECISION.md) | 13 |
| v250 | [报告](NATIVE_NEWTON_DEVELOPMENT_V250.md) | [决定](validation/v250/NEXT_DECISION.md) | 21 |
| v249 | [报告](NATIVE_NEWTON_CONTROL_V249.md) | [决定](validation/v249/NEXT_DECISION.md) | 29 |
| v248 | [报告](FAS_NATIVE_RECOVERY_V248.md) | [决定](validation/v248/NEXT_DECISION.md) | 37 |
| v247 | [报告](FAS_BOUNDED_DEVELOPMENT_V247.md) | [决定](validation/v247/NEXT_DECISION.md) | 45 |
| v246 | [报告](REAL_FAS_CONTROL_V246.md) | [决定](validation/v246/NEXT_DECISION.md) | 53 |
| v245 | [报告](FAS_COMPATIBILITY_V245.md) | [决定](validation/v245/NEXT_DECISION.md) | 61 |
| v244 | [报告](VISCOUS_REGIONAL_BOUNDARY_V244.md) | [决定](validation/v244/NEXT_DECISION.md) | 67 |
| v243 | [报告](FOREBODY_SPATIAL_COMPARISON_V243.md) | [决定](validation/v243/NEXT_DECISION.md) | 75 |
| v242 | [报告](FOREBODY_FIXED_WALL_PLC_V242.md) | [决定](validation/v242/NEXT_DECISION.md) | 83 |
| v241 | [报告](INTERNAL_VOLUME_STREAMTUBES_V241.md) | [决定](validation/v241/NEXT_DECISION.md) | 91 |
| v240 | [报告](FOREBODY_PLANAR_EULER_V240.md) | [决定](validation/v240/NEXT_DECISION.md) | 99 |
| v239 | [报告](NATIVE_MUSCL_RAMP_V239.md) | [决定](validation/v239/NEXT_DECISION.md) | 107 |
| v238 | [报告](GENERAL_HLLC_COMPATIBILITY_V238.md) | [决定](validation/v238/NEXT_DECISION.md) | 115 |
| v237 | [报告](PHASE_STEP_ACCEPTANCE_V237.md) | [决定](validation/v237/NEXT_DECISION.md) | 123 |
| v236 | [报告](RESTART_BRANCH_V236.md) | [决定](validation/v236/NEXT_DECISION.md) | 129 |
| v235 | [报告](SINGLE_DOMAIN_NATIVE_ENTRY_V235.md) | [决定](validation/v235/NEXT_DECISION.md) | 137 |
| v234 | [报告](HYBRID_MULTIZONE_COMPATIBILITY_V234.md) | [决定](validation/v234/NEXT_DECISION.md) | 143 |
| v233 | [报告](HYBRID_REGIONAL_CLOSURE_V233.md) | [决定](validation/v233/NEXT_DECISION.md) | 151 |
| v232 | [报告](HYBRID_PHASE_EULER_DEVELOPMENT_V232.md) | [决定](validation/v232/NEXT_DECISION.md) | 157 |
| v231 | [报告](HYBRID_PHASE_PRIMITIVE_RECOVERY_V231.md) | [决定](validation/v231/NEXT_DECISION.md) | 165 |
| v230 | [报告](HYBRID_PHASE_EULER_ROCKET_V230.md) | [决定](validation/v230/NEXT_DECISION.md) | 171 |
| v229 | [报告](HYBRID_PHASE_EULER_NATIVE_V229.md) | [决定](validation/v229/NEXT_DECISION.md) | 177 |
| v228 | [报告](HYBRID_PHASE_CONTINUATION_V228.md) | [决定](validation/v228/NEXT_DECISION.md) | 183 |
| v227 | [报告](HYBRID_PHASE_NOSE_HEATING_V227.md) | [决定](validation/v227/NEXT_DECISION.md) | 189 |
| v226 | [报告](HYBRID_PHASE_NOSE_EDGE_V226.md) | [决定](validation/v226/NEXT_DECISION.md) | 195 |
| v225 | [报告](HYBRID_PHASE_WAVES_V225.md) | [决定](validation/v225/NEXT_DECISION.md) | 201 |
| v224 | [报告](HYBRID_PHASE_LAMINAR_CHAIN_V224.md) | [决定](validation/v224/NEXT_DECISION.md) | 207 |
| v223 | [报告](HYBRID_SU2_FROZEN_BALANCE_V223.md) | [决定](validation/v223/NEXT_DECISION.md) | 213 |
| v222 | [报告](HYBRID_SU2_ROCKET_PHASE_V222.md) | [决定](validation/v222/NEXT_DECISION.md) | 219 |
| v221 | [报告](HYBRID_SU2_COMPATIBILITY_V221.md) | [决定](validation/v221/NEXT_DECISION.md) | 225 |
| v220 | [报告](HYBRID_SU2_PHASE_AIR_V220.md) | [决定](validation/v220/NEXT_DECISION.md) | 231 |
| v219 | [报告](HYBRID_SU2_ROUTE_V219.md) | [决定](validation/v219/NEXT_DECISION.md) | 237 |
| v218 | [报告](HYBRID_INDEPENDENT_REFERENCE_DECISION_V218.md) | [决定](validation/v218/NEXT_DECISION.md) | 243 |
| v217 | [报告](HYBRID_LOCAL_STENCIL_EVIDENCE_V217.md) | [决定](validation/v217/NEXT_DECISION.md) | 249 |
| v216 | [报告](HYBRID_COMPLETE_REGION_DECISION_V216.md) | [决定](validation/v216/NEXT_DECISION.md) | 255 |
| v215 | [报告](HYBRID_NATIVE_FILE_BRIDGE_V215.md) | [决定](validation/v215/NEXT_DECISION.md) | 261 |
| v214 | [报告](HYBRID_NATIVE_NK_OPERATOR_V214.md) | [决定](validation/v214/NEXT_DECISION.md) | 267 |
| v213 | [报告](HYBRID_FROZEN_NATIVE_OPERATOR_V213.md) | [决定](validation/v213/NEXT_DECISION.md) | 273 |
| v212 | [报告](HYBRID_SEQUENTIAL_REGIONS_V212.md) | [决定](validation/v212/NEXT_DECISION.md) | 279 |
| v211 | [报告](HYBRID_REGIONAL_NK_MEMORY_V211.md) | [决定](validation/v211/NEXT_DECISION.md) | 285 |
| v210 | [报告](HYBRID_PHASE_AWARE_EILMER_V210.md) | [决定](validation/v210/NEXT_DECISION.md) | 291 |
| v209 | [报告](HYBRID_ANISOTROPIC_ROUTE_DECISION_V209.md) | [决定](validation/v209/NEXT_DECISION.md) | 297 |
| v208 | [报告](HYBRID_REGIONAL_LAYER_COMPATIBILITY_V208.md) | [决定](validation/v208/NEXT_DECISION.md) | 303 |
| v207 | [报告](HYBRID_REGIONAL_VISCOUS_ENTRY_V207.md) | [决定](validation/v207/NEXT_DECISION.md) | 309 |
| v206 | [报告](HYBRID_PHASE_AWARE_AIR_V206.md) | [决定](validation/v206/NEXT_DECISION.md) | 315 |
| v205 | [报告](HYBRID_STANDARD_AIR_ROCKET_V205.md) | [决定](validation/v205/NEXT_DECISION.md) | 321 |
| v204 | [报告](HYBRID_PRESSURE_DOMAIN_V204.md) | [决定](validation/v204/NEXT_DECISION.md) | 327 |
| v203 | [报告](HYBRID_BOUNDED_THERMO_V203.md) | [决定](validation/v203/NEXT_DECISION.md) | 333 |
| v202 | [报告](HYBRID_NIST_AIR_V202.md) | [决定](validation/v202/NEXT_DECISION.md) | 339 |
| v201 | [报告](HYBRID_ROCKET_TEMPERATURE_LIMIT_V201.md) | [决定](validation/v201/NEXT_DECISION.md) | 345 |
| v200 | [报告](HYBRID_ROCKET_OPENFOAM_READY_V200.md) | [决定](validation/v200/NEXT_DECISION.md) | 351 |
| v199 | [报告](HYBRID_SYNCHRONOUS_CONTROL_V199.md) | [决定](validation/v199/NEXT_DECISION.md) | 357 |
| v198 | [报告](HYBRID_OPENFOAM_ENTRY_V198.md) | [决定](validation/v198/NEXT_DECISION.md) | 363 |
| v197 | [报告](HYBRID_NEWTON_KRYLOV_V197.md) | [决定](validation/v197/NEXT_DECISION.md) | 369 |
| v196 | [报告](HYBRID_GEOMETRIC_TIME_SCALE_V196.md) | [决定](validation/v196/NEXT_DECISION.md) | 375 |
| v195 | [报告](HYBRID_ROCKET_EULER_DEVELOPMENT_V195.md) | [决定](validation/v195/NEXT_DECISION.md) | 381 |
| v194 | [报告](HYBRID_ROCKET_EULER_PREPARATION_V194.md) | [决定](validation/v194/NEXT_DECISION.md) | 385 |
| v193 | [报告](HYBRID_ROCKET_EULER_MESH_V193.md) | [决定](validation/v193/NEXT_DECISION.md) | 389 |
| v192 | [报告](HYBRID_EULER_REFINEMENT_V192.md) | [决定](validation/v192/NEXT_DECISION.md) | 393 |
| v191 | [报告](HYBRID_EULER_CONTROL_V191.md) | [决定](validation/v191/NEXT_DECISION.md) | 397 |
| v190 | [报告](HYBRID_REFERENCE_ENTHALPY_V190.md) | [决定](validation/v190/NEXT_DECISION.md) | 401 |
| v189 | [报告](HYBRID_ROUTE_SELECTION_V189.md) | [决定](validation/v189/NEXT_DECISION.md) | 405 |
| v188 | [报告](EILMER_RESOURCE_FEASIBILITY_V188.md) | [决定](validation/v188/NEXT_DECISION.md) | 415 |
| v187 | [报告](EILMER_WORKSPACE_STORAGE_V187.md) | [决定](validation/v187/NEXT_DECISION.md) | 419 |
| v186 | [报告](EILMER_DENSE_FULL_DOMAIN_V186.md) | [决定](validation/v186/NEXT_DECISION.md) | 425 |
| v185 | [报告](EILMER_DENSE_STORAGE_V185.md) | [决定](validation/v185/NEXT_DECISION.md) | 429 |
| v184 | [报告](EILMER_EXPLICIT_FULL_DOMAIN_V184.md) | [决定](validation/v184/NEXT_DECISION.md) | 433 |
| v183 | [报告](EILMER_EXPLICIT_ENTRY_V183.md) | [决定](validation/v183/NEXT_DECISION.md) | 437 |
| v182 | [报告](EILMER_NATIVE_OPERATOR_REVIEW_V182.md) | [决定](validation/v182/NEXT_DECISION.md) | 441 |
| v181 | [报告](EILMER_REAL_NK_V181.md) | [决定](validation/v181/NEXT_DECISION.md) | 445 |
| v180 | [报告](EILMER_FULL_DOMAIN_DEVELOPMENT_V180.md) | [决定](validation/v180/NEXT_DECISION.md) | 449 |
| v179 | [报告](EILMER_RUN14_CONTROL_V179.md) | [决定](validation/v179/NEXT_DECISION.md) | 453 |
| v178 | [报告](EILMER_NATIVE_INTERFACE_V178.md) | [决定](validation/v178/NEXT_DECISION.md) | 457 |
| v177 | 查完整快照 | 查完整快照 | 461 |
| v176 | [报告](AEROTHERMAL_EVIDENCE_LEDGER_V176.md) | [决定](validation/v176/NEXT_DECISION.md) | 463 |
| v175 | 查完整快照 | [决定](validation/v175/NEXT_DECISION.md) | 467 |
| v174 | [报告](NACA_NATIVE_IDEAL_GAS_V174.md) | [决定](validation/v174/NEXT_DECISION.md) | 473 |
| v173 | [报告](NACA_SPATIAL_STEP_V173.md) | [决定](validation/v173/NEXT_DECISION.md) | 479 |
| v172 | [报告](NACA_ADAPTIVE_DEVELOPMENT_V172.md) | [决定](validation/v172/NEXT_DECISION.md) | 485 |
| v171 | [报告](NACA_SPATIAL_WALL_V171.md) | [决定](validation/v171/NEXT_DECISION.md) | 491 |
| v170 | [报告](NACA_SEGMENTS_MESH_V170.md) | [决定](validation/v170/NEXT_DECISION.md) | 499 |
| v169 | [报告](NACA_LAMINAR_ENTRY_V169.md) | [决定](validation/v169/NEXT_DECISION.md) | 505 |
| v168 | [报告](WALL_HEAT_ROUTE_V168.md) | [决定](validation/v168/NEXT_DECISION.md) | 511 |
| v167 | [报告](SST_ENERGY_CONSTRAINT_V167.md) | [决定](validation/v167/NEXT_DECISION.md) | 517 |
| v166 | [报告](FULL_DOMAIN_DEVELOPMENT_V166.md) | [决定](validation/v166/NEXT_DECISION.md) | 527 |
| v165 | [报告](CONSERVATIVE_STEP_V165.md) | [决定](validation/v165/NEXT_DECISION.md) | 537 |
| v164 | [报告](IDEAL_AIR_ADMISSIBILITY_V164.md) | [决定](validation/v164/NEXT_DECISION.md) | 547 |
| v163 | [报告](NATIVE_PROPERTY_EXCEPTION_V163.md) | [决定](validation/v163/NEXT_DECISION.md) | 555 |
| v162 | [报告](NATIVE_TERMINATION_OBSERVATION_V162.md) | [决定](validation/v162/NEXT_DECISION.md) | 563 |
| v161 | [报告](IDEAL_AIR_LINEAR_CONTROL_V161.md) | [决定](validation/v161/NEXT_DECISION.md) | 569 |
| v160 | [报告](IDEAL_AIR_FULL_DOMAIN_V160.md) | [决定](validation/v160/NEXT_DECISION.md) | 575 |
| v159 | [报告](AXISYMMETRIC_IDEAL_AIR_V159.md) | [决定](validation/v159/NEXT_DECISION.md) | 581 |
| v158 | [报告](PUBLIC_THERMAL_PREPARATION_V158.md) | [决定](validation/v158/NEXT_DECISION.md) | 587 |
| v157 | [报告](IDEAL_AIR_MAIN_WORKBENCH_V157.md) | [决定](validation/v157/NEXT_DECISION.md) | 593 |
| v156 | [报告](IDEAL_AIR_CHECKPOINTS_V156.md) | [决定](validation/v156/NEXT_DECISION.md) | 599 |
| v155 | 查完整快照 | [决定](validation/v155/NEXT_DECISION.md) | 605 |
| v154 | [报告](IDEAL_AIR_ENGINE_V154.md) | [决定](validation/v154/NEXT_DECISION.md) | 613 |
| v153 | [报告](FAST_IDEAL_AIR_V153.md) | [决定](validation/v153/NEXT_DECISION.md) | 619 |
| v152 | [报告](ACC_SBLI_CONTINUATION_V152.md) | [决定](validation/v152/NEXT_DECISION.md) | 625 |
| v151 | [报告](ACC_SBLI_INITIALIZATION_V151.md) | [决定](validation/v151/NEXT_DECISION.md) | 631 |
| v150 | [报告](ACC_SBLI_CONVERGENCE_V150.md) | [决定](validation/v150/NEXT_DECISION.md) | 638 |
| v149 | [报告](SST_JACOBIAN_PUBLIC_BENCHMARK_V149.md) | [决定](validation/v149/NEXT_DECISION.md) | 645 |
| v148 | [报告](SST_FLOW_FLUX_V148.md) | [决定](validation/v148/NEXT_DECISION.md) | 650 |
| v147 | [报告](SST_NATIVE_UPDATE_V147.md) | 查完整快照 | 655 |
| v146 | [报告](SST_STARTUP_RECONSTRUCTION_V146.md) | 查完整快照 | 661 |
| v145 | [报告](GENERAL_HYBRID_INITIALIZATION_V145.md) | 查完整快照 | 667 |
| v144 | [报告](BASE_CAP_ENVELOPE_V144.md) | 查完整快照 | 675 |
| v143 | [报告](TETGEN_PLC_EXTERIOR_V143.md) | 查完整快照 | 679 |
| v142 | [报告](TETGEN_FIXED_REPAIR_V142.md) | 查完整快照 | 683 |
| v141 | [报告](CONTINUOUS_PRISM_ENVELOPE_V141.md) | 查完整快照 | 687 |
| v140 | [报告](GRADED_LAYER_TERMINATION_V140.md) | 查完整快照 | 691 |
| v139 | [报告](MMG_SEED_OPTIMIZATION_V139.md) | 查完整快照 | 695 |
| v138 | [报告](MMG_LOCAL_REMESH_V138.md) | 查完整快照 | 699 |
| v137 | [报告](CGAL_LOCAL_CONTRACT_V137.md) | 查完整快照 | 703 |
| v136 | [报告](SU2_EXPLICIT_TRANSITION_V136.md) | 查完整快照 | 707 |
| v135 | [报告](SU2_PYRAMID_TRANSITION_V135.md) | 查完整快照 | 713 |
| v134 | 查完整快照 | 查完整快照 | 719 |
| v133 | 查完整快照 | 查完整快照 | 725 |
| v132 | 查完整快照 | 查完整快照 | 733 |
| v131 | 查完整快照 | 查完整快照 | 739 |
| v130 | 查完整快照 | 查完整快照 | 747 |
| v129 | [报告](SU2_DNS_HEAT_COMPARISON_V129.md) | 查完整快照 | 755 |
| v128 | [报告](SU2_SST_SPATIAL_POCKET_V128.md) | 查完整快照 | 759 |
| v127 | [报告](SU2_SST_CFL_COMPARISON_V127.md) | 查完整快照 | 762 |
| v126 | [报告](SU2_SST_MULTISTEP_V126.md) | 查完整快照 | 765 |
| v125 | [报告](SU2_SST_PHASE_V125.md) | 查完整快照 | 768 |
| v124 | [报告](SU2_SST_INTERFACE_V124.md) | 查完整快照 | 771 |
| v123 | [报告](SU2_HEOS_WORKBENCH_V123.md) | 查完整快照 | 774 |
| v122 | [报告](SU2_ACTUAL_COOLPROP_V122.md) | 查完整快照 | 777 |
| v121 | [报告](SU2_COOLPROP_CAPABILITY_V121.md) | 查完整快照 | 782 |
| v120 | [报告](SU2_ACTUAL_NEAR_WALL_V120.md) | 查完整快照 | 786 |
| v119 | [报告](SU2_ACTUAL_CONTINUATION_V119.md) | 查完整快照 | 790 |
| v118 | [报告](ACTUAL_ROCKET_FULL_DOMAIN_V118.md) | 查完整快照 | 794 |
| v117 | [报告](SU2_FROZEN_ENERGY_DECISION_V117.md) | 查完整快照 | 798 |
| v116 | [报告](SU2_COLLAR_STUDY_V116.md) | 查完整快照 | 802 |
| v115 | 查完整快照 | 查完整快照 | 805 |
| v114 | [报告](SU2_CLASSIC_THERMAL_BOUNDARY_V114.md) | 查完整快照 | 807 |
| v113 | [报告](SU2_HYBRID_ORIGINAL_SURFACE_V113.md) | 查完整快照 | 809 |
| v112 | [报告](SU2_MULTIGRID_CHECKPOINT_UI_V112.md) | 查完整快照 | 811 |
| v111 | [报告](SU2_NORMALIZATION_CHECKPOINT_V111.md) | 查完整快照 | 813 |
| v110 | [报告](SU2_CHECKPOINT_DIAGNOSTICS_V110.md) | 查完整快照 | 815 |
| v109 | [报告](SU2_ROCKET_STARTUP_SECTIONS_V109.md) | 查完整快照 | 817 |
| v108 | [报告](SU2_INITIALIZATION_MESH_V108.md) | 查完整快照 | 819 |
| v107 | [报告](EMBEDDED_WORKBENCH_V107.md) | 查完整快照 | 821 |
| v106 | [报告](SU2_ROCKET_INTERFACE_V106.md) | 查完整快照 | 825 |
| v105 | [报告](EMBEDDED_SU2_V105.md) | 查完整快照 | 827 |
