# v224：同相态空气的附着层流热链与工作台出口

2026-11-08。本轮已实现并测试**具名外缘包→同模型恢复/参考焓→层流积分→有符号热流→工作台JSON**，总目标active。没有继续整箭CFD，没有新真实箭体外缘或物理合格面积。v223及更早的未收敛/资源/物性失败保持。

## 1. 实际程序与调用位置

新`calcPhaseAwareLaminarHeating.m`通过`calcSurfaceHeating`的显式`cfg.heating.sidewallModel='phaseAwareLaminarIntegral'`分派，输入来自`mesh.phaseAwareEdge`。该包由程序内部外流生产器提供；日常最终使用目标仍不要求用户导入外部CFD。**当前生产器尚未接通**，本轮没有把旧NASA7结果或v222未收敛壁状态改标签冒充新外缘。

旧一键工况入口`runAeroHeatSingleCase`在选择此新模式时明确拒绝，防止它先跑旧NASA7外流再静默进入新气体热链。已有默认模式不变。直接热流API和主工作台数据API已贯通，不等同于新模式已在五步骤一键流程完整可用。

工作台入口`buildAerothermalViewerData`识别新结果并交给`buildPhaseAwareLaminarViewerData`，不混入mesh中可能保留的旧压力对象。输出8字段：q、Cf、参考性质动量厚度、恢复温度、参考温度、给定Tw、未知固体Tw、未知压缩性能量厚度。所有方法/来源/原因按字段保存，有限候选数值不增加物理合格面积。现有WebGL页面可“打开结果JSON”读取[制造平板控制](validation/v224/evidence/matlab01/workbench.json)；这不是箭体新解。

## 2. 同模型计算方法

气体固定`rah-phase-aware-dry-air-of13-v1`，同一Cp节点表及60—2000 K/p-T气相低密度联合域；对外hs(298.15 K)=0、es=hs−RT。内部SU2的+300000 J/kg坐标禁止进入工程外缘合同。

`phaseAwareAirState`新增可选`linearTable`输运：在已绑定Cp节点上取Lemmon2004稀薄mu/k，再线性插值，对应`rah-lemmon-dilute-icoTabulated-of13-v1`。默认`lemmonAnalytic`保持旧解析数值。`phaseAwareAirReferenceState`传递该选择，返回参考密度和输运身份。未重算物性全网格/原生setter控制或改变相态范围。

层流恢复因子采用显式工程近似r=√Pr_edge，参考Pr来自参考状态，而非强制0.72。复用Riley NASA-TM-107838附录H.6/H.7：

`h_aw = h_e + r U_e²/2`；`h_ref = 0.5(h_e+h_w)+0.22(h_aw−h_e)`。

同一焓表有界反解T_aw与T_ref；参考密度、mu、nu、k和Pr全部来自同一物性。参考性质Thwaites型近似保留既有形式：

`theta_ref² = 0.45/U_e^6 × ∫ nu_ref U_e^5 ds`；`lambda = theta_ref²/nu_ref × dU_e/ds`；`Cf = 2(0.2227+1.57 lambda−1.8 lambda²)/Re_theta`。

沿真实前驱链用梯形积分，速度梯度取最多4步割线。仅−0.09<lambda≤0.1使用当前拟合。Reynolds比拟`St=Cf/(2Pr_ref^(2/3))`，`q=St rho_e U_e(h_aw−h_w)`，正q进入壁面。负热流和接近零的合法热驱动保留，不裁剪。Cf/q均属于沿用的参考性质工程估算，不因气体一致就获得高超三维精度保证。

来源继承[v190参考焓](HYBRID_REFERENCE_ENTHALPY_V190.md)、[v189方法选择](HYBRID_ROUTE_SELECTION_V189.md)、[相态气体](HYBRID_PHASE_AWARE_AIR_V206.md)及新接口的实际源。当前未新增论文实验对照，Blasius只检查低阶平板公式极限，不能替代目标M5/20°热实验验证。

## 3. 输入与拒绝范围

合同`rah-phase-engineering-edge-v1`携带原坐标/面连接/排序SHA、逐面原编号、p/T/三分量U、来源SHA/角色/方法、当前有效/附着/层流/准二维掩码、特征长度与自由流平均自由程，以及绑定同一来源的前驱路径。路径必须严格上游，弧长闭合且不短于面中心直距，三维速度不能反向；标量掩码不能广播成全面有效。

`prescribed-control`只是给定控制；`unqualified-research`全包拒绝计算；`accepted-edge`须另有绑定相同源/几何的逐面接受记录。**本模块并未独立重建和证明输入生产器的接受证据**，输出始终physicalQualification=false，不能靠手填acceptedFace取得物理资格。

Kn采用`max(l_free,l_edge)/L`、阈值0.01，复用项目特征长度筛选思路，不是局部梯度Kn或过渡流解。未知/域外外缘、未知壁温、非附着层流、强三维、稀薄、分离、关联式越域及数值失败会阻断整条下游。种子必须是真实零距离起点；没有入口厚度模型时不接受非零距离种子。根面热流保持NaN，不能用0掩盖前缘奇异性。

`referenceMomentumThicknessM`明确是参考性质积分量。尚未取得速度/密度/焓剖面，故压缩性质量加权动量厚度、位移厚度、能量厚度和CFD入口剖面保持未提供，`profileReady=false`。一个Thwaites厚度不能唯一决定这些量。后体/底面/分离再附仍需要适用的区域模型/黏性反馈，当前模块不强行覆盖。

## 4. 本轮实际验证

使用81个沿平板路径的三角采样面，U=1000 m/s、p=10000 Pa、Te=220 K，给定Tw300/1000 K；这是**显式强制层流制造控制**，不预测其真实转捩。另有沿程变速/变温/非均匀Tw。没有把制造值画回原箭体或宣称形成新全模型解。

| 检查 | 实际结果 |
|---|---|
| 经过主`calcSurfaceHeating`分派 | 通过，80个非根面返回有限候选q |
| 平板Cf与0.664/√Re对照 | 相对差5.5545e−5，即0.0055545%；门槛2e−4 |
| MATLAB/Python独立语言的参考焓/反算/输运/厚度/Cf/q复算 | 三组最大缩放差8.5427e−16；门槛3e−9 |
| Tw1000 K负热流 | 保留，范围约−199616至−22317.8 W/m² |
| Tw取恢复温度的零热驱动 | q绝对值<1e−6 W/m²，无人工归零 |
| 分离、域外温度、未知温度、非层流、稀薄、未收敛来源 | 6组阻断/下游传播通过 |
| 原面/几何、能量坐标、严格上游图、掩码、断链、缺接受记录、非零起点、反向速度 | 9组畸形输入拒绝通过 |
| 主`buildAerothermalViewerData`分派与现有前端纯接口 | 8字段、原面积、NaN/给定/候选、旧压力隔离通过 |
| 修改边界的兼容性 | 默认解析输运数值相同；旧层流测试通过 |

[MATLAB结果](validation/v224/evidence/matlab01/report.json)、[独立核对](validation/v224/evidence/audit01/report.json)、[前端测试](validation/v224/evidence/audit01/frontend_test.txt)。独立语言复算仍共享公式及来源数据，是软件/关联式校核，不是实验或真实三维误差评估。

MATLAB唯一运行38.390秒、exit0、峰648654848 B；原环境locale/perl警告保留于原日志，不重跑换日志。独立Python/前端1.906秒，全部进程输入稳定/组清空。采样1CPU/4GiB；本轮无CFD/网格/原生查询/引擎构建/浏览器操作。CLI收口关闭。完整MATLAB替换文件直接提交仓库，不打工程包。29保护和并发工程候选不改、不复制其未提交源码。

## 5. 下一项是实际外缘生产器

新热链现在能算，但仍缺真实模型的同气体外缘。下一直接把已有内部激波/膨胀工程外流扩展到显式相态Air生产器，核对正常激波/弱斜激波/等熵膨胀的同模型质量、动量、能量和熵，再将**有适用依据的真实原面/路径**交给本接口。不能只把NASA7的p/T标签改成新气体，也不能从v222研究壁值直接读取Te。

复用原几何、区域和流向算法，附着工程适用范围保持；强三维/回流/底面仍由必要黏性及反馈负责。原面覆盖按面积报告新增候选与资格，源/模型/范围保持可查。当前没有新CFD/原生查询/构建准入。[下一执行点](validation/v224/NEXT_DECISION.md)规定完整推进顺序，避免停在新的物性小测试上。
