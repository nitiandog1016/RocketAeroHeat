# v64 局部变输运与SA候选合同

状态：**当前有界本构温度版本已完成定向MATLAB、独立Python及Node数据核对，真实公开网格达到原稳态标准。** 归档身份以`docs/validation/v64/summary.json`和同哈希的`v64_release_verification.json`为准。v63交付包及其证据冻结，不覆盖。当前真实箭体生产热流覆盖66.5626327424%，新增合格区域面积0，全面积目标未完成。

## 本次推进的真实方程

保留同一轴对称有限体积、HLLC内部通量、HLLE开放储库、Newton黏性应力、Fourier导热与压力/环向源。本次把层流输运从冻结入口常数改成局部Sutherland模型：

`mu(T)=muRef*(T/Tref)^(3/2)*(Tref+S)/(T+S)`，`k(T)=cp*mu(T)/Pr`。

EOS、守恒能量和声速继续使用定gamma/R；没有把NASA7比热混入原常比热总能量定义。NASA TMR官方[SA说明](https://tmbwg.github.io/turbmodels/spalart.html)和[运行说明](https://tmbwg.github.io/turbmodels/noteonrunning.html)采用定比热气体、Pr=.72和Sutherland输运，此受限基准切片与其一致。它不等于已认证全温域真实空气、化学/辐射或湍流。

`createTmrAxisymmetricAirGas.m`显式固定muRef=1.716e−5Pa·s、Tref=491.6°R×5/9、S=198.6°R×5/9。不得默换成常见273.15K/110.4K后声称逐字复现官方参数。generic测试可以使用其它明确参数，其角色仅为公式/算子校核。

保留WIND样例M7.11、总压364.8psi、总温1600°R和原静态储库。其导出Re/m约5.57254e6，与图示5.706e6约差−2.33887%；这些来自不同给定组合，不能暗调rho/p使两组条件同时“精确满足”。温域40..2000K只是数值接受域，未认证Sutherland在全区间的实验精度。

## 输运取样和守恒

- 内部面：取两侧有界线性温度迹的平均；每单元统一Barth–Jespersen式alpha，由本单元、面邻居和本次真实边界T的上下界确定。所有面左右迹（含周期坐标像）及平面源点共同约束alpha。没有用全局温域端点裁面温度，也没有裁单元Q；壁面严格用输入Tw。
- 原WLS热/黏性梯度和黏性功速度不作上述限制；该步骤只用于mu/k的本构取样，不是完整NS状态正性方法。保存原/有界面温度、左右迹、源温度、上下界和alpha，Newton每次完整残差重新计算。若最终本构温度仍超定义域则拒绝。
- 环向应力源：用平面质心的同一有界温度重构取mu。压力/环向源仍按平面面积；不能改成体积中心或r加权源。
- 扩散频率：取face、owner和other/实际边界迹温度所对应mu/k的最大值，除两侧最小密度；单元环向项使用cellMu。它是当前显式步长估计，不是任意非线性状态的无条件稳定性证明。
- 原四量耦合Newton对完整新残差作差分，系数对温度的变化纳入Jacobian。支持拓扑不扩到全局CFL或scale，不用旧常输运Jacobian冒充新方程。
- 没有显式transport的旧gas仍走原常数表达式，保留完整info逐值重放；显式gas.mu/k是参考温度的值，不能在显示中当作每个面元当前系数。

## 验证和科研界面

独立温点、本构墙热/剪切、斜网格环向应力、频率、真实保守FV均值的无源非等温Couette空间参考，以及有色/逐列/方向Jacobian检查服务于新算子。只做必要空间检查，不追加旧制造场长推进。守恒账保留有单位缺陷、逐分量物理预算和包含内部共享面消去的显式浮点界；漏掉柱源负对照必须失败。

真实推进从v63实际已收敛的1000格TMR状态起步：它只作新方程初场，旧壁热不能复用。保持原2e−6及连续三步、BE缺陷和累计四量账要求。新结果接入原三维研究页，逐面显示p/q、给定Tw、输运参考和Tw处mu/k；未知Te/Tr/求解Tw保持NaN。默认五步骤工作台和真实生产例子不自动改换模型。

首次实际重启PID52024在首个新残差处被温域保护拒绝，未推进。只读重算冻结Q发现cell T=80.005..703.915K，但原无限制面迹−1143.273..2803.329K，6个内部面超域；平面源迹79.936..703.492K。由此增加上述显式有界本构重构，未扩大接受温域，亦未将未解字段补值。此前无限制变输运core PID14836仅为历史，当前版本已重新测试并真实求稳态。

### 实际终态及验证边界（2026-10-20）

| 核对 | 实际证据 | 结论 |
| --- | --- | --- |
| 当前算子、旧常物性重放和Jacobian | MATLAB PID31620，exit0，源码稳定 | μ/k、墙热/剪切、柱源、频率、完整四量账与旧完整info重放通过 |
| 有界本构温度与真实TMR推进 | MATLAB PID27864，exit0，源码稳定 | 新方程19接受/6拒绝/24783次完整残差调用，原残差129926.557558降至1.89366780475e−8；原2e−6连续3步不变 |
| 独立核心与真实TMR算子 | Python PID32084/55260，exit0 | 由Q直接WLS重建、共享HLLC/HLLE通量/源/系数/频率独立复算；累计BE四量账7.12138e−15；未重演所有中间BE状态 |
| 新三维研究数据、旧导出回归 | MATLAB PID53828，exit0 | 40源墙面→3840共享顶点三角面，28项篡改拒绝；旧常物性JSON除来源绑定外不变、HTML逐字节不变 |
| JSON、几何、两截面、CSV和本构来源 | Node PID40104，exit0 | 两个真实开口、真曲面/弦面面积及积分分列；未知Te/Tr/求解Tw均null；未验收浏览器/WebGL运行 |
| 公开墙面测点诊断 | MATLAB PID41492，exit0 | 24独立压力点、27热流点原区间映射，16项拒绝；包含包内证据默认回退；Python PID33824独立原数据/映射/统计通过 |

无源非等温Couette只作空间网格参考，没有长推进：12/24/48格墙热RMS为8.9262%/4.9707%/2.6441%，剪切RMS为8.9001%/4.9563%/2.6365%。MATLAB有色/逐列Jacobian相对差3.6923e−11，方向差6.3513e−11；独立Python相应4.1705e−11/8.2768e−11。不是工程热流精度或实验验收。

归档前对已冻结MAT再执行独立核对：核心Python PID9968、真实TMR Python PID11876均由父运行器观测实际exit0/脚本稳定，未再推进长稳态。早先独立PID32084/55260为同脚本的首次成功记录；交付独立JSON绑定最后实跑身份。

真实TMR新稳态signed q范围2689.04798..112915.82158 W/m²。真旋转曲面积0.726522289308 m²，总热10177.7837152 W；显示弦面0.725947045224 m²，积分10168.9108228 W。q未因面积差被缩放。Tw=311 K时μ=1.89769076448e−5 Pa·s、k=.0264800342889 W/(m K)，与μref/kref分开。

公开数据诊断保持原始27行压力，末3行重复不增加权重；27热流点保留signed计算值与正加热实验幅值身份。按逐点相对偏差统计，压力平均绝对偏差20.5873%、相对RMS40.1293%；热流平均有符号偏差−56.9875%、平均绝对偏差58.2772%、相对RMS63.4477%。图保留实际粗网格的阶梯/热流起伏，不平滑。±10%是来源给出的估计，不是人为置信区间；层流/湍流、入口剖面、Re、出口和网格/外域差异仍存在，**本对照未通过实验精度验收**。

完整包`example_axisymmetric_v64/tmr/web/index.html`为新三维研究页；`example_axisymmetric_v64/comparison/tmr_wall_diagnostic.png`及CSV为诊断对照。默认运行`main`仍进入统一五步骤工作台。归档前独立运行和清单验证结果以发布核验JSON为准。

### 可复用的完整文件和运行

`axisymmetricFVTransport.m`负责显式本构；`axisymmetricFVTemperatureTrace.m`负责有界本构取样；`createTmrAxisymmetricAirGas.m`固定来源常数。原残差、完整Jacobian、RK2/隐式推进及TMR导出已直接更新。所有变更MATLAB文件按完整文件交付。

```matlab
addpath(pwd,fullfile(pwd,'tests'),fullfile(pwd,'internalFlow'), ...
    fullfile(pwd,'regionalFlow','axisymmetric'));
testAxisymmetricVariableTransport;
testAxisymmetricVariableJacobian;
testAxisymmetricTemperatureTrace;
studyTmrAxisymmetricSutherland;  % 实际重求局部输运稳态，不复用旧墙热
testTmrAxisymmetricSutherlandView;
testTmrAxisymmetricWallComparison;
```

测试初场优先读工作区`outputs/v63_validation`，完整包可回退到冻结`validation/v63`。查看和比较测试接受显式证据文件及输出目录，不要求重跑长稳态；包内隔离测试验证可移植性。

## 唯一湍流候选及后续实现边界

已从官方方程与源代码说明固定一个候选：`SA-noft2-nonnegative-AJS2012-v1`，非负分支、无ft2，采用TMR Note1(c)的2012 S修正，Prt=.90，不加专用高Mach压缩性修正。这里固定的是我方可复现候选；WIND样例只写`spalart`，未明确全部ft2/限制器默认，不能声称已精确复刻其版本。

给定入口/回流`nuTilde/nuInfinity=5`不是`muT/mu=5`。第五个变量若取`W=rho*nuTilde`，保守式须包含`−(nu+nuTilde)/sigma * grad(rho)·grad(nuTilde)`密度修正；轴对称共享面已有r权，不能额外重复加径向标量几何源。壁距须到真实壁边的最短距离；流体入口、外域、回流及壁面条件分别明确。

这些SA条目目前是来源/方程合同，未实现第五量输运、湍流壁热或实验验收。下一步实现同一共享质量通量、完整源/扩散及五量耦合接受账，随后核对公开入口剖面和p/q；未完成入口、网格/外域/边界敏感性前不晋升真实区域。零攻角外扩成功也不能自动认证大攻角31°收缩、后体或三维横流。
