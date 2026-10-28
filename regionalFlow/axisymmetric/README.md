# 贴体轴对称区域NS数值模块

本库沿用X主轴、SI及默认来流−X，在子午面上实际推进四个无旋NS守恒量`[rho,rho*ux,rho*ur,rho*E]`，v65新增`rho*nuTilde`第五输运量与SA耦合。它是接向公开区域基准的内部求解路线，不导入外部CFD解。当前已测常物性环隙解析稳态、公开五档原网格几何/静止场，以及1000格Mach7.11常物性/局部Sutherland层流和SA五变量数值稳态；**未完成TMR湍流实验验证、真实收缩/后体或大攻角三维区域验收，新增生产已解面积0。**

## 模块用途

| 完整文件 | 用途 |
| --- | --- |
| `loadTmrMach7Grid.m` | 从冻结原始ZIP读取五档双块公开网格、逐字节绑定、去重共享接口、转为X=−sourceX的SI节点；不读取流场解 |
| `buildAxisymmetricFVGeometry.m` | 凸四边形真实每弧度体积/面权重、平面和体积质心、共享面及显式i周期配对、WLS梯度矩阵 |
| `axisymmetricFVGas.m` | 明确gamma/R/mu/Pr/接受温域，导出cp/cv/k；不自动替换为空气模型 |
| `axisymmetricFVTransport.m` | 按显式常物性/Sutherland合同求局部μ/k，拒绝温域外输入，不裁温度 |
| `axisymmetricFVTemperatureTrace.m` | 逐格统一alpha限制仅用于本构μ/k的面/源温度迹，不限制原黏性/热梯度或Q |
| `createTmrAxisymmetricAirGas.m` | 冻结TMR官方Rankine参考常数、定cp及Pr，不调密度匹配另一组Re |
| `axisymmetricFVConservative.m` | 将`[rho ux ur T]`明确输入打包成四守恒量 |
| `axisymmetricFVState.m` | 从四量恢复rho/ux/ur/T/p/a；非法密度/温度拒绝而不裁值 |
| `axisymmetricNavierStokesResidual.m` | 一阶HLLC对流、Newton黏性与Fourier导热、环向应力/压力源、逐墙q/功和完整四量账 |
| `advanceAxisymmetricNavierStokes.m` | 真实全局步长SSP-RK2、声学与扩散限制、整步拒绝、边界和源的两阶段积分、独立稳态停止 |
| `axisymmetricFVJacobian.m` | 完整残差的稀疏有色有限差分，包含梯度模板、边界与柱坐标源，失败扰动和调用成本显式记录 |
| `advanceAxisymmetricSteady.m` | 全局后向Euler伪时间及四量耦合Newton，逐步/累计BE账、闭域质量约束与原残差稳态验收 |
| `exportAxisymmetricAnnularStudy.m` | 将已收敛周期环隙的真实墙面旋转并接入原WebGL合同；重放空间残差绑定气体/几何/墙值，未知温度保持NaN |
| `exportTmrAxisymmetricStudy.m` | 将真实收敛TMR圆柱/外扩墙映射为共享圆环三维研究页，保留原壁面编号、算法和signed q，真/弦面积及热积分分列 |
| `axisymmetricSAModel.m` | SA-noft2非负分支与AJS2012修正，局部源/涡黏性及固定常数；不是独立输运求解 |
| `axisymmetricSAWallDistance.m` | 到全部实际wall原节点线段的最短距，保留来源；起伏周期墙明确拒绝 |
| `axisymmetricSADiscretization.m` | 共享实际质量通量上的上风N、SA扩散、梯度平方和完整密度修正，NS有效mu/k |
| `axisymmetricRANSResidual.m` | 真五量耦合空间算子与完整每弧度五量账；前四量仍原EOS |
| `exportTmrAxisymmetricSAStudy.m` | 仅接受真实五量稳态的三维研究导出，给定N/nu、墙muT0、分子系数和signed q各自溯源 |

## 离散合同

全部体积和通量按每弧度计算，完整旋转体的面积/热功率另乘`2*pi`。每单元保存平面面积A、`V=integral(r dA)`和两种不同质心。直边面使用`integral(r*n ds)`及r加权采样点；满足`sum(Sx)=0`、`sum(Sr)=A`，所以常压静止不会出现虚假径向力。

空间算子为`V*dQ/dt = -sum[(Fc-Fv)*surfaceArea] + [0,0,(p-tauTheta)*A,0]`。压力和环向应力在平面质心求积，速度散度保留`dux/dx+dur/dr+ur/r`。共享面只算一次，以相反符号作用于两侧；径向动量带真实柱坐标源，不能用Cartesian的“仅边界四量相消”检查替代。移动壁的总能量外流为`qIntoWall-uWall·tauNormal`，热量和机械功分开。

Q是有限体积守恒均值；`EOS(Q均值)`并不自动等于真实平均温度或质心点值。测试分别比较独立积分的Q均值及连续剖面，未混淆三者。对流当前一阶，黏性梯度为WLS加完整位移约束，没有宣称完整NS二阶精度。

## 明确边界和运行

必须给出`boundary.iMin/iMax/jMin/jMax`。支持固定网格的等温无滑移切向运动墙、成对i周期、正交首层的r=0正则轴、全特征超声速入口/出口，以及显式零周向速度储库的特征边界。特征边界复用`internalFlow/resolveNSCharacteristicBoundary3D.m`的法向解，但四量接口不报告五特征数量。`prescribedTrace`只用于有完整解析迹的数值校核，不能充当缺失的真实亚声速入口。

轴边界的偏斜首层目前明确拒绝，需要独立Neumann梯度模板后才扩展。环隙稳态并不检验r=0推进。当前支持局部Sutherland输运与SA-noft2五量耦合稳态，EOS仍为定gamma/R；无周向速度、无三维横流、无转捩/化学/辐射/烧蚀，湍流模型尚未通过实验精度验收；不能把一个固定全局横流分量直接填成全周同号ur来代表攻角。

在工程根目录：

```matlab
maxNumCompThreads(1);
addpath(pwd,fullfile(pwd,'tests'),fullfile(pwd,'internalFlow'), ...
    fullfile(pwd,'regionalFlow','axisymmetric'));
testAxisymmetricFVGeometry;
testAxisymmetricNSResidual;
testAxisymmetricOpenBoundary;
testAxisymmetricNSAdvance;       % 实际推进8/16格到稳态，约一分钟级，不能用步数上限代替
testAxisymmetricAnnularView;     % 读取上一步真实证据，不重新长推进
testTmrMach7Grid;
testTmrAxisymmetricGeometry([],4:-1:0); % 五档公开网格，256000格只做静止空间校核
testTmrAxisymmetricStartup;      % 1000格20步，明确未收敛/未对实验
```

`testTmrAxisymmetricGeometry`的输出目录请省略参数以采用缺省目录；检查全部五档可写`testTmrAxisymmetricGeometry(fullfile(pwd,'outputs','v62_validation','tmr_geometry'),4:-1:0)`。

真实环隙导出读取`mesh`和`result`，不读取当前待编辑工况：

```matlab
s=load(fullfile(pwd,'outputs','v62_validation','advance','axisymmetric_advance_evidence.mat'));
source=s.solutions{end};
exportAxisymmetricAnnularStudy(source.mesh,source.result,fullfile(pwd,'outputs','annular_web'));
```

完整包的`example_axisymmetric/web/index.html`可离线查看。它的两片圆柱墙在周期端口开口，没有鼻部或端盖。原热流密度不为弦面面积差而缩放，真曲面积/功率和实际绘图三角面/功率分别输出，探针与CSV保持原q。默认`main`仍打开原五步骤科研工作台，本数值例子没有改变生产场量。

## 已测结果与后续

v63在同一离散上新增耦合稳态推进，原环隙残差调用988/1655，冻结历史RK2正确调用33517/160317（初态1次、每步2次；旧3N估算撤回）。公开1000格常物性层流Mach7.11实际达到原2e−6连续三步标准，末态残差1.6569e−8、累计四量账2.93e−15。开放面显式HLLE四量储库跨过旧特征边界声速障碍；HLLC完全单向分支不再检查未使用星态，混合区全部原正性护栏保留。当前显式TMR对照仍触发混合区保护，不报告该算例的加速比。见[完整范围与实际终态](../../docs/AXISYMMETRIC_STEADY_DRAFT.md)。

```matlab
testAxisymmetricFVJacobian;
testAxisymmetricImplicitStep;
testAxisymmetricRiemannBoundary;
testAxisymmetricHLLCUpwind;
testAxisymmetricSteadyAnnular;       % 原稳态阈值和10个连续步不变
testAxisymmetricSteadyAnnularView;   % 真实稳态、伪时间和未知温度身份
testTmrAxisymmetricContinuation;     % 复用已保存真实稳态，检查父链与累计账
testTmrAxisymmetricSteadyView;       % 3840面研究页、源面绑定及18项篡改拒绝
testAxisymmetricRK2WorkCount;        % 实际分析器观察7步15次完整残差调用
```

v63完整包的`example_axisymmetric_v63/tmr/web/index.html`展示实际收敛的TMR截短墙（40源墙面/3840显示三角面，两端开口），并不是带鼻部和底面的完整箭体。探针、X截面周向曲线及CSV保留原signed q、压力迹、给定Tw、算法与研究身份；Te/Tr/求解Tw保持未知。真曲面积0.726522289308m²、真总热4921.35966445W，与弦面0.725947045224m²/4917.08293649W分列。MATLAB、Node数据与拓扑核对通过；本阶段没有浏览器视觉或WebGL运行验收。

检查点试验完整文件为`tests/studyTmrAxisymmetricSteady.m`、`studyTmrAxisymmetricContinuation.m`与`studyTmrAxisymmetricExplicitComparison.m`。续算只调整Newton/CFL等工作参数，不允许放宽稳态/BE/守恒门槛，禁止覆盖父MAT。生产热流面积保持66.5626327424%，不因数值稳态增加真实箭体区域资格。

8/16格真实环隙在原`2e-6 1/s`残差门槛及10个连续步下稳态，墙热相对RMS误差0.372508%/0.0833112%，热功闭合约1e−10。人工气体R=1、mu=.5，不是空气或实验验证。Mach7.11原网格20步只验证内部短推进，时间1.17936e−7s，远小于0.000799148s穿越时间；有符号负热流是启动阶段的当前值，不是稳态环境。

v64局部Sutherland真实新稳态：19接受/6拒绝，原残差1.89367e−8，累计四量账7.12e−15；μref/Tref/S与墙面μ(Tw)/k(Tw)分开，40源墙面的新三维映射/两截面/CSV核对通过。公开p/q诊断偏差显著，层流与湍流身份明确；实验精度未通过。具体数值、复用接口和证据见[变输运阶段说明](../../docs/AXISYMMETRIC_VARIABLE_TRANSPORT.md)。

v65已实际实现上述第五量，并在原稳态门槛下达到1.0612e−7，累计35接受步/10拒绝/47673次残差、五量账8.96e−15。完整五量算子/Jacobian/BE端点独立复算及三维数据核对通过；相应完整文件与测试见[SA阶段说明](../../docs/AXISYMMETRIC_SA_COUPLING.md)。RK2接口仍是四量NS，五量目前由耦合BE伪时间推进，不能声称已实现五量物理时间历史。

下一步先核对真实SA粗网格两个肩部负墙热的非正交Dirichlet梯度与近壁网格，再检查公开21点入口u/Tt剖面、p/q和网格/外域/出口敏感性。实际大攻角横流、31°收缩、后体入口和分离/再附仍需独立三维和区域边界验收，不能靠零攻角外扩数值稳态增加覆盖。累计热量继续后置。
