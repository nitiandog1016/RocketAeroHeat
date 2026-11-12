# RocketAeroHeat工程算法与适用边界报告

编制日期：2026年11月12日。对象：项目负责人及计算方法开发人员。工程算法源码基线为v231提交66e52d708240762b6a3756de8b15044cec9d5dda，并纳入编制期间主线新增的v232完成记录。独立工作区材料另行标注，准确来源版本见文件清单。

本报告整理项目已经收集、实现或使用的主要气动热工程算法，给出实际公式、输入、适用边界与证据状态。工程核心包括外流近似、热物性、附着边界层、驻点与球鼻热流、底压与代表底热，以及分离再附和突起候选方法。成熟Euler及NS/RANS只说明它们与工程链的分工，不展开通用CFD离散内核。

结论：项目已有多种受限区域方法和数值校核，但尚无适用于当前真实M5、20度攻角箭体的完整全表面合格解。当前有效路线是低成本全局外流、工程边界层和必要的局部黏性反馈。公式存在、程序返回有限数值、数值复算通过、实验符合与真实区域资格是不同层级。

本轮只检查源码和既有资料并编制报告，没有运行MATLAB、CFD、引擎构建或新的物理验证。文献页码和访问状态沿用项目已有阅读记录；本轮未重新下载和逐页审核全部论文。并发的engineering54、harness_tasks及未交付工作区内容不视为已合入生产结果。

## 1 阅读规则与符号

### 1.1 状态和适用范围的读法

“主链”表示已有调用路径，不等于通过真实工况物理认证；“可选区域”表示必须显式选择或独立调用；“条件研究”表示输入压力、入口或几何仍依赖给定条件；“文献候选”表示资料已经收集，但尚无完整可用模块。

本报告将范围分成三层：来源的理论与实验条件、代码实际门限、本项目实际验收证据。代码中的Kn、Re、温度、角度及厚度门限通常用于筛选，不能合起来宣称为一个已经全面实验验证的矩形有效域。未取得原文定义时保留资料缺口，不以常见同名公式补成项目已经核实的公式。

### 1.2 坐标和物理量

主轴X、鼻尖朝正X，零攻角来流速度朝负X，全部使用SI。外法向记为n_out；壁法向热梯度使用从壁面指向流体的坐标y；正热流表示进入壁面。下标∞为自由流，e为边界层外缘，w为壁面，aw为绝热恢复状态，星号为参考性质，s为驻点或源文特定下标。不同小节独立定义相同字母，特别是球面极角θ、动量厚度θ_m和激波角β_s不得混用。

$$
q_\infty=\tfrac12\rho_\infty U_\infty^2,\qquad C_p=\frac{p-p_\infty}{q_\infty},\qquad H_0=h(T)+\tfrac12|\boldsymbol U|^2.
$$

热流q_w单位W/m²；壁剪应力τ_w单位Pa；比焓h和总焓H_0单位J/kg；比熵s单位J/(kg·K)；μ单位Pa·s；k单位W/(m·K)；ν=μ/ρ单位m²/s。C_f=τ_w/(ρ_e U_e²/2)，St使用质量通量与焓差定义，不能与温差换热系数混为一谈。

给定T_w、壁节点流体温度、求解固体壁温、T_e和T_aw分别记录。当前单时刻主要输出q_conv(T_w,given)；没有材料及独立热平衡时不能从热流反推出固体壁温。真实零热流和负热流允许存在，未知量保留NaN或null。

## 2 方法总览和当前状态

| 编号 | 方法族 | 当前角色 | 关键限制 |
| --- | --- | --- | --- |
| A01 | 定比热和NASA7空气 | 旧工程主链物性 | 固定组分；NASA7分支及低温近似须显式区分 |
| A02 | 相态标准Air与线性输运 | v224至v230显式新分支 | 60至2000 K且满足p与T联合域；没有解离模型 |
| A03 | 修正牛顿压力 | 主链迎风压力 | 背风和底Cp默认是给定条件，不是分离压力解 |
| A04 | 正激波与弱斜激波 | 独立内核及局部外流 | 附体弱支、同物性守恒；局部公式不自动构成三维场 |
| A05 | Taylor–Maccoll切锥 | 主链局部锥近似及经典锥 | 基本解为零攻角直圆锥，不含肩部和横流 |
| A06 | Prandtl–Meyer膨胀 | 主链局部凸转角 | 超声速简单波；超温域、分离和三维干扰拒绝 |
| A07 | 球鼻驻点熵外缘与Billig几何 | 可选区域及几何预览 | 拟合球鼻前缀；不能沿全后体固定驻点熵 |
| A08 | 准流线和压力驱动流管 | 主链图路径及条件研究 | 图路径不是流管；v96不生产未知压力 |
| A09 | 恢复焓与参考焓 | 主链及新相态链 | 同一Cp、焓零点与输运；恢复因子仍是工程近似 |
| A10 | Thwaites层流积分 | 主链及v224新分支 | 附着准二维；参考厚度不能充当CFD剖面 |
| A11 | Head与Ludwieg–Tillmann | 主链湍流工程法 | 原不可压卷吸形式；分离后无再附解 |
| A12 | 固定Re转捩与平板回退 | 主链筛选和对照 | 不是稳定性或高超声速转捩预测 |
| A13 | Sutton–Graves与余弦延拓 | 驻点估算及旧鼻部默认 | 余弦分布是项目低阶延拓，无独立局部精度保证 |
| A14 | Levy–Lees局部相似BVP | 可选球鼻、质量厚度及新相态鼻部 | 附着薄层、给定相容外缘；不含沿程非相似与真实三维横流 |
| A15 | v61非相似层流 | 冻结条件研究 | 积分账和入口问题未过，不计生产面积 |
| A16 | Kinney型参考流长 | 可选侧壁方法对照 | 仅替换原链接受面，入口有离散假设，新面积0 |
| A17 | SAND93底压与代表底热 | 独立关联及受限底压映射 | 底热是代表标量，不是全底分布 |
| A18 | 声速流线尾锥压力 | 条件压力对照 | 直尾锥0至21度；不产热流或分离泡 |
| A19 | 熵层吞入和共同激波层 | 条件研究及历史原型 | 几何、驻点、质量与动量闭合未形成真实整箭资格 |
| A20 | Lees/Kemp、热峰、腔体、突起关联 | 文献候选 | 原始定义、输入和实际空间责任不能省略 |
| A21 | 成熟Euler及NS/RANS | 混合路线的外流及困难区生产器 | 当前真实运行尚未通过稳态和热资格 |

## 3 空气模型与连续流筛选

### 3.1 定比热与旧NASA7固定组分空气

定比热分支采用理想气体与显式γ、R。NASA7分支按N2、O2、Ar质量分数0.755、0.232、0.013混合，程序使用的各组分分子量与系数见airThermalProperties.m。各段系数a_i是组分系数，R_i为组分气体常数。

$$
\begin{aligned} c_p/R_i&=a_1+a_2T+a_3T^2+a_4T^3+a_5T^4,\\ h/(R_iT)&=a_1+a_2T/2+a_3T^2/3+a_4T^3/4+a_5T^4/5+a_6/T,\\ s^\circ/R_i&=a_1\ln T+a_2T+a_3T^2/2+a_4T^3/3+a_5T^4/4+a_7.\end{aligned}
$$

$$
\rho=\frac p{RT},\qquad s(T,p)=s^\circ(T)-R\ln(p/p_{ref}),\qquad a^2=\gamma(T)RT,\quad \gamma=\frac{c_p}{c_p-R}.
$$

代码域为80至2500 K。150至200 K使用多项式低温外推；80至150 K在150 K锚定振动冻结近似：h=h_150+c_p,150(T−150)，s°=s°_150+c_p,150 ln(T/150)。这是稀薄、保持气相的固定组分近似，并非已经完成露点或真实气体相态审查。没有解离、电离、化学非平衡或辐射。

旧生产airThermalProperties在1000 K有有限系数积分接缝。v95的surfaceFlowAirThermalProperties使用continuous-project-cp-v1，保持Cp而连续匹配h、s积分常数。两者总态必须显式转换，不能只凭“NASA7”名称互换。[R01、R11]

Sutherland旧输运与固定Pr热导率为：

$$
\mu(T)=1.458\times10^{-6}\frac{T^{3/2}}{T+110.4},\qquad k=\frac{\mu c_p}{Pr}.
$$

旧主链Pr默认0.72，可配置。固定Pr和Sutherland不代表相态新链的线性Lemmon输运。算法名称相同也必须核对具体物性和能量参考。

### 3.2 新相态标准Air模型

phaseAwareAirState采用模型rah-phase-aware-dry-air-of13-v1。固定摩尔组分N2/O2/Ar为0.7812/0.2096/0.0092，R=287.0491267875601 J/(kg·K)。3882节点的Cp(T)分段线性，焓和标准熵按同表解析积分；h_s(298.15 K)=0，e_s=h_s−RT，s(298.15 K,100 kPa)=0。[R02]

$$
\begin{aligned}c_p(T)&=a_i+b_iT,\\ h(T)-h(T_i)&=a_i(T-T_i)+\tfrac12b_i(T^2-T_i^2),\\ s^\circ(T)-s^\circ(T_i)&=a_i\ln(T/T_i)+b_i(T-T_i).\end{aligned}
$$

域为60≤T≤2000 K和0<p≤p_max(T)。T_c=132.6312 K，p_c=3785020 Pa，低温露点线按Lemmon2000表10及项目固定资料计算。

$$
p_{max}=\min\{10^6,\;10^4(T/100)^3\}\ {\rm Pa};\qquad T\le T_c:\quad p_{max}\leftarrow\min\{p_{max},0.1p_{dew}(T)\}.
$$

$$
\ln\frac{p_{dew}}{p_c}=\frac{T_c}{T}\left[-0.1567266\vartheta^{1/2}-5.539635\vartheta+0.7567212\vartheta^{5/2}-3.514322\vartheta^4\right],\quad\vartheta=1-T/T_c.
$$

0.1露点距离和1 MPa等界限是项目事前气相与低密度工程保护，不是文献给出的统一精度定理。高温来源包含预测数据，不能把整个域称为实验认证域。

低密度输运来自Lemmon与Jacobsen2004。以下μ单位Pa·s、k单位W/(m·K)，温度用K；28.9586是输运原式的g/mol数值，不能与EOS组分分子量静默合并。

$$
\begin{aligned}z&=\ln(T/103.3),\\ \Omega&=\exp(0.431-0.4623z+0.08406z^2+0.005341z^3-0.00331z^4),\\ \mu_0&=2.66958\times10^{-8}\frac{\sqrt{28.9586T}}{0.36^2\Omega},\\ k_0&=0.001308(10^6\mu_0)+0.001405(T_c/T)^{-1.1}-0.001036(T_c/T)^{-0.3}.\end{aligned}
$$

新工程链linearTable在同一Cp节点上计算μ_0、k_0后线性插值；默认lemmonAnalytic直接求原解析式。Pr(T)=c_pμ/k随温度变化。新相态链不自动替代旧NASA7默认链。

SU2内部为避免能量参考与成熟恢复流程冲突使用e_int=e_phys+300000 J/kg。对外恢复物理量时必须减去该平移；逐节点能量残差也必须先转换，再求范数：R_E,phys=R_E,int−300000 R_ρ。能量坐标平移本身不产生物理改善或收敛资格。[R16]

### 3.3 Kn连续流筛选

$$
\ell=\frac\mu p\sqrt{\frac{\pi RT}{2}},\qquad Kn=\frac{\max(\ell_\infty,\ell_e)}L<0.01.
$$

主链优先使用给定连续流特征长度或核验鼻半径，缺少时可用横向半径替代并记录来源；不是用三角形尺寸作为物理尺度。门限只能进一步收紧。新相态链同样按明确长度筛选。

平均自由程模型和Kn通过只说明选定尺度上的工程筛选；没有局部梯度Kn、滑移/温跳、过渡流或DSMC解。未知外缘和域外流向占位不能被标为真实稀薄流，也不能因占位Kn较小放行热流。上游未知和真实Kn超限分别记录并传播。[源码：airMeanFreePath.m、assessContinuumRegime.m]

## 4 压力与工程外流

### 4.1 修正牛顿压力

令μ_n=−n_out·U_hat∞。迎风面μ_n>0采用以下压力；μ_n≤0处的背风Cp和非迎风尾底Cp由配置给定，默认0，必须标为输入条件。

$$
C_p=C_{p,max}\max(\mu_n,0)^2,\qquad p=p_\infty+q_\infty C_p,\qquad C_{p,max}=\frac{p_{02}-p_\infty}{q_\infty}.
$$

calcCpMax以定比热正激波后等熵驻止求p_02。即使选择热流NASA7模式，此压力模块也不会自动变成变Cp压力模型。M≥1只是计算定义域；方法主要服务高超声速迎风工程估算，不解析遮挡、激波干扰、黏性压力、凹面或三维背风分离。[源码：modifiedNewtonianPressure.m、calcCpMax.m；R01]

给定背风Cp=0表示p=p∞，不是求得零压力。选受限SAND93底压时，拒绝底面保持未知，不退回给定Cp。压力计算、压力给定和压力未知分别统计。

### 4.2 正激波与同物性守恒

定比热正激波为：

$$
\begin{aligned}\frac{p_2}{p_1}&=\frac{2\gamma M_1^2-(\gamma-1)}{\gamma+1},\quad \frac{\rho_2}{\rho_1}=\frac{(\gamma+1)M_1^2}{(\gamma-1)M_1^2+2},\\ M_2^2&=\frac{1+(\gamma-1)M_1^2/2}{\gamma M_1^2-(\gamma-1)/2},\quad p_{02}=p_2[1+(\gamma-1)M_2^2/2]^{\gamma/(\gamma-1)}.\end{aligned}
$$

NASA7及新相态分支不用“平均γ替换”代替完整能量守恒，而求压缩根。n为激波法向，t为切向：

$$
\rho_1U_{n1}=\rho_2U_{n2},\quad p_1+\rho_1U_{n1}^2=p_2+\rho_2U_{n2}^2,\quad h_1+U_{n1}^2/2=h_2+U_{n2}^2/2,\quad \boldsymbol U_{t2}=\boldsymbol U_{t1}.
$$

要求Mn1>1、压缩支、熵增和两侧热力域合法。新相态v225显式去掉零强度平凡根，在量热合法区间寻根，最终再审p/T联合域。11保存状态完成独立守恒核对，但完整MATLAB套件的失败与终态缺口保持；不能写成完整套件已重跑通过。[R03]

### 4.3 弱斜激波与局部切楔

$$
\tan\delta=2\cot\beta_s\frac{M_1^2\sin^2\beta_s-1}{M_1^2(\gamma+\cos2\beta_s)+2},\quad M_{n1}=M_1\sin\beta_s,\quad M_2=\frac{M_{n2}}{\sin(\beta_s-\delta)}.
$$

δ为壁面压缩转角；求附体弱支并保持切向速度。变Cp用同一法向守恒和向量转角搜索，不能把墙面转角当激波穿越角。超过最大附体转角拒绝；新相态若因温域无法取得完整极线，返回“极线不完整”，不能当作已证明脱体。

基本模型是二维直楔局部近似。旧solveExternalFlow对三角面局部角度选择切楔，不能由此推出真实连续三维外流或跨肩部总压继承。真实鼻尖可能处于脱体/亚声速区，不能从鼻尖开始强行用附体弱波。[源码：solveObliqueShock.m、solveThermallyPerfectObliqueShock.m、solvePhaseAwareWave.m]

### 4.4 Taylor–Maccoll切锥

对零攻角、无黏、定比热、直圆锥附体弱支，速度以sqrt(2H_0)归一化，θ为从锥轴量起的极角，令A=(γ−1)(1−v_r²−v_θ²)/2：

$$
\frac{dv_r}{d\theta}=v_\theta,\qquad \frac{dv_\theta}{d\theta}=\frac{v_\theta^2v_r-A(2v_r+v_\theta\cot\theta)}{A-v_\theta^2},\qquad v_\theta(\theta_c)=0.
$$

从候选激波角跃迁后积分到锥面，射击满足无穿透。solveThermallyPerfectConicalShock为显式变Cp对照，不以定γ式冒充同物性解。旧外流将局部曲率合格面近似成锥面属于额外工程假设；直锥基本解不支持有攻角完整锥流、钝鼻熵层、肩部或收缩相互作用。

### 4.5 Prandtl–Meyer膨胀与顺序继承

$$
\nu(M)=\sqrt{\frac{\gamma+1}{\gamma-1}}\tan^{-1}\sqrt{\frac{\gamma-1}{\gamma+1}(M^2-1)}-\tan^{-1}\sqrt{M^2-1},\quad\nu(M_2)=\nu(M_1)+|\delta|.
$$

定Cp情况下总温、总压保持；变Cp时沿同总焓求速度、用同熵求压力，并按简单波积分转角：

$$
U(T)^2=2[H_0-h(T)],\quad s(T,p)=s_1,\quad |\delta|=\int_{T_2}^{T_1}\frac{c_p(T)\sqrt{M(T)^2-1}}{U(T)^2}\,dT.
$$

仅超声速局部凸转角；近无限Mach角极限、温域外、回流、强横流和分离不能由简单波继续。v225的marchPhaseAwareWavePath只支持同平面三维向量表达的二维折线，逐段继承p/T/U；它不是三维流管。

旧主链各面可能独立从自由流重算膨胀。原336个膨胀温域拒绝根及5483面下游缺口说明局部算法超域，不证明真实全局外流降到该低温，更不证明该处真实分离。不能通过降低物性下限或重播种填补。[R03、R17]

### 4.6 旧外流的受限脱体闭合

旧solveExternalFlow还保留局部切楔脱体时的工程闭合：以修正牛顿压力及局部法向激波总状态恢复热力量。当前限定光滑前体前10%轴长、凸且受限曲率的候选面；逐面仍保留独立方法ID。其具体筛选和恢复以源码为准，不是通用脱体激波层解。

该分支不生产共同激波形状、相邻流线熵输运或真正三维压力；热流后续仍是受限附着工程估算。新相态链不因旧分支有限而自动接受它。[源码：solveExternalFlow.m；R01]

## 5 球鼻外缘与驻点热流

### 5.1 驻点熵球鼻外缘

已核验球面中，θ从逆来流方向量起。正激波后的驻点态为p_02、T_02；球鼻压力及同熵/总焓恢复为：

$$
p_e=p_\infty+(p_{02}-p_\infty)\cos^2\theta,\quad s^\circ(T_{02})-s^\circ(T_e)=R\ln(p_{02}/p_e),\quad U_e=\sqrt{2[h(T_{02})-h(T_e)]}.
$$

NASA SP-8062第4.2.1节支持钝鼻用驻点熵恢复的思路，同时指出后体需要分析较弱斜激波来源熵。p_e的修正牛顿分布是本项目明确选择，不是报告对任意鼻形的精确解。

solveSphericalNoseEdge与新solvePhaseAwareNoseEdge均要求重新核验球面前缀、完整面元、法向、特征和拓扑。合成入射角最多30度、极角默认60度且最大80度是项目筛选；新分支还检查偏置驻点至目标面的完整拟合短弧。速度方向是局部切向近似，不能自动成为已验证三维流管。

v226新相态原面外缘及v227/v228热候选对应1306面、0.003518440977521075 m²，占全表面0.1790893%；它们是工程候选，新物理资格0。鼻部以外不得沿用同一驻点熵填完整后体。[R04、R05]

### 5.2 Billig球鼻激波几何预览

$$
\frac\Delta{R_n}=0.143\exp(3.24/M_\infty^2),\qquad \frac{R_c}{R_n}=1.143\exp\left[\frac{0.54}{(M_\infty-1)^{1.2}}\right].
$$

estimateBowShockSurface采用Billig1967球形钝鼻分支，远场斜率取Mach角。代码仅零攻角、3≤M≤17、近似定比热空气及r≤R_n的鼻端预览，不求全箭激波、压力或热流。不能将图上激波面当作具有物理匹配资格的外缘生产器。[R18]

### 5.3 Sutton–Graves驻点关联

$$
q_{s,ref}=1.83\times10^{-4}\sqrt{\rho_\infty/R_n}\,U_\infty^3,\qquad q_s=q_{s,ref}\frac{h(T_0)-h(T_w)}{h(T_0)-h(T_{ref,w})}.
$$

第一式系数对应代码的SI组合，q单位W/m²，默认参考壁温300 K；系数和参考壁温可配置。第二式是项目追加的同焓差比例壁温修正，不能全部称为Sutton–Graves原始关联。定Cp极限变成温差比。

要求超声速、可信物理鼻半径和Kn∞,Rn<0.01。未核验球鼻的给定半径仅允许合成角≤5度；通过球面和偏置驻点深度检查时可到≤30度。还检查候选区唯一性、非平封盖和法向对齐。该式提供单点估算，不含催化、解离、辐射及局部球面分布验证。原始完整经验工况包络未由本轮重建，因此不能把M>1程序检查说成全超声速实验有效域。[源码：calcStagnationHeating.m]

### 5.4 驻点余弦空间延拓

$$
q_w(\theta)\approx q_s\cos\theta.
$$

calcBluntNoseHeating在核验球鼻上作低阶延拓；不同给定壁温时使用同焓差/等效Cp修正换热系数。要求入射角≤30度，轴向和相对偏置驻点极角均≤80度，球面顶点误差≤3%、法向偏差≤20度，排除拓扑问题。

这是旧默认鼻模型cosineExtension，通常补原侧壁未解处。它没有解局部激波层或边界层，不应标成Levy–Lees解或局部实验认证。显式选局部相似法时，拒绝球鼻不能退回余弦填充。[R01；源码：calcBluntNoseHeating.m]

## 6 附着边界层与工程热流

### 6.1 图准流线和起点责任

traceSurfaceStreamlines按共享边、有向速度和Dijkstra成本构造前驱树，提供路径距离和真实几何步长。每个面沿自己的前驱积分；路径连通不证明流体质量守恒，也不提供相邻流线横向尺度。

起点零积分、零厚度是模型入口假设。不能把每个断连位置重新设置成自由流/零厚度，以获得更多有限面。失效传播、真实流动分离和数值断连必须分别说明。复杂带突起几何还需独立检查遮挡与特征保留。

### 6.2 恢复焓与Eckert参考焓

$$
h_{aw}=h_e+rU_e^2/2,\quad r_l=\sqrt{Pr_e},\quad r_t=Pr_e^{1/3},\quad h^*=\tfrac12(h_e+h_w)+0.22(h_{aw}-h_e).
$$

由同一h(T)反求T_aw与T*，再用p_e求ρ*和输运。Riley NASA-TM-107838附录H.6/H.7为项目已核对的原页来源。恢复因子仍是层流/湍流工程估计，并非实际绝热壁BVP解。[R06]

变Cp时参考焓反算与线性温度平均不等价。旧温度式T*=(T_e+T_w)/2+0.22(T_aw−T_e)只在定Cp时与参考焓式等价。当前airReferenceTemperature默认referenceEnthalpy；linearTemperature保留显式旧版对照。Kinney可选模块仍直接使用线性参考温度，见第6.8节，不能因名称写“参考焓”便认为它已经切成同一公式。

$$
\left.\frac{dT^*}{dT_w}\right|_{edge,aw}=\frac{c_p(T_w)}{2c_p(T^*)}.
$$

导数用于固定外缘/恢复状态的壁温耦合；端点温域必须逐个合法，不能用一个合法平均温度掩盖无效壁温或恢复温度。新相态链采用同一Cp表与linearTable输运，不能混用旧NASA7的H_0或S_0。

### 6.3 Thwaites型层流积分

$$
\theta_{ref}^2=\frac{0.45}{U_e^6}\int_0^s\nu^*U_e^5\,ds,\quad Re_\theta=\frac{U_e\theta_{ref}}{\nu^*},\quad\lambda=\frac{\theta_{ref}^2}{\nu^*}\frac{dU_e}{ds}.
$$

$$
l(\lambda)=0.2227+1.57\lambda-1.8\lambda^2,\qquad C_f=\frac{2l(\lambda)}{Re_\theta},\qquad -0.09<\lambda\le0.1.
$$

calcBoundaryLayerState沿前驱链梯形积分，速度梯度默认取最多4步割线。λ≤−0.09只标疑似层流分离；λ>0.1是拟合域外。它原为附着准二维工程闭合，代入参考ν不等于完整可压缩变换，更不是三维横流边界层。

v224 calcPhaseAwareLaminarHeating将同样的积分结构接入显式新相态合同，根面不补零，6类失效沿下游阻断。制造平板与Blasius极限、独立语言复算已过，但尚无真实前体可信外缘和入口全链资格；physicalQualification=false。输出referenceMomentumThicknessM仅为参考性质厚度，profileReady=false，不能据此给δ*、能量厚度或CFD速度/温度剖面。[R07]

### 6.4 固定阈值转捩

$$
Re_s=\frac{\rho_eU_es}{\mu_e},\qquad \chi_i=\left(Re_{s,i}\ge Re_{crit}\right)\ \lor\ \chi_{parent(i)},\quad Re_{crit}=5\times10^5\ \text{默认}.
$$

calcTransition在合法外流前驱上保持已转捩标签；热流输出缺失不必丢失标签，但标签不能跨外流/路径/连续流失效传播。阈值可配置，0/1间歇因子是类别选择，不是自然转捩区模型。

没有扰动谱、粗糙度、第二模态、横流稳定性或e^N分析；不能由Re_s≥阈值宣布真实高超声速转捩位置已经确定。外部并发候选有入口和状态继承修订，见第11节，未自动覆盖当前根目录公式。

### 6.5 Head卷吸与Ludwieg–Tillmann湍流法

令H=δ*/θ_m，H_1为卷吸形状量，E为无量纲卷吸率。当前实际源码为：

$$
\frac{d\theta_m}{ds}=\frac{C_f}{2}-(H+2)\frac{\theta_m}{U_e}\frac{dU_e}{ds},\qquad \frac{d(U_e\theta_mH_1)}{ds}=U_eE.
$$

$$
\begin{aligned}C_f&=0.246\,10^{-0.678H}Re_\theta^{-0.268},\quad E=0.0306(H_1-3)^{-0.6169},\\ H_1&=3.3+0.8234(H-1.1)^{-1.287},&&H\le1.6,\\ H_1&=3.3+C_h(H-0.6778)^{-3.064},&&H>1.6,\\ C_h&=\frac{0.8234(1.6-0.6778)^{3.064}}{(1.6-1.1)^{1.287}}.\end{aligned}
$$

高H分支C_h是项目为H=1.6接缝连续化而设置的系数，不是应默认归为原论文数值。段内U、ν线性插值，RK4步长加倍控制；H≤1.1、数值失败、Reθ域外分别处理；H≥2.4只标疑似准二维分离，整条下游无再附解。

当前代码采用300≤Reθ≤10^6；内侧300至600及8×10^5至10^6与平板状态smoothstep衔接，衔接同时作用θ/H/Cf。w=smoothstep(a)smoothstep(b)，a=(Re_min−300)/300，b=(10^6−Re_max)/(2×10^5)，a、b限制在0至1；smoothstep(t)=t²(3−2t)。带宽、H=1.4启动和回退均为工程选择。

重要来源边界：项目方法归属为Head及Ludwieg–Tillmann，但目前并发来源复核未取得这些具体幂律和Re门限的原始方程页。因此这里明确是“当前已实现公式”，不宣称所有系数已经逐项核对首发原文。不可压卷吸形式仅换参考ν，不自动取得高超声速密度加权厚度、完整能量或三维分离资格。[R08；源码：marchHeadBoundaryLayerSegment.m、calcTurbulentBoundaryLayerState.m]

### 6.6 固定H与局部平板对照

$$
H=1.4,\quad\frac{d(\theta_mU_e^{H+2})}{ds}=\frac{C_f}{2}U_e^{H+2},\quad C_{f,l}=\frac{0.664}{\sqrt{Re_s^*}},\quad C_{f,t}=\frac{0.0592}{(Re_s^*)^{1/5}},\quad\theta_{flat,t}=\frac{0.037s}{(Re_s^*)^{1/5}}.
$$

fixedShape用于复现固定形状因子积分；缺少连续上游时的局部平板启动和LT域外回退都必须记录来源。默认Head不是此固定H积分。平板回退不恢复真实不利压梯、横流、分离或激波相互作用适用性。

### 6.7 Reynolds比拟与有符号热流

$$
St=\frac{C_f}{2Pr^{2/3}},\qquad q_w=\rho_eU_eSt[h(T_{aw})-h(T_w)].
$$

旧calcSurfaceHeating用同物性的割线等效Cp把焓差写成换热系数乘温差；新相态v224直接保持焓差。密度在此为外缘ρ_e，不能与Kinney的参考ρ*混用。正、负、零热驱动均保留；T_w未知或焓域外不计算热流。Pr修正和摩擦关联属于工程估算，没有自动包含激波—边界层相互作用、热化学或催化影响。

### 6.8 Kinney型参考流长和附着面热流对照

令G=ρ_eμ_eU_e，I=∫G ds。背风加权长度S_w=I/G，迎风取几何准流线距离；项目在背风附近5度内smoothstep切换。当前源码明确使用线性参考温度：

$$
T^*=\frac{T_e+T_w}{2}+0.22(T_r-T_e),\quad \rho^*=\frac{p_e}{RT^*},\quad Re^*=\frac{\rho^*U_eS}{\mu^*}.
$$

$$
C_{f,l}=\frac{0.664}{\sqrt{Re^*/M_{f,l}}}\frac{T_e}{T^*},\quad C_{f,t}=\frac{0.37}{[\log_{10}(Re^*/M_{f,t})]^{2.584}}\frac{T_e}{T^*},\quad q_w=\rho^*U_e\frac{C_f}{2Pr^{2/3}}[h(T_r)-h(T_w)].
$$

默认Mangler因子M_f,l=3、M_f,t=2；它们是简化工程因子，不是已计算的一般三维流管变换。来源Kinney2007 II.I，印刷4至5页；原文使用CEA，本实现沿用NASA7/定Cp及常Pr。[R09]

只替换原工程链已接受侧壁，Thwaites/Head仍承担适用筛选。v104真实11880个侧壁面全部继承711个离散根的零积分入口假设；productionQualified=false、historyCompatible=false，新增面积0。两种方法面积加权热流差45.2129%表示模型分歧，不是误差改善。此公式与当前主链参考焓的差别必须保留，不在报告中擅自统一。

## 7 有实际剖面的层流方法

### 7.1 Levy–Lees局部相似坐标与边界条件

NASA TP-3271式6至7、17、25、54、56至61为项目局部相似剖面依据。r_t(s)为轴对称横向度量；球面r_t=R_n sin(s/R_n)。为避免与比焓混淆，此处把源文部分使用的h度量统一记r_t。

$$
\xi=\int_0^s\rho_e\mu_eU_er_t^2ds,\quad K=\frac{\sqrt{2\xi}}{\rho_eU_er_t},\quad\beta=\frac{2\xi\,dU_e/ds}{\rho_e\mu_eU_e^2r_t^2},\quad \frac{dp_e}{ds}+\rho_eU_e\frac{dU_e}{ds}=0.
$$

$$
F=u/U_e,\quad D=\rho_e/\rho=T/T_e,\quad l=\frac{\rho\mu}{\rho_e\mu_e},\quad y=K\int_0^\eta D\,d\eta.
$$

外缘同总焓、同EOS及沿程压力/速度梯度须先相容，不能同时给不一致p和U。边界为f(0)=0、F(0)=0、T(0)=T_w；远端F→1、T→T_e。beta要求0至1，因此本实现不处理任意不利压梯分离。局部相似省去沿程非相似项，不是全Navier–Stokes。[R10]

### 7.2 新相态变Pr五状态BVP的实际方程

令H=H_0−offset>0、g=[h(T)+U_e²F²/2−offset]/H，m=lF′，N为总焓扩散通量变量，撇号表示d/dη。新solvePhaseAwareLaminarProfile的实际五状态方程可以写成：

$$
\begin{aligned}f'&=F,\quad F'=m/l,\\ m'&=-fF'+\beta(F^2-D),\\ N&=(l/Pr)g'+(1-1/Pr)(U_e^2/H)Fm,\\ g'&=(Pr/l)[N-(1-1/Pr)(U_e^2/H)Fm],\\ T'&=\frac{Pr}{lc_p}(HN-U_e^2Fm),\quad N'=-fg'.\end{aligned}
$$

Pr=c_p(T)μ(T)/k(T)随温度变化，N中的耗散/变Pr项不能省略。程序以z=ln[(T−60)/(2000−T)]代替温度推进，达到端点饱和即拒绝，不裁温；终解各p/T仍另查联合域。

$$
q_w=\frac{\mu_eH}{K}N(0),\qquad\tau_w=\frac{\mu_eU_e}{K}m(0),\qquad C_f=\frac{2\tau_w}{\rho_eU_e^2}\quad(U_e>0).
$$

这里y由壁面向流体，q=k∂T/∂y|w，正号表示热进入壁。驻点不以U_e=0硬除求Cf；取有限极限K=sqrt[ν_e/(2a)]、a=(dU_e/ds)_0>0、beta=1/2，并检查p_e,ss=−ρ_ea²。球鼻前缀必须核对这项二阶曲率，点值U=dp/ds=0不足以证明兼容。

### 7.3 物理厚度与有限外缘

$$
\begin{aligned}\delta&=K\int_0^{\eta_e}D\,d\eta,\quad \delta_m=K\int_0^{\eta_e}F\,d\eta,\\ \delta^*&=\delta-\delta_m,\quad \theta_m=K\int_0^{\eta_e}F(1-F)\,d\eta.\end{aligned}
$$

这些分别对应实际密度加权的质量/位移/动量厚度，来自真实剖面积分；不是把Thwaites参考厚度乘一个比例。有限η_e要求速度、静温、密度、总焓、剪切及热梯度都接近外缘；数值域η_max不是物理边缘。

旧solveAxisymmetricLaminarMassProfile使用NASA7/定Cp、Sutherland及给定Pr；新相态使用同Cp和线性μ/k、变Pr，二者应分别标明。能量offset只是数值变量平移，输出保留物理焓零点。

### 7.4 旧零攻角球鼻热边界与新偏置球鼻候选

旧solveSphericalNoseLaminarHeating要求零合成攻角、球面前缀、附着层流、Kn<0.01、Re_s低于配置转捩阈值、δ/R_n≤0.05、80至2500 K、最多16种给定壁温。默认2度角表；恢复温度可由同一BVP的q(T_w)=0根取得，失败根不以总温代替。热壁小幅速度过冲可以是附着支，不能一概用F≤1删掉。

新v227/v228采用偏置拟合球鼻旋转近似，完整短弧r_t=Rsinθ、32点Gauss积分ξ、2度相邻角站插值；代码要求T_w严格在60与2000 K之间、beta在0至1、Re_s<5×10^5、δ/R<0.05，及剖面/联合域门限。实际原箭体31角站/1306面是候选，不是全箭轴对称或真实三维非相似解。新T_aw仍为sqrt(Pr)恢复估计，不是该BVP求得的绝热壁温。[R05、R10]

v227曾误取外缘m导致壁剪错误；当前取保存壁列Y(3,1)，v228实际30度站回归τ=121.5200727 Pa、q=257888.6681 W/m²，与独立温度BVP最大相对差3.305×10^−9。此为方程数值复算；未重算全部31站，也没有新增实验或真实区域资格。

### 7.5 v61非相似层流

$$
\begin{aligned}\partial_s(r_t\rho u)+\partial_y(r_t\rho v)&=0,\\ \rho u\partial_su+\rho v\partial_yu+\partial_sp_e-\partial_y(\mu\partial_yu)&=0,\\ \rho u\partial_sH+\rho v\partial_yH-\partial_y(k\partial_yT+\mu u\partial_yu)&=0,\quad H=h+u^2/2.\end{aligned}
$$

regionalFlow/nonsimilar采用物理s/y薄层方程，连续式由壁面无吹吸积分v，不额外强置外缘v=0。给定真实入口剖面、整段T_w历史，同物理法向网格；BDF1/BDF2及隐式Newton推进。要求零攻角轴对称、连续附着层流；不处理三维横流、回流和薄层法向曲率高阶项。

已冻结为条件研究，不接默认生产。40度局部相似入口至60度的保存段虽有实际数值与面元映射，但独立N121首站动量积分缺陷约−4.973462%；小Newton残差和末站热流接近不能代替积分账。physicalSolution=false、eligibleForSurfaceCoverage=false，新生产面积0。[R12]

## 8 底压 底热与尾锥

### 8.1 SAND93基本量和底压

Lamb与Oberkampf的SAND93-0280是项目已进行原页和公式复算的底关联来源。e为肩部当地外缘，s为鼻端到肩部母线长，θ_c数值以度输入。

$$
Re_{e,s}=\frac{\rho_eU_es}{\mu_e},\quad C_e^4=\left[1+\frac{2}{(\gamma-1)M_e^2}\right]^{-2},\quad A=p_b/p_\infty.
$$

$$
\begin{aligned}A_{cone,l}&=3.05\left(\frac{M_e^2}{\sqrt{Re_{e,s}}}\right)^{0.6}/C_e^4,\\ \frac{p_{b,cyl,t}}{p_e}&=0.05+\frac{0.967}{1+(\gamma-1)M_e^2/2},\\ A_{cone,t}&=\left(\frac{p_e}{p_\infty}\right)^2\left[0.025+\frac{0.906}{1+(\gamma-1)M_e^2/2}\right]^J,\quad J=\frac{1.7}{\ln(21/\theta_c)}.\end{aligned}
$$

分别对应报告式8b、10、11。肩部M_e和p_e不能用任意自由流替代。底压近似均匀由低速回流静压假设支持，不解析底唇局部细节。[R13]

### 8.2 代表底热

$$
St_{b,l}=10^{-4}(AC_e^4)^{1.8},\quad B_b=\sqrt{A(p_\infty/p_e)^2},\quad St_{b,t}=8.3\times10^{-4}B_b^{1.85},\quad q_{b,rep}=St_b\rho_\infty U_\infty(H_\infty-h_w).
$$

对应式17、18；平方根作用于整个A(p∞/p_e)²。St以自由流质量通量和总焓差定义，允许负/零热驱动。这里q_b,rep是代表量，资料未建立其为全底平均、中心值或某径向基准的完整映射，不能赋给每个底面。

| 分支 | 当前独立接口门限 | 仍需保留的来源限制 |
| --- | --- | --- |
| 圆锥层流底压及代表底热 | 半角9至10度；M∞为2至20；M_e>1；Re_e,s≥10^5 | 主要9度锥数据；低Re不可靠 |
| 圆锥湍流底压及代表底热 | 半角6至12度；M∞为2至20；M_e>1 | 主要单一9度钝锥；其它角度为受限估算 |
| 圆柱湍流底压 | M_e为1.2至15 | 不含对应圆柱底热闭合 |
| 圆柱层流 | 拒绝 | 当前无闭合，不套圆锥公式 |

共同要求零攻角/侧滑、孤立平圆底、无尾撑/喷流、鼻钝度R_n/R_b≤0.3、定Cp空气γ=1.4、同总焓/肩部状态、Kn<0.01。80至2500 K是代码筛选，不保证定Cp真实精度。厚湍流边界层被源研究排除，而独立接口没有自动证明肩部薄层，调用者须提供相应依据。

### 8.3 真实网格底压的更窄接入

solveMeshBasePressure并非把全部独立分支都接入网格。当前只接受同轴尖头直圆锥、9至10度半角、零攻角、闭合单一二流形及完整平圆底；显式无尾撑/无喷流，肩部Taylor–Maccoll弱解，经验转捩仍层流，Re≥10^5，同剖面物理δ/R_b≤0.05。还核对真实三角面积、底圈圆拟合和同一气体/总态。

压力模型可映射到通过审核的真实底三角，热流和外缘/恢复温度仍未知。默认主链底压是prescribed，只有显式sand93Laminar才调用；不支持目标箭体、球鼻锥柱、尾锥或α20度工况。δ/R_b≤0.05是项目保守门限。[R14]

### 8.4 底面径向幂函数的资料边界

$$
q_b(r)\propto(1-r/R_b)^{1.2},\qquad 0\le r/R_b\le0.8.
$$

项目在NASA NTRS19810011647印刷155页核对到该趋势，引用Rabinowicz1958，不是已读到Bulmer1974的完整模型。仅为比例式，振幅与SAND93代表量对应关系未知；适用半径到0.8只占圆底64%面积，外侧36%及底唇没有解，不能外推到1并补零。

Rabinowicz1957报告明确径向测位主要为中心和r/R_b=2/3，保留尖杆促湍、侧撑与激波管条件。两测位不能唯一确定全圆底函数。Bulmer DOI10.2514/3.62179及报告SLA-74-230有书目，但原公式全文在既有记录中未取得。[R15]

### 8.5 声速流线尾锥压力候选

$$
\nu(M_s)=\theta_b,\quad\nu(1)=0,\qquad\frac{p_{candidate}}{p_c}=\left[\frac{1+(\gamma-1)/2}{1+(\gamma-1)M_s^2/2}\right]^{\gamma/(\gamma-1)}.
$$

estimateSonicBoattailPressure复算NASA TN D-2761的声速流线近似。p_c是收缩前当地柱面静压，M_s是边界层内声速流线膨胀代理，不是给真实外缘强置M=1。γ=1.4，直尾锥角0至21度；18/21度已有分离迹象，公式不解分离泡或底区。

数学复算及15个原图实验符号的条件压力比较已测；未接默认生产，不给q、T_e、T_aw或固体T_w。当前约31度真实收缩不在其代码范围内，不能只替换角度沿用。[R19]

## 9 熵来源 流管与共同激波层研究

### 9.1 压力驱动的曲面流向和横向度量

$$
\mathcal D_s\boldsymbol t=-\frac{(I-\boldsymbol t\boldsymbol t^T)\nabla_\Sigma p}{\rho U^2},\qquad h_\beta=(\boldsymbol n\times\boldsymbol t)\cdot\frac{\partial\boldsymbol r}{\partial\beta}.
$$

v94至v96在给定光滑曲面与相容连续压力及其两向导数下，积分流向并用相邻标签恢复横向尺度。h_β不是轨迹间距的模，也不是Dijkstra距离。每个入口标签保留同一H_0和S_0，真正横向导数由相邻标签重建。

$$
\nabla_\Sigma T=\frac{T}{c_p}\left(\nabla_\Sigma S_0+R\frac{\nabla_\Sigma p}{p}\right),\qquad\nabla_\Sigma U=\frac{\nabla_\Sigma H_0-c_p\nabla_\Sigma T}{U}.
$$

用于光滑、非驻点、附着补片。输入压力是给定或其它已解模块提供；本模块不求未知压力，不证明完整三维连续式，不解边界层/q。聚焦、邻线停止、激波或分离后结束有效前缀，不能靠单条辅助射线填值。v95显式使用连续NASA7积分约定，新相态接入仍需相容转换。[R11]

### 9.2 CR112180质量吞入与熵匹配

源报告式65至69以近似激波和边界层质量相等选择实际穿越激波的流线。Γ为激波与来流夹角，h_β为横向度量，δ_m=δ−δ*：

$$
\begin{aligned}\frac{d\sin\Gamma}{dS}&=\frac{(p_b-p_{sh})h_\beta}{U_\infty\Psi_{sh}},\quad \Psi_{sh}\approx\rho_\infty U_\infty\int\sin\Gamma\,h_\beta dS,\\ \Psi_{BL}&=\rho_eU_eh_\beta\delta_m,\qquad\Psi_{sh}(\Gamma_{cross})=\Psi_{BL}.\end{aligned}
$$

源文p_b是激波形状推进的辅助压力，不是可直接显示的底压或壁压。本文式中S为源激波坐标，不能与比熵符号混同。源文δ_m=5.27θ_m层流、δ_m=θ_m(10.3−H)湍流对应其特定变换定义，不能直接乘项目Thwaites参考θ。

项目recoverShockEntropyEdgeState、recoverEntropyLayerEdgeState等已提供有来源输入下的热力恢复与匹配核验，轴对称球锥条件研究仍有限制。v52鼻端耦合曾导致压力曲率与局部相似应变率约1.9951%失配，不能以s=0点值零残差通过。当前不构成真实20度后体外缘生产器。[R20]

### 9.3 v97至v100共同激波层原型

v97共同曲面向量跃迁、v98受限曲面薄层输运、v99一般凸图投影质量映射，均是独立研究。一般图z=f(q)的投影质量系数a_m满足：

$$
\nabla_q\cdot(a_m\nabla f)=1,\qquad\frac{da_m}{d\tau}=1-\operatorname{tr}(\nabla_q^2f)a_m,\quad\frac{dq}{d\tau}=\nabla f,\quad\frac{d\ln\Psi}{d\tau}=1/a_m.
$$

它选择同一流线入口熵，不能替代壁面匹配。v99真实45列固定形变体形误差约170至188微米，未过10微米原门槛；显式一阶压力轴线端点约低0.42%至0.47%，完整能量恢复出现约55至61 m/s虚假余速。不能以继续加曲面自由度或裁速度消除。

v100转向轴线正则的法向动量/质量耦合。向内法向n_i、切向s_t下，曲率和相邻列项需要保留：

$$
\frac{dp}{dn_i}=-\rho\left(\frac{\kappa u^2}{h_g}+\frac{u}{h_g}\left.\frac{\partial v}{\partial s_t}\right|_{n_i}+v\frac{\partial v}{\partial n_i}\right),\qquad h_g=1-n_i\kappa.
$$

当前v100仍为局部轴对称条件研究，高阶独立积分、缩放及外侧域有未通过项。v101至v103激波拟合/ALE接口也属于历史研究资产；按当前目标不继续自研通用流场内核。它们不授予真实壁压、外缘或热面积。[R21]

## 10 已收集但尚不能直接接入的关联

### 10.1 Lees与Kemp球鼻分布

Tauber NASA TM-101077印刷7至8页式42至43记载以下带变换参数的分布；这里只转录项目已保存的草稿，不将其当作已实现函数：

$$
\frac{q_w}{q_{ws}}=\Gamma\sqrt{\frac{\rho_w\mu_w U'_e/(\beta t_e)}{\rho_{ws}\mu_{ws}a/(\beta_s t_s)}}.
$$

ρ_wμ_w为当地壁态而非外缘；t_e、t_s、Γ是原热力变换参数，Γ=1对应该来源的Lees对照。配套beta还依赖含Γ²的沿程积分。全部变换定义尚未与现有内核逐项闭合，因此不省略参数写成通用高焓式。

配套压力梯度参数的项目草稿如下，j=1表示轴对称，L为原文参考长度。此处t_s及Γ严格保留源定义，不与流管标签或激波穿越角混用：

$$
\beta=\frac{2t_s U'_e}{\rho_w\mu_wU_e^2(r_t/L)^{2j}\Gamma^2}\int_0^s\rho_w\mu_wU_e(r_t/L)^{2j}\Gamma^2\,ds.
$$

来源为薄、附着、轴对称零攻角冷壁层流。不能用一条冷壁q/q_s分布反造恢复温度或热壁BVP。Tauber图5的24个圆点草稿涉及Kemp等解离空气实验，源工况和催化性仍未齐，不作为当前固定组分M5/Tw300实验误差。[R22]

### 10.2 激波干扰和再附热峰

Tauber式60及63至64在项目来源审查中已记录：

$$
\frac{C_{h,peak}}{C_{h,up}}=\left(\frac{p_{peak}}{p_{up}}\right)^{0.85},\qquad C_h=\frac{q_w}{h_{aw}-h_w}.
$$

$$
St_{max,l/l}=0.664(Re'_x)^{-1/2}(Pr')^{-0.667},\qquad St_{max,l/t\;or\;t/t}=0.0576(Re'_x)^{-1/5}(Pr')^{-0.667}.
$$

0.85是来源给出的特定拟合，审查记录其它实验约0.73至0.85，受流态和壁温影响；不能当作任意三维热峰定律。带撇号的Re、Pr使用原局部参考定义，不能直接换自由流Re。需要上游已解边界层/流态、实际压力跳、参考焓、再附位置；当前没有完整位置/长度闭合。[R22]

这些只给峰值或增强，不给热峰宽度、分离长度、再附位置与全区分布；不能涂满所有分离面元，也不能用压力峰唯一求出整区热场。

### 10.3 腔体和层流分离平均关联

$$
St_{ave}=0.48Re_L^{-0.4}(L/H_c)^{0.2},\qquad St=2.07/Re_\delta.
$$

前者是Tauber式65的平面湍流腔体平均，H_c为腔深；后者式66来自Lamb1980平面超声速层流分离/腔体关联。流态、几何、Re和热参考必须按原文恢复，不能用于轴对称圆底或任意回流区。当前未实现完整区域模块。

Chapman NACA TN-3792的能量积分只能给总体平均热率，原文明确不能给局部分布。Meng NASA TN-D-7508采用近尾区NS与其它区边界层、Telenin条带/多重打靶，是特定二维数值方法，需热边界层入口及下游近尾迹剖面；不是一条可直接粘贴的底热相关式。[R22]

Tauber式67另记录圆锥湍流代表底热的比例式(Nu/Pr)_b/(Nu/Pr)_c=35.5(Re_b/Re_c)^2.2。b与c及Re、Nu的原始参考定义必须随源资料恢复；它同样不是径向分布，本项目未把它与SAND93拼接成全底模型。

### 10.4 短突起增强

Mazaheri与Wood2008的短圆钝突起关联针对k_p/δ<0.3，输出相对光滑基体热流的唇部峰值增强；此处k_p是突起高度，不是热导率。

$$
A_q=\frac{q_{lip,peak}}{q_{smooth}},\qquad k_p/\delta<0.3.
$$

项目已核对该定义与条件，但没有完成可用于任意构型的完整系数和输入闭合，本报告不编填增强函数。需要已解基体q、边界层厚度与焓剖面。不能套锐鳍、悬挑、高凸台、突起群，也不能赋整个突起根部/侧面/背风面热流。原运行存在约67.78 K低温条件，旧80 K物性链不能裁到80 K后称为同工况。[R23]

### 10.5 尚未完整核实的扩展方法

Zoby/Moss/Sutton完整热流闭合已有书目和被引用关系，但既有资料记录原全文访问未完成，不能从二手OCR实现缺符号积分。CBAero/MINIVER公开说明可作路线比较，不等于本项目复刻其完整算法或获得同样验证。

BL3D的可压缩三维附着层流、增强PNS、黏无黏相互作用均是候选。BL3D原应用为后掠翼，未完成本项目逐式代码审查；常规PNS空间推进受推进方向超声速约束，增强全局迭代可处理部分上游影响但代价和边界须另验收。不能仅更换方法名称宣称解决真实三维分离。[R23]

Fay–Riddell、Detra等常见驻点名称不作为本报告“项目已用方法”补入：当前核查的主要实现是Sutton–Graves及局部相似BVP。若以后采用，应建立各自源公式、热化学、催化及壁面条件合同，不能套用别的驻点法名称。

## 11 当前混合路线与并发候选

### 11.1 Euler和局部黏性求解的分工

全局Euler应给连续p/T/ρ/三维U、H_0/S及激波/流向，供相容外缘和流管审核；Euler没有黏性壁热或真实边界层剖面。局部NS/RANS应求强相互作用区p/q/壁剪、分离再附与热峰，需可信上游剖面、特征外边界、周向/下游缓冲和必要反馈。

$$
\partial_t\boldsymbol Q+\nabla\cdot\boldsymbol F_c=\nabla\cdot\boldsymbol F_v,\quad \boldsymbol Q=(\rho,\rho\boldsymbol U,\rho E),\qquad q_w=k\left.\frac{\partial T}{\partial y}\right|_w.
$$

Euler取F_v=0；它的壁热通量0是无黏方程省略项，不能记录为壁热“已求零值”。RANS湍流热导、模型Pr_t与近壁处理必须按固定版本声明；本报告不把SA-neg、SST、Wilcox相互替代，也不把单个模型收敛当实验真值。

v230相态SU2 Euler完成30更新后异常，仅0/1完整场，当时std::terminate的原始原因未知。随后v231同输入捕获到RAH_PHASE_ENERGY_ENTROPY_DOMAIN及ID17337量热内能域外状态；将具名相态域拒绝接回原生Solution_Old恢复路径，合法/拒绝控制通过。该拒绝接口修复不裁温、不放物性域；完整异常调用栈未取得，不泛称排除所有其它潜在故障。

编制末期v232主线记录真实800更新正常完成、保存18场，末50步原面压力变化0.000904%已稳定；但五式−8及净质量/物理能量1e−8门槛未过，末预处理仍有8点原生回退，空间/外域与热资格未过。更新后壁压与更新所用阶段RHS/cache的5.67174%差属于跨阶段比较，不能当成独立同阶段动量守恒误差；v230历史10.2848%也须保留这一取样边界。最低温尾部点不能自动认作回退节点。新外缘/热资格仍0，原面13字段/217272值用于研究展示。v222完整黏性场200更新同样未闭合稳态与完整收支。本报告未运行这些任务，只引用主线既有完成记录。[R16]

### 11.2 工程层与CFD入口必须相容

输入包至少绑定原面、几何/工况/源版本、气体及输运、物理能量参考、p/T/三分量U、H_0/S、真实起点/前驱/弧长、相邻流线尺度及适用掩码。边界层入口还需实际u(y)、T(y)、密度/质量厚度、湍流身份及必要变量。

壁邻流体状态不能自动叫边界层外缘；单个积分θ/H不能唯一恢复全部剖面；把区域边界固定p∞或堵住回流不能代替外部反馈。全局高Mach不保证每个边界法向/近壁都超声速。强三维、收缩、后体和底面仍需专用黏性反馈；原3128面/0.317970330711786 m²前体是下一目标，不是已恢复面积。[R07、R16、R17]

### 11.3 本次工作区可见的独立候选

本次git检查发现harness_tasks及engineering54相关并发资产。定向阅读显示，三候选calcBoundaryLayerState、calcTransition、marchHeadBoundaryLayerSegment有独立集成恢复报告，27/27项软件记录通过、补丁仅检查而未应用。这不能覆盖根目录当前算法，也不授真实箭体物理资格。

Head来源审查明确区分原始卷吸思想、当前数值拟合、工程连续化和程序门限；具体幂律原始页式仍未核实。该信息已反映在第6.5节。候选路径、测试状态和准确SHA随本报告来源清单单列，保持“已收集或已测候选”身份，未执行合并、MATLAB或新计算。[R08、R24]

## 12 验证口径与区域选择建议

### 12.1 当前数字应如何使用

| 记录 | 当前可引用的事实 | 不能推出的结论 |
| --- | --- | --- |
| 旧工程q计算面积66.5626327424% | 有受限公式输出的面积，缺热约0.6569202998 m² | 不是已获独立物理认证的覆盖率 |
| 旧工程p计算50.2710271898% | 与给定背风/底压分开统计 | 不是全压场已解 |
| 同面p/q计算交集约43.0999% | 同一面同时有计算值的统计 | 不代表外缘和所有温度相容 |
| v227/v228鼻部1306面 | 0.1790893%全表面积的局部工程热候选 | 不与旧链重复相加，不授新热资格 |
| v232原面13字段 | 更新50/400/800及保存18场的研究证据 | 压力稳定不等于稳态或外缘合格，Euler不提供q |
| 公式复算误差很小 | 指定模型/数据下的软件与方程一致性 | 不等于实验误差或真实三维精度 |

### 12.2 面积 功率与资格的统一口径

$$
A_{computed,f}=\sum_{i\in C_f}A_i,\qquad\eta_f=\frac{A_{computed,f}}{\sum_iA_i},\qquad Q_{known}=\sum_{i\in C_q}q_iA_i.
$$

集合C_f按变量独立确定，给定量不算计算；研究候选与合格结果分开。所有原面均保留在分母；重叠区域每面每变量一个主结果，不能重复计面积/功率。Q_known只包含有效面，不能称整箭总热功率。完整合格应另外审来源、同物性、守恒、网格/外域、入口/接口及公开实验，而非只检查finite。

### 12.3 按区域选择方法

| 区域或任务 | 可复用算法 | 进入可信结果还需具备 |
| --- | --- | --- |
| 核验球鼻驻点 | Sutton–Graves或同相态局部相似BVP | 半径/偏置驻点、同气体、真实工况与热边界验证 |
| 球鼻前缀分布 | 驻点熵外缘与局部相似 | 完整短弧、压力曲率、薄层/附着、非相似/三维误差评估 |
| 平顺附着前体 | 连续全局外缘＋工程积分/剖面 | 可信壁邻—外缘、入口熵/焓、流管度量与转捩声明 |
| 强收缩与再附 | 成熟局部NS/RANS＋必要反馈 | 上游热剖面、外域/回流、近壁、稳态或统计与热验证 |
| 平圆底压 | 受限SAND93映射或黏性解 | 真正匹配前体/肩部、操作声明与薄层条件 |
| 全底热流 | 专用回流/剪切层/壁热求解 | 全空间分布及同工况测点、底唇/尾迹和守恒 |
| 小圆钝突起峰值筛查 | 来源匹配的增强关联 | 已解光滑基体与剖面，不把峰值扩成全表面 |
| 一般翼片及突起群 | 有依据的局部黏性区域 | 真实外边界、遮挡/干扰、几何与网格特征保留 |

报告后的实际研发优先级仍按PROJECT_GOALS及最新v232/NEXT_DECISION执行；本报告不批准额外数值运行，不改变任何门限或原失败身份。新的方法应先有完整来源与相容输入，再按真实区域逐个验收并接入科研工作台。

### 12.4 已收集公开基准与匹配条件

| 资料 | 可支持的算法对照 | 必须匹配及当前证据边界 |
| --- | --- | --- |
| NACA1135表II | 定Cp正激波等解析公式 | 独立转录复算属于解析校核，不是热流实验 |
| NASA TM-101075及TMR轴对称Mach7压缩拐角 | 黏性压力、热流、入口与分离相互作用 | 原实验当地M7.05与调整计算入口7.11区分；Tw311 K；不是20度攻角箭体 |
| NASA TN D-962 | M6尖/钝15度锥压力和热流 | 主分布为0及45/60/70/80/90度；不能虚构20度测点，薄壁瞬态热流须匹配热边界 |
| Penland NACA TN-3861 | M6.86圆柱横流压力 | 实际14.9/30.2/44.9/60.4/75.7/90度；无20度实测，q∞与横流动压分母不同 |
| Rabinowicz1957 HRP41 | 半球柱底面中心与2/3半径热流 | 促湍杆、侧撑、激波管Mach身份及试验时间不能省略；不是全半径曲线 |
| NASA TN D-1372独立来源笔记 | 超声速湍流平板及附着突起增强 | 工作区笔记记载已取得官方OCR，M2.65/3.51/4.44、三档δ；PDF未完整取得，逐点表列对应未恢复 |
| Run14理想N2控制 | 同几何层流压力与热峰链路 | 原控制未收敛；N2气体与Tw/入口匹配，不能转授空气M5/三维湍流资格 |

TN D-1372在本次工作区已有2026年10月30日来源笔记，更新了早期路线文档中“全文访问失败”的检索状态。该笔记报告20构型、Re/ft约1.3×10^6至4.7×10^6、δ为0.70/1.50/6.00英寸；当前只有OCR条件登记，未可靠恢复逐点热流/壁温列，也没有本项目实验对比。h/h_0是同测位有/无突起换热系数比，不能把不同阶段存在约15%至20%偏移的绝对h混用。这里引用独立来源记录，不称本轮又审核了原始扫描页或完成物理验收。[R25]

## 13 来源与代码索引

下列索引提供项目文献入口、已保存的审查位置和实际程序。本次没有重新验证全部远程URL可达性；准确本地文件SHA、大小和读取时身份见同目录reports/20261112_algorithm_inventory/source_manifest.json。公式以当前源码/明确来源为准，历史文档中的旧下一步不覆盖v232状态。独立工作区来源笔记的小型原文副本保存在source_notes，仅作报告归属证据，不复制其候选程序或工程包。

| 索引 | 一手来源或项目记录 | 代码和证据定位 |
| --- | --- | --- |
| R01 | 五步骤方法边界及旧工程链 | WORKFLOW_METHOD_SCOPE.md；airThermalProperties.m；calcSurfaceHeating.m |
| R02 | Lemmon2000标准Air；Lemmon与Jacobsen2004输运 | HYBRID_PHASE_AWARE_AIR_V206.md；phaseAwareAirState.m；data/thermo/phase_aware_dry_air_v1.json |
| R03 | 同物性激波与膨胀结构及v225 | HYBRID_PHASE_WAVES_V225.md；solvePhaseAwareWave.m；marchPhaseAwareWavePath.m |
| R04 | NASA SP-8062第4.2.1节 印刷30页 | SPHERICAL_NOSE_EDGE.md；HYBRID_PHASE_NOSE_EDGE_V226.md |
| R05 | 新相态局部相似鼻部及壁列回归 | HYBRID_PHASE_NOSE_HEATING_V227.md；HYBRID_PHASE_CONTINUATION_V228.md |
| R06 | Riley NASA TM-107838附录H.6/H.7 印刷117页 | HYBRID_REFERENCE_ENTHALPY_V190.md；airReferenceTemperature.m；phaseAwareAirReferenceState.m |
| R07 | 同相态附着层流热链 | HYBRID_PHASE_LAMINAR_CHAIN_V224.md；calcPhaseAwareLaminarHeating.m |
| R08 | Head ARC R&M3152及LT原始拟合页待核 | marchHeadBoundaryLayerSegment.m；harness_tasks/head_closure_review_20261108/SOURCE_REVIEW.md |
| R09 | Kinney2007 II.I 印刷4至5页 | KINNEY_ACREAGE_V104.md；calcWeightedReferenceLength.m；calcKinneyAcreageHeating.m |
| R10 | NASA TP-3271式6至7 17 25 54 56至61 | SPHERICAL_NOSE_HEATING.md；solveAxisymmetricLaminarMassProfile.m；solvePhaseAwareLaminarProfile.m |
| R11 | Riley与DeJarnette NTRS20040121205印刷4页式13至14 | PRESSURE_DRIVEN_SURFACE_FLOW_V94.md；THERMAL_PRESSURE_SURFACE_FLOW_V95.md；LABELLED_PRESSURE_SURFACE_FLOW_V96.md |
| R12 | TP-3271薄层方程及项目v61冻结审核 | NONSIMILAR_BOUNDARY_LAYER_DRAFT.md；V61_ENERGY_AUDIT.md；regionalFlow/nonsimilar/ |
| R13 | Lamb与Oberkampf SAND93-0280式8b 10 11 17 18 | AXISYMMETRIC_BASE_CORRELATION.md；solveAxisymmetricBaseCorrelation.m |
| R14 | 真实底面及肩部薄层审核 | MESH_BASE_PRESSURE.md；solveMeshBasePressure.m；solveSharpConeBaseFlow.m |
| R15 | Rabinowicz1957 HRP41；NASA19810011647印刷155页；Bulmer1974书目 | BASE_RADIAL_HEATING_RESEARCH.md |
| R16 | 相态SU2固定实现及真实运行和域拒绝修复 | HYBRID_SU2_PHASE_AIR_V220.md；HYBRID_PHASE_EULER_NATIVE_V229.md；HYBRID_PHASE_EULER_ROCKET_V230.md；HYBRID_PHASE_PRIMITIVE_RECOVERY_V231.md；HYBRID_PHASE_EULER_DEVELOPMENT_V232.md；validation/v232/NEXT_DECISION.md |
| R17 | 混合路线和旧膨胀根缺口 | PROJECT_GOALS.md；HYBRID_ROUTE_SELECTION_V189.md |
| R18 | Billig1967 DOI10.2514/3.28969 | estimateBowShockSurface.m |
| R19 | NASA TN D-2761印刷7至9页及图9至10 | SONIC_BOATTAIL_PRESSURE.md；NASA_BOATTAIL_PRESSURE_COMPARISON.md；estimateSonicBoattailPressure.m |
| R20 | NASA CR-112180式65至74；NASA TN D-2970 | NOSE_AFTERBODY_MATCHING_RESEARCH.md；ENTROPY_MASS_COUPLING.md |
| R21 | 共同激波层与轴线正则研究 | SHOCK_LAYER_TRANSPORT_V98_DRAFT.md；GENERAL_SHOCK_SURFACE_MASS_V99.md；COUPLED_SHOCK_PATCH_V100.md |
| R22 | Tauber NASA TM-101077式42至44 60 63至67；Chapman TN-3792；Meng TN-D-7508 | REGIONAL_HEATING_SOURCE_REVIEW.md；tests/data/experimental_hemisphere/ |
| R23 | NASA TM-107838；CR-4531；短突起NTRS20080024047；TMR及TN-3861 | REGIONAL_ALGORITHM_ROADMAP_2026-10-17.md；REGIONAL_ALGORITHM_RECHECK_V92.md；ENGINEERING_SOURCE_TRIAGE_V93.md；CYLINDER_CROSSFLOW_SOURCE.md |
| R24 | 工作区独立候选恢复报告 未合入根目录 | harness_tasks/integration_recovery_20261110/REPORT.md及candidates/三文件 |
| R25 | Burbank等 NASA TN D-1372 1962 独立OCR来源记录 | docs/research/opencode_protuberance_d1372.md；d1372_verified_conditions.csv；本轮未复核原扫描页 |

原始入口如下。除DOI和条目索引外，链接文本中的页式身份均对应上述本地记录，不表示本轮又完成一次原页验证。

- Lemmon2000标准Air：https://trc.nist.gov/refprop/Documents/Air.pdf
- Lemmon与Jacobsen2004输运：https://trc.nist.gov/refprop/FAQ/NAO.PDF
- NASA SP-8062：https://ntrs.nasa.gov/api/citations/19710021703/downloads/19710021703.pdf
- NASA TP-3271：https://ntrs.nasa.gov/api/citations/19930006245/downloads/19930006245.pdf
- Riley NASA TM-107838：https://ntrs.nasa.gov/citations/19920016036
- Riley与DeJarnette曲面流向：https://ntrs.nasa.gov/api/citations/20040121205/downloads/20040121205.pdf
- Kinney2007：https://ntrs.nasa.gov/citations/20070017947
- SAND93-0280：https://doi.org/10.2172/10111104
- CR-112180：https://ntrs.nasa.gov/api/citations/19730005231/downloads/19730005231.pdf
- Tauber1989：https://ntrs.nasa.gov/api/citations/19890017745/downloads/19890017745.pdf
- 短圆钝突起2008：https://ntrs.nasa.gov/citations/20080024047
- Chapman层流分离：https://ntrs.nasa.gov/api/citations/19930084514/downloads/19930084514.pdf
- Meng尾迹数值方法：https://ntrs.nasa.gov/api/citations/19740005592/downloads/19740005592.pdf
- Rabinowicz原报告：https://authors.library.caltech.edu/records/df7cq-am353
- Bulmer径向底热书目：https://doi.org/10.2514/3.62179
- Zoby完整热流闭合书目：https://doi.org/10.2514/3.57788
- NASA TN D-1372资料：https://ntrs.nasa.gov/citations/19630001117

## 附录 已有热响应与后续流态的边界

项目还保留壁体导热、辐射/背面边界及历史累计功能，它们是给定气动热边界后的响应算法，不是恢复当前未解外流/热流的方法。基本热平衡为：

$$
\rho_sc_s(T)\frac{\partial T}{\partial t}=\frac{\partial}{\partial x_s}\left(k_s(T)\frac{\partial T}{\partial x_s}\right),\quad q_{net}=q_{conv}-\epsilon\sigma(T_w^4-T_{sur}^4)-q_{back},\quad E_A=\int q_w\,dt.
$$

具体材料分层、初温、背面边界、辐射配置及其验证须分别引用calcWallThermalResponse、solveLayeredWallConduction、completeWallEnergyBalance等模块，不能把给定Tw改称该方程的解。本报告不展开其时间离散；累计热量、长时程响应与轨迹仍按目标后置。

滑移/温跳、过渡流、自由分子流、化学反应、催化、辐射流场和烧蚀未因已有Kn或物性筛选而实现。任何跨流态拼接都需单独来源与连续/稀薄极限验证，不从当前连续流关联外推。
