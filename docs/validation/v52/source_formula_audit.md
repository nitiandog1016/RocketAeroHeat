# CR-112180 球锥激波质量表：原页核对与鼻端闭合

日期：2026-10-12。性质：只读文献研究和从原方程推导，未修改生产程序，未验证真实球锥激波/热流。源文献：F. R. DeJarnette, NASA CR-112180 (October 1972), *Calculation of Heat Transfer on Shuttle-Type Configurations Including the Effects of Variable Entropy at Boundary Layer Edge*。本地 PDF：`C:/Users/bo yunpeng/Documents/Codex/2026-10-12/nose_matching_research/NASA_CR_112180.pdf`。PDF 共99页；以下 PDF 页码为从1开始的文件页码。

## 原页证据和坐标合同

- 印刷页13 / PDF18，式14：真正的激波面质量坐标为 `Psi_sh = rho_inf U_inf integral sin(Gamma) h_sh dS_sh`，单位是每单位周向坐标的质量流率。轴对称情况下 `h_sh=r_sh`、`sin(Gamma)dS_sh=dr_sh`，因此准确地有 `Psi_sh = rho_inf U_inf r_sh^2/2`。这里的 `r_sh` 不是体表半径。
- 印刷页24 / PDF29，式64：先明确替换激波坐标度量 `(h_sh,S_sh)` 为体表流线度量 `(h,S)`；式65 是
  `d(sin Gamma)/dS = (p_b-p_sh) h / (U_inf Psi_sh)`。
  原符号带横线的 Gamma 是激波与来流速度的夹角，鼻端为90度；不能当作壁面转角或项目攻角。
- 同页，式66 原印刷式清楚给出近似符号：
  `Psi_sh = rho_inf U_inf integral sin Gamma h_sh dS_sh ~= rho_inf U_inf integral sin Gamma h dS`。
  因此在体表积分所得的是 **body-associated mass-angle approximation**，不是三维激波节点坐标。这一近似足以作为原文熵匹配输入；禁止把 `h=r_body` 又解释成实际激波半径。
- 印刷页25 / PDF30，式67：
  `p_b = p_inf + (p_e-p_inf)(p2-p_inf)/(p_stag-p_inf)`。
  段末明确写 `p_b is the pressure to be used in eq. (65) only`。p2 是正激波后静压；p_stag 是正激波后的驻点压力；p_e 是物理面压力模型。p_b 为积分辅助压力，不能冒充 p_e、表面压力或 UI 物理压力。

以上页已渲染并视觉核对：`shock_start-29.png`、`shock_start-30.png`。

## 不需额外曲率输入的零攻角球鼻闭合

原文附录 B 印刷页42 / PDF47 开头说明：一般未知激波鼻部曲率会要求求整个亚声速/跨声速区；但若给出符合式42的鼻部压力分布，曲率可由该分布计算。式B1 选用

`p_b = p_inf + (p2-p_inf)sin^2(Gamma_b)`。

`Gamma_b` 是体表切线与来流的角度，鼻端为90度。零攻角球鼻，令 S 从鼻尖沿体表向下游递增，Rn 为球鼻半径，phi=S/Rn，则

`h(S)=r_body(S)=Rn sin(phi)`，

`p_e(S)=p_inf+(p_stag-p_inf)cos^2(phi)`，

`p_b(S)=p_inf+(p2-p_inf)cos^2(phi)`。

与光滑直锥相切时 `phi_join=pi/2-alpha`，圆半径与切线均连续，物理压力在接缝两侧同为 `p_inf+(p_stag-p_inf)sin^2(alpha)`，辅助压力相同趋势也连续；其导数在曲率跃变处可不连续。此模型输入是公开工程压力假设，并不等于求解了无粘压力 PDE。

设完整正/斜激波前向热物性模型的鼻端展开为

`p_sh(Gamma)=p2-D cos^2(Gamma)+O(cos^4 Gamma)`。

定比热理想气体由附录 A 式A9（印刷40 / PDF45）给出 `D/(rho_inf U_inf^2)=2/(gamma+1)`；D 在源文献记为 Delta p。避免与源式59 的另一无量纲 D 混用。

写 `cos(Gamma)=k S+O(S^3)`。式66 与 `h=S+O(S^3)` 给出

`Psi_sh = rho_inf U_inf S^2/2 + O(S^4)`，

`sin(Gamma)=1-k^2 S^2/2+O(S^4)`。

把上述展开代入式65，并比较 S 的一阶项，得到

`k^2 = 2(p2-p_inf) / [Rn^2(rho_inf U_inf^2+2D)]`。

等价地，令 `x=S/Rn`、`a0=(p2-p_inf)/(rho_inf U_inf^2)`、`d0=D/(rho_inf U_inf^2)`、`A=a0/(1+2d0)`，则

`sin(Gamma)=1-A x^2+O(x^4)`，

`Psi_sh/(rho_inf U_inf Rn^2)=x^2/2-(A+1/6)x^4/4+O(x^6)`。

所以不存在要人为猜测的鼻端 k；几何 Rn、工况和上述工程压力假设已足够闭合。数值启动点必须由此局部展开确定，并核对启动点缩小/积分容差收敛；不能只把 Psi=epsilon 而让 Gamma保持90度，那会选择平凡/错误启动。

也可以用 `t=(S/Rn)^2` 积分：在 t=0 明确指定上述方程的连续延拓，`d(sinGamma)/dt=-A`、`d[Psi/(rho U Rn^2)]/dt=1/2`。对于离开0的站，仍按原式65–66求导；极近鼻端采用同一展开避免两个大压力之差的舍入。

对于非定比热 NASA7 固定组分模型，数值计算同模型的 `D=-d p_sh/d(cos^2Gamma)|0` 是工程热物性推广。源文献使用理想气体或特定平衡空气相关式；NASA7 固定组分推广不得写成已完整复现源报告的平衡空气模型。导数须有步长收敛与有效物性域检查。

## 附录 B 与真正激波鼻部几何的关系

印刷页45 / PDF50 式B16（已视觉核对）为

`Rn^2/R_sh^2 = [(p2-p_inf)/(rho_inf U_inf^2)](1-Delta/R_sh)^2 / [D/(rho_inf U_inf^2)+1/(B+1)]`。

零攻角球鼻 B=1。因此

`(R_sh-Delta)^2=Rn^2[d0+1/2]/a0`，

`k=1/(R_sh-Delta)`。

这正好验证由式65的体表坐标启动推导。它不是 `k=1/R_sh`。附录 B 印刷44 / PDF49 式B11–B12 的局部法向对应关系是 `r_body ~= (1-Delta/R_sh) r_sh`。所以真实激波的质量坐标与式66体表近似质量坐标，鼻端即具有这个映射差异；不能用近似表反推精准激波半径并叠到三维 UI。

若后续需要真正激波鼻端位置，可增加原式57（印刷22 / PDF27）与附录 B 的闭合。B=1 时式57化为

`c=Delta/R_sh = sqrt[(rho_inf U_inf^2/(2D))(rho_inf/rho2)] / [sqrt(1+E)+sqrt(1+rho_inf U_inf^2/(2D))]`，

`E=(Delta h/U_inf^2)(rho2/rho_inf)(rho_inf U_inf^2/D)-1`，

理想气体 `E=(M_inf^2-1)^2/[2M_inf^2(1+(gamma-1)M_inf^2/2)]`。

该表达式是从原式57的B=1线性积分直接化简，等价于式61（印刷23 / PDF28）；Delta h 是近正激波静焓导数系数，不是总焓差。由已知 L=1/k=R_sh-Delta，得 `R_sh=L/(1-c)`、`Delta=c R_sh`。下游真正三维激波还需对应体/激波弧长度量，不能将式65的表直接称 shockmesh。

此真实几何扩展可保持为后续步骤，不是当前 body-associated Gamma/Psi 内核的前置依赖。

## 应保留的物理和验证边界

- 限定零攻角、轴对称、共轴球鼻与直锥；只计算迎风、附着、有效正/斜激波域。圆柱后体、有攻角、非圆截面、肩部、尾迹与分离区没有因此自动闭合。
- gamma_inf（连续的激波角函数）不能通过每面元切楔/切锥重新求根或无来源平滑来拼接。需要用同一前向法向跃迁与切向速度连续关系得到 p_sh。
- 一个可用于独立校核的理想气体渐近角是 `sin^2(Gamma_inf)=sin^2(alpha)+cos^2(alpha)/M_inf^2`，源于直锥上的常数 p_b 与 `p_sh=p_b` 的式65平衡。这是本工程模型的渐近极限，不是 Taylor–Maccoll 精确锥流。
- 质量表需要单调 Psi，角度应不增并保持法向来流Mach>1；边界层质量要求超出表的既有范围时返回真实未解状态。切勿把源式66精确激波积分与体表近似在不同步骤混用。
- CR-112180 印刷37 / PDF42 结论原文：“This distribution was found to yield accurate shock shapes for a sphere, but it has not been checked for other body shapes.” 因此有球体激波证据不能冒称球锥波形实验验证。报告钝圆锥/椭圆锥的热流对照是复合方法证据，仍需本项目工况与测点重做对照。
- 原源方程为近似激波几何与熵层吞入算法；当前只实现质量表不足以证明全表面压力、温度、热流或界面梯度已联合匹配。

原页图：`appendixB-47.png`、`appendixB-49.png`、`appendixB-50.png`、`appendixB-51.png`、`eq57_61-27.png`、`eq57_61-28.png`、`appendixA-44.png`、`appendixA-45.png`、`conclusions-42.png`。原报告中的beta表示周向/流管坐标，项目角度命名应继续使用Gamma避免混淆攻角。

## 2026-10-13 驻点代数耦合审查：质量/热力闭合不能代替动量闭合

对 `solveEntropyMatchedStagnation.m` 的只读审查确认，设 `Bp=(ps-pinf)/Rn^2`、`p0=ps+C*Psi+O(Psi^2)`，质量吞入 `Psi=rhos*a*delta_m(a)*s^2` 给出

`a^2=2*Bp/rhos+2*C*a*delta_m(a)`。

此代数推导正确；驻点局部相似剖面无量纲形状独立于 `a`，因此质量厚度正比 `a^(-1/2)`、热流正比 `sqrt(a)`。但 TP-3271 印刷3–6 / PDF5–8 的源局部相似方程采用 `A_pressure=-beta`，要求外缘沿程动量相容。当前物理压力近鼻端为 `pe=ps-Bp*s^2+O(s^4)`，耦合速度为 `Ue=a*s+O(s^3)`，故

`pe_s+rhos*Ue*Ue_s=(rhos*a^2-2*Bp)*s+O(s^3)`，

`relativeMomentumMismatch=1-a_fixed^2/a^2=2*C*delta_m/a`，

`A_pressure+beta=C*delta_m/a`，其中 `a_fixed=sqrt(2*Bp/rhos)`、`beta=.5`。

这项失配在 `s→0` 的相对极限仍有限。驻点点值 `Ue=pe_s=0` 无法证明相容，也不能靠减小数值启动点消除。源熵层横向涡量/相互作用尚未求解，因此新耦合剖面仅为带条件的研究诊断；默认 `maximumMomentumMismatch=1e-6` 不放宽，主场保持未解 NaN，诊断另存，不能计入全面积解。典型 `Tinf=220 K,pinf=5000 Pa,M=4,Rn=.1 m,Twall=300 K` 工况的定比热失配约 `0.0199507`，NASA7 约 `0.0207709`；通过代数残差检验不改变这一物理结论。

## NASA7 总压鼻端导数的独立守恒参考

以 `w=cos^2(Gamma)`、`N=Uinf*sqrt(1-w)`、`N_w=-Uinf/2` 表示法向上游速度。独立用温度未知量求法向质量、动量、能量守恒，得到鼻端静温 `T2`、法向速度 `v2`、静压 `p2`。在 `w=0` 隐式微分可解二阶线性系统：

```text
[ cp2             v2                   ] [ T2_w ]   [ N*N_w                           ]
[ N*R/v2          N*(1-R*T2/v2^2)       ] [ v2_w ] = [ -(R*T2/v2+v2-2*N)*N_w           ]
```

完整总焓对应的总温 `T0` 与 `w` 无关。由 `p0=p2*exp((s(T0)-s(T2))/R)`，直接得到

`B=dp0/dw=p0*(N_w/N+(1-cp2/R)*T2_w/T2-v2_w/v2)`。

该参考不复用被测 `recoverShockEntropyEdgeState` 或 Richardson 导数。2026-10-13 用 Python/SciPy、独立 NASA7 多项式与温度守恒根计算上述典型工况：

```text
Uinf = 1189.3621820118547 m/s
T0   = 894.7849544902645 K
T2   = 866.6208249435045 K
v2   = 250.874925823623 m/s
p2   = 93375.58002330553 Pa
p0   = 105754.25191288511 Pa
cp2  = 1114.1957574995502 J/(kg K)
T2_w = -623.7032707758015 K
v2_w = -49.282174065100264 m/s
B    = 187213.83038021036 Pa
```

独立中心差分 `dw=1e-5` 得 `187213.83036536278 Pa`，与隐式导数相对差约 `7.93e-11`。当前保存的被测模块导数 `187213.8290835816 Pa` 与隐式参考相对差约 `6.93e-9`。这是固定组分 NASA7 守恒/关联式数值复算，不能视为平衡化学空气或实验验证。
