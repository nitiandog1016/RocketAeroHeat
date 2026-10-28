# 鼻部—后体熵匹配：v52 独立模块及验证边界

2026-10-13。当前为**已测试的独立研究阶段**。默认网格、外流、生产热流和科研工作台沿用 v51；本版没有增加默认算例的已求解面积。全面积计算、鼻部到后体完整闭合以及实验比较仍未完成。

## 模块实际计算的内容

| 完整 MATLAB 文件 | 计算内容 | 当前状态与边界 |
|---|---|---|
| `solveAxisymmetricEntropyShockTable.m` | 解析相切球锥几何、内部压力模型、近似激波穿越角 Γ 与每弧度质量坐标 Ψ | 零攻角、物理 Mach≥3、15°–45°球锥；按 CR-112180 式64–67 使用物面关联近似，不输出真实三维激波面 |
| `recoverEntropyLayerEdgeState.m` | Ψ=ρe Ue h δm 的表内根与同熵外缘 T、ρ、U、总压 | 输入压力、度量、质量厚度须逐量溯源；数值唯一性不是数学证明；多根、表外、域外、错工况拒绝 |
| `solveAxisymmetricLaminarMassProfile.m` | 同一速度/总焓剖面的物理 δ、δ*、θ、δm、热流与剪应力 | TP-3271 的 Levy–Lees 局部相似层流 BVP；沿程非相似项未求解，0≤β≤1、附着单调支；严格检查梯度合同 |
| `solveEntropyMatchedStagnation.m` | 鼻端质量吞入导致的应变率代数诊断与条件剖面估计 | **不能视为已求解驻点模型**：外缘动量极限失配时主热流/厚度保持 NaN，条件估计另存；始终不计入表面覆盖 |
| `standaloneAirThermalScope.m` | 独立内核真实物性范围说明 | 与生产侧壁回退策略分开，避免来源描述暗示独立内核也会回退或填值 |

以上均不接入 `solveExternalFlow.m` 或默认 GUI。坐标为主轴 X、鼻尖 +X、零攻角速度 −X；单位 SI。壁温是给定边界，不是本模块求得的壁温。

## 激波质量表与压力的区别

来源：[NASA CR-112180 原报告](https://ntrs.nasa.gov/api/citations/19730005231/downloads/19730005231.pdf)，印刷24–28、39–40、42–46页，式64–69、A9、B16。原报告明确式66采用物面关联近似，其球体检查不能自动证明球锥波形精度。

物理压力采用同热物性正激波驻点总压的修正牛顿模型 `pe=pinf+(ps-pinf)*mu²`。式65所需辅助压力 `pb=pinf+(p2-pinf)*mu²` 用正激波静压 p2，单独存储，不标为物面或外缘物理压力。

鼻端展开使用 `t=S/Rn`、`q=t²`、`psi=Psi/(rhoinf*Uinf*Rn²)`：

```text
u=sinGamma=1-A*q+O(q²)
psi=q/2-(A+1/6)*q²/4+O(q³)
A=[(p2-pinf)/Q]/[1+2D/Q], Q=rhoinf*Uinf²
D=2Q/(gamma+1)  （定比热）
```

NASA7 的 D 来自同热物性前向激波压力导数，保留 Richardson 收敛误差。这是固定组分热物性推广，不是原报告平衡空气复现。实际鼻端 Ψ=0、Γ=90°保留，不添加质量常数。

## 物理质量厚度的同剖面定义

来源：[NASA TP-3271](https://ntrs.nasa.gov/api/citations/19930006245/downloads/19930006245.pdf)，印刷3–6页，式6–7、17、25、54、56–61。当前将沿 ξ 非相似项置零，不能称原报告的完整非相似求解器。

```text
xi=integral rhoe*mue*Ue*h² ds
K=sqrt(2xi)/(rhoe*Ue*h), F=u/Ue, D=rhoe/rho
delta=K*integral D deta
deltaStar=K*integral (D-F) deta
thetaPhysical=K*integral F*(1-F) deta
massThickness=K*integral F deta=delta-deltaStar
```

所有积分采用明确的有限外缘，不用原 Thwaites 参考温度厚度，也不乘固定5.27。外缘同时检查速度、总焓、静温、密度、剪切及热梯度；ηmax只是数值边界。近零热梯度采用显式绝对/相对数值容差，保留极小热流是否已被数值精度解析的状态，不裁剪物理热流。

驻点 `h~s,U~a*s,xi~rho*mu*a*s⁴/4` 给 `K=sqrt(nu/(2a)), beta=.5`。原始 NASA7 焓零点可为负；数值 g 的一致移位不会改变温度、热流或物理总焓。

## 已发现的鼻端耦合缺口

内部热力与质量展开产生代数关系：

```text
C=4*A*(dp0/dcos²Gamma)/(rhoinf*Uinf*Rn²)
a_fixed²=2*(ps-pinf)/(rhos*Rn²)
a²=a_fixed²+2*C*a*delta_m(a)
```

同一局部相似条件剖面的 δm 正比 a^(-1/2)。该代数关系可以求根，却不保证源边界层方程在同一压力曲率下闭合。真实压力参数极限是 `Apressure=-(ps-pinf)/(rhos*a²*Rn²)`；耦合后与假定 `-beta=-.5` 有有限失配，不能用 s=0 的两个零梯度隐藏。

M4、220 K、5000 Pa、Rn=0.1 m、壁温300 K的定比热诊断得到应变率7180.36011 s⁻¹，固定熵值7108.37273 s⁻¹。条件热流约202192.347 W/m²、条件质量厚度0.000323461947 m；外缘动量极限相对失配约**1.9951%**。NASA7也存在同类失配。默认严格梯度门限下，主场保持未解，以上值仅在 `conditionalEstimate` 中保存；`eligibleForSurfaceCoverage=false`。

下一步需要有来源的熵层横向速度/涡量及非相似耦合，或另选具备独立适用性的区域工程闭合。不能放宽动量门限、平滑速度接缝或以条件值补全默认云图。

## 实际验证证据

对应 MATLAB 测试及保存日志位于 `docs/validation/v52/`，完整数值数据附本地工程包：

- 内部激波表：鼻端展开、定比热直锥解析不变量、几何尺度协变、独立物理质量积分、启动长度/压力表/容差收敛及拒绝。质量积分相对误差9.21×10⁻⁷；直锥不变量误差1.51×10⁻⁷。有限长锥不能当成无穷远；该工程极限不是精确切锥解。
- 熵匹配核心与真实生成表接口：定比热/NASA7、守恒、数组/来源合同、缩放、多根和域外拒绝。接口试验制造了质量厚度，不能证明物理边界层耦合。
- 独立局部相似剖面：Blasius/Falkner–Skan壁面斜率、Pr=1 Crocco关系、NASA7物理焓、物理坐标积分、驻点极限、静温/密度外缘与ηmax收敛。它只证明限定的局部相似方法。
- 驻点代数诊断：独立完全气体总压导数、NASA7温度未知量激波根与守恒隐式导数、有限小s热力展开、严格动量拒绝、Kn筛选及条件结果不计覆盖。NASA7独立导数相对差6.93×10⁻⁹。

本轮没有重新测试整个旧工程或新增 UI 功能，也没有实验精度通过。公开球锥数据已开始转录，详见 [实验来源审计](EXPERIMENTAL_SPHERE_CONE_DATA.md)：56个印刷热流值、4个缺测位置；名义上游47.34 K低于当前NASA7范围，不能裁温使其通过。

## 运行独立验证

在项目根目录的 MATLAB 中运行：

```matlab
addpath('tests');
testAxisymmetricEntropyShockTable;
testShockEntropyEdgeState;
testEntropyLayerEdgeState;
testGeneratedEntropyMassMatching;
testAxisymmetricLaminarMassProfile;
testEntropyMatchedStagnation;
```

测试通过不代表全面积求解完成。总目标保持 active；底面/尾迹、分离再附、真实鼻部到后体匹配和公开实验对比继续推进。
