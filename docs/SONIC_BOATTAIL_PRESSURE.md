# 短尾锥声速流线压力候选

状态：**独立模块的关联式数学复算，以及对原图15个实验符号的条件压力比较已测试；未接默认生产链，未完成完整NS/RANS实验或热流验证。** 来源/几何草稿见[原报告核对记录](reviews/NASA_TND2761_SOURCE_CONTRACT_2026-10-24.md)，后续实际比较及工作台见[公开尾锥压力对照](NASA_BOATTAIL_PRESSURE_COMPARISON.md)。阶段工程包与包内实际重生以发布核验JSON为准。本模块不增加生产合格面积，也不提供尾锥热流。

`estimateSonicBoattailPressure.m`复算NASA TN D-2761印刷7—9页讨论、图9/10展示的声速流线膨胀近似。输入是收缩前当地柱面静压与直尾锥半角，采用显式声明的定比热空气γ=1.4；这里的声速流线属于厚湍流边界层内的近似，不是将真实外缘设成Mach1。

先从`nu(Ms)=beta`、`nu(1)=0`求声速流线膨胀后的Mach代理，再计算

```text
p_candidate / p_c = [(1 + (gamma-1)/2) /
                      (1 + (gamma-1)/2 * Ms^2)]^(gamma/(gamma-1))
```

程序直接接受声速端点，不修改旧的`solvePrandtlMeyerExpansion.m`，也不通过把1改成1+eps绕过其超声速输入合同。范围限于原报告直尾锥角0—21°；这一数值范围不是全部角度和所有几何均已物理认证的声明。原报告讨论了18/21°的分离迹象，本模块没有求解分离泡或底区。

## 使用与结果身份

```matlab
r = estimateSonicBoattailPressure(1000, [6 12 18]);
% r.candidatePressureRatio 约为 [0.68899886 0.51324127 0.37976554]
% r.candidateWallPressurePa 为给定柱压乘上述候选压力比。
```

1000Pa仅为该数学演示的给定输入，不是实验测量或本程序上游求解值。实际比较必须使用同工况当地柱面参考，不能拿任意远场压力替代。输出保留原报告URL、算法名称、给定输入身份、适用限制及`productionQualified=false`。声速流线Mach只作为代数代理发布；`heatFluxWm2`、`edgeTemperatureK`、`recoveryTemperatureK`、`wallTemperatureSolvedK`全部为NaN。

## 已执行的校核

2026-10-24在单线程MATLAB引擎50516（首次）和45836（最终）中执行原报告0:3:21°的曲线数学复算，并核对声速零转角、未知量身份及五类非法输入拒绝。第一次独立积分发现15/18°压力相对误差略超2e-12；原源和失败进程保留，将求根精度从1e-12收紧至1e-14后重新实算，没有放宽独立校核门槛。

最终`sonic02`与`sonic_independent02`均实际exit0且源码在运行前后保持。Python使用Prandtl-Meyer微分式的数值积分和Brent求根，而不是复制MATLAB的反正切闭式；最大Mach绝对差约2.89e-15，压力比绝对差约1.94e-15，所有未知热场序列保持null。该结果证明数学复算一致，**不是实验精度**。

复现入口：

```matlab
addpath('tests');
reproduceSonicBoattailPressure('outputs/new_sonic_reproduction');
```

本轮实际证据为`outputs/v79_validation/sonic02/`、`sonic_pressure_independent.json`和对应进程/源码记录。当前独立Python研究脚本读取`outputs/v79_validation/<试验名>`，例如`python outputs/v79_validation/audit_sonic_pressure.py sonic02`；禁止覆盖已有报告，重做完整研究时须保留旧证据或使用新的证据目录/名称。

下一步是完成原实验符号的逐点提取、读图误差和条件合同，再比较这种压力关联及程序内部黏性解。不能据此求出或认证真实20°攻角箭体的压力/热流全区域，也不能把一条有数值的压力曲线算作新增热流覆盖。
