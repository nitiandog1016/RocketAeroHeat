# v174 定比热NACA独立入口合同审查

2026-11-06。仅负责新入口的独立结构、实验锚点及代数核对，不修改旧Lemmon/OTH入口，不启动CFD或构建。当前v173 adaptive+节点后处理限步支线保持关闭，不能把新模型作为其冷重启或根因修复。

## 最小适配边界

新caseKind采用`naca-tn3344-constant-cp-spatial-wall-v1`，新物理合同为`ideal_gas_contract.json`/`manifest.idealGas`。初始化独立具名，`priorSuccessfulUpdates=0`；local0对应本系列1次成功更新。`RESTART_SOL=YES`仅用于让引擎读取新生成9列初场，并不表示继承历史数值轨迹。不能携带旧来源步1500、旧系列1501或`samePhysicsAndNumerics=true`。

旧`prepareNACASpatialCase.py`与`su2NACASpatialInterface.py`均硬绑定旧case、schema6、Lemmon核、OTH参考、旧`source_data().freeStream`、`IdealAir()`初始化及v171计划。这些旧断言保持原样。新入口应实现独立来流/初态、配置、模型合同及审核，只复用已过的几何、表读写、原节点身份、有限段功率和有明确模型边界的保存机制。

父任务将新建`prepareNACAConstantCpCase.py`、`su2NACAConstantCpInterface.py`及`nacaConstantCpModel.py`；本审查不改其共享文件。原生包优先准确v171空间Tw高精度二进制，不启用v173节点限步包。所选物理模型必须与包的既存用途元数据分别记录：含项目适配的二进制不能称未修改官方发行包，不能篡改原`runtime.json`来让旧Lemmon模型声明变成定比热。新合同应绑定准确exe、父运行包及原生IDEAL_GAS分派证据。

## 必须同步闭合的运行接线

| 位置 | 现有陷阱 | 新入口要求 |
| --- | --- | --- |
| `su2ThermalRuntime.resolve` | spatial-wall flavor会送入旧IdealAir验证器，旧验证拒IDEAL_GAS | 先按严格新caseKind分派，不放宽旧模型断言；返回实际冻结Tw表 |
| `runEmbeddedSU2.validate_resource_budget` | 旧具名资源上限按caseKind分派 | 新case的预算只能来自新事前计划，不能借旧v171授权 |
| `runEmbeddedSU2`的物性回执/容量 | 目前只读取`ideal_air_contract.json`，容量取`state.thermodynamicModel.checkpointContract` | 新`ideal_gas_contract.json`须实际填入同等回执/容量路径；不能只增加validate函数 |
| `PHASE_FLAVORS` | 容量监督与自动清点以runtime flavor启用，曾发生遗漏 | 复用既有v171 flavor仍要证明分派正确；若新descriptor则明确登记，禁止离线尺寸冒充事中cap |
| `su2IdealAirCheckpoints.inventory` | 未识别新case时落到SST通用列与`manifest.idealAir` | 新case必须早分派至新层流清点，不借旧SST列或旧offset |
| `nacaSavedFrameStudy` | 调用上面的清点，异常需显式开关 | 正常/异常完整帧身份保持；末有限历史不等于终止全场 |
| `su2NACALaminarInterface.inlet_state_check` | 内部明确使用IdealAir和`totalEnthalpy_OTH_Jkg` | 新入口另实现定比热Mach向量/速度/总焓校验 |
| 旧`load_state` | 字段结构可复用，但60–2000K筛选不是新低焓资格 | 新模型读取与字段身份明确；[150,360]K另列开发域，不clip、不冒称物理解 |
| 旧测量/外缘分析 | `column_edge`写死Lemmon/OTH标签，`transport`重用旧物性核 | 新审核不得调用旧measurement_analysis；先做20段真实功率，外缘/归约缺项明确保留 |

新9列初场继续使用`PointID,x,y,z,Density,Momentum_x,Momentum_y,Momentum_z,Energy`，Energy明确为体积总能量rhoE，不是比内能或总焓。保留原全局编号、mesh02、130站空间Tw和完整原壁面；热流正号入壁，原加热实验可能为负。新包若输出相同具名字段，只说明文件结构可复用，能量物理身份仍必须独立。

## 实验锚点的独立重建

固定gamma=1.4、R=8.31451/0.02896546、M=1.97、Re_D=2.28e6、T0=570°R=316.6666666666667K及D=.1016m，Sutherland与Pr=.72保持。T∞、rho、p及U是由模型推导的量，不是必须逐字保留的独立实验测量。

```text
cv = R/(gamma-1); cp = gamma*cv
Tinf = T0 / (1+(gamma-1)*M^2/2)
a = sqrt(gamma*R*Tinf); U = M*a
mu = 1.716e-5*(Tinf/273.15)^1.5*(273.15+110.4)/(Tinf+110.4)
rho = Re_D*mu/(U*D); p = rho*R*Tinf
h0 = cp*T0; e = cv*T
```

已执行独立50位Decimal代数检查`tests/verifyNACAIdealGasAnchorsV174.py`，监督`independent_anchor01`为0.156秒exit0，无新CFD/构建。所得：

| 量 | 独立结果 |
| --- | ---: |
| R | 287.04912678756006 J/(kg K) |
| cp / cv | 1004.6719437564603 / 717.6228169689002 J/(kg K) |
| T∞ | 178.28523385392623 K |
| a / U | 267.6700748601525 / 527.3100474745005 m/s |
| μ | 1.2022255553899439e−5 Pa s |
| rho | .511635944608275 kg/m³ |
| p | 26183.79867211039 Pa |
| h0=cpT0 | 318146.1155228791 J/kg |

Mach/ReD/T0/EOS独立回环最大2.22e−16。原Lemmon/OTH总焓443281.1266472082J/kg不直接放入这个模型，差别包含能量零点，不能把两绝对焓之差当作热力学近似误差。模型差应使用同参考温度Δh。原`contract.json`及测量文件没有修改。

## 独立初态与可区分小检查

沿用几何C1壁距混合但明确只是初猜：b=smoothstep(min(d/.001,1))，u_x=−bU∞，p=p∞；

```text
T(b) = (1-b)*Tw + b*T0 - b^2*(T0-Tinf)
rho = p/(R*T)
rhoE = rho*(cv*T + |u|^2/2)
```

本基准Tw>T0>T∞，故dT/db=T0−Tw−2b(T0−T∞)<0，解析初态在[T∞,Tw]内；不需要截温或读取旧rhoE。独立检查对Tw两端、6个b值共12态，T/p/全焓守恒量往返最大4.44e−16。错误的单位Mach向量−1与错误沿用OTH总焓两个示例均有明显辨别力。**这不是生产新validator的拒绝测试，也不是模型差/实验资格。**准确源与报告在`outputs/v174_validation/independent_anchor01`。

待新入口落盘后最少应检查：正确新合同可解析；旧OTH标签/旧rhoE、错误R或gamma、旧来流p/T/U、单位Mach向量−1、非零prior/来源步、错误模型/包、SST或二维轴线、Tw缺失/错表头/覆盖不足、初场列/编号改变、未登记资源与cap、阶段/保存合同不符均在启动前拒绝。复用旧拒绝规则，不为凑数量重复全部历史测试。各拒绝测试须对真实新validator有区分力，不能只验证测试自己写的条件。

## 模型域与物理资格边界

父计划事前给出的发展比较门槛为cp/k/同参考Δh≤2%、声速≤1%、相同明确edge/recovery状态下Nu/√Re模型因子≤3%。它们是开发模型误差额度，不是原文实验标准，也不是原约15%Nu与1%Re probable error的统计合成。

原输入T∞—T0及实测Tw决定主比较范围；[150,360]K是额外开发域，不保证最终全场范围。新审核应分别保存严格能量/EOS/黏度/热流一致性和实际温域是否满足；超域不裁温，不把IDEAL_GAS代数可以继续计算称物理适用。旧非收敛场温区只能作诊断。

给定Tw、壁流体T、未解固体Tw和独立首邻y+继续分列。后6段中心、原输运表、真实局部外缘、网格/外域和公开归约仍缺，不因分段功率生成而称comparisonReady。原−8、热梯度2e−10、功率2e−9和近壁门槛保持；任何低焓层流通过不授予真实M5/20°、分离再附、底面或突起资格。
