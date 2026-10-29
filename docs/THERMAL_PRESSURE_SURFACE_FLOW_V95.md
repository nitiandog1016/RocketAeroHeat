# v95：同物性压力恢复与曲面流线衔接

本版已完成六组MATLAB核验、独立Python复算、旧定比热回归、三维工作台映射和实际浏览器下载核对，随源码交付。它让已有内部NASA7切锥状态能够接入v94表面流线；仍是光滑附着补片的条件流线模块，不能据此增加整箭合格热流面积。

## 模块及使用范围

| 文件 | 用途 |
| --- | --- |
| `createSurfaceFlowGas.m` | 显式声明定比热或固定组分变比热、温域、参考温压及积分约定；域外不静默回退 |
| `surfaceFlowAirThermalProperties.m` | 保留原cp逐值结果，连续匹配高低温段的焓、标准熵积分常数 |
| `surfaceFlowThermodynamicState.m` | 提供同一状态的cp、h、s及偏导，用于Euler方向方程和沿程一致性检查 |
| `surfaceFlowTotalStateFromAnchor.m` | 从实际p/T/rho/U转换一次总状态约定，保存原锚点和总状态差量 |
| `bindIsentropicSurfacePressure.m` | 从连续压力及两向导数、单一p0/T0恢复相容Te、rho、U及解析导数 |
| `marchPressureDrivenSurfaceFlow.m`、`pressureSurfaceFlowRHS.m` | 在原五射线RK4和横向尺度方法中显式支持变cp，原定比热路径保留 |
| `exportThermalPressureSurfaceFlowStudy.m` | 读取已通过的保存结果，接入同侧栏五步骤、三维站点拾取、曲线、探针和CSV |

这些文件位于`regionalFlow/surface/`；中文模块注释说明每段职责。原调用接口、主轴X、鼻尖+X、零攻角来流−X和SI单位保持。

适配器要求一块光滑附着补片共享同一总焓和熵。输入包括曲面的第一/第二导数、连续p及其两向导数、明确入口方向、来源和适用声明。禁止同时传入另一套T/rho/U、逐站总压字段、失败锚点或不匹配的物性约定。只支持声明温域中的非驻点`0 < p < p0`；激波、分离、驻点和光滑性边界不能跨越。

本版不求未知压力，不证明完整外流连续方程，不解激波层、边界层、热流、恢复温度或壁温。两个平面是解析压力验证输入；锥面压力来自实际内部Taylor–Maccoll。零攻角直锥不能替代真实20°箭体的横流、收缩、后体和分离再附。

## 1000K分段问题及处理

最初两次MATLAB任务实际exit1，保留原进程、源码和诊断MAT。直接把原NASA7的h/s接入连续压力恢复，在跨1000K时使横向动量误差由约1e−11跳至`7.7951085564e−7`，超过原`2e−7`核对门槛。单位流向漂移约`1.74e−12`，不足以解释这一跳变。

原多项式在1000K的高段减低段差约为：cp `−2.02629242e−5 J/(kg K)`，h `−0.0234989397 J/kg`，标准熵 `−0.00113042636 J/(kg K)`。有限位数积分常数不完全衔接，使压力—熵恢复与局部cp导数不严格配套。

新流线显式采用`continuous-project-cp-v1`：

- cp逐值保持生产库的结果；没有修改`airThermalProperties.m`。
- 1000K以上h加`0.02349893969949335 J/kg`，标准熵加`0.0011304263616693788 J/(kg K)`，使积分连续。
- 80–150K保留原振动冻结近似及连续低温积分，150–200K仍有原外推；固定组分，不含解离、电离或化学反应。默认变比热声明范围80–2500K。
- 未放宽动量、状态方程、沿程总焓/熵或流线有效前缀门槛；旧生产总状态必须显式转换后才能使用。

固定总状态下，恢复关系为`S°(T) − R ln(p/p_ref) = S0`及`h(T) + U²/2 = H0`。导数使用同一cp：`dT = R T dp/(cp p)`、`dU = −dp/(rho U)`。不逐站重置p0，也不裁温度、动能或压力以通过检查。

## 实际内部切锥锚点

工况：M∞=5，T∞=220K，p∞=1200Pa，零攻角15°直锥。真正调用`solveThermallyPerfectConicalShock`，没有制造锥面压力。

| 锥面物理量 | 实际输入锚点 |
| --- | ---: |
| p | 4439.071465653153 Pa |
| T | 334.8996703819758 K |
| rho | 0.04617637810928356 kg/m³ |
| U | 1406.958421345186 m/s |

旧p0/T0为`593390.6492226454 Pa / 1240.453897745665 K`；连续积分约定转换后为`593392.9468888028 Pa / 1240.453877789539 K`。差量分别约`2.29766615736 Pa`和`−1.99561259e−5 K`。原物理锚点逐值保留，差量属于热力总状态约定转换，不是新增外流解；恢复时的微小求根舍入单列。

## 验证结果和证据范围

| 检查 | 结果与边界 |
| --- | --- |
| MATLAB六组 | 定比热平面/旋涡、NASA7跨1000K/150K平面、真实内部锥流两个标签通过；另含六项拒绝和总压重置失败前缀 |
| 独立热力与轨迹 | Python独立多项式、连续积分、Brent恢复和横向动量ODE复算通过；两平面动量相对误差约`1.55e−11`、`2.71e−11` |
| 独立切锥 | 从法向激波守恒根与SciPy Taylor–Maccoll重新复算，最大状态相对差`8.63e−13`；是同模型数值复核，不是实验验证 |
| 定比热回归 | 原16组/80射线通过；32672个保存对象逐值相同，仅排除MAT头和16个新增物性模式字段，没有放宽容差 |
| 历史v94输入 | 用明确Git来源核对旧冻结输入，标明不是当前MATLAB重算；当前JS也直接通过旧数据复核 |
| 三维映射 | 3模型、4流管、124站点、120线段逐值对应保存MAT；完整锥柱闭合、定向、正体积及Euler特征通过 |
| Node工作台 | 15次流程、15次无关/过时消息拒绝、93行CSV通过；未知值不补零，条件反算字段和逐变量方法列通过 |
| 实际浏览器 | 五步骤同页切换、低温末站、三维锥面点选第2管第18站、31行真实CSV下载及数值/方法逐项完全相等通过 |

浏览器探针分列：输入p、条件反算Te/rho/U、未计算Tr、未设置给定Tw、未知热流/求解Tw。压力方法与外缘状态方法各自展示。首选视角使两条锥面流管可见，网格与结果阶段沿用模型相机。当前新标签页未记录警告/错误；不覆盖此前历史浏览器记录或据此宣称所有环境控制台均无错误。

所有实际MATLAB任务均已终态；最终UI各次导出只读取同一`thermal03`结果，没有重新积分NASA7。修复前的标签、相机、方法列版本及原失败保留于验证索引，不能把旧图当作最终验收。

## 从仓库复现

小型核心MAT及原失败MAT位于`tests/data/thermal_surface_flow/`，共约387KB；没有工程ZIP、大矩阵或完整隔离工程副本。独立审核无需MATLAB，依赖Python的NumPy、SciPy：

```text
python tests/auditThermalPressureSurfaceFlow.py tests/data/thermal_surface_flow docs/validation/v95/raw/thermalContinuous95_process.json --failure tests/data/thermal_surface_flow/failure/thermal_case1_before_invariant.mat --report outputs/thermal_frozen_audit_new.json
```

新报告路径必须尚不存在。已交付的冻结复核见`docs/validation/v95/frozen_input_independent_audit.json`。

在MATLAB重新计算核验，或只从冻结MAT重导出工作台：

```matlab
setupSurfacePressureFlow; addpath('tests');
testThermalPressureSurfaceFlow('outputs/thermal_check_new');
exportThermalPressureSurfaceFlowStudy( ...
    'tests/data/thermal_surface_flow/thermal_pressure_surface_flow_evidence.mat', ...
    'outputs/thermal_view_new');
```

两个输出目录都必须尚不存在；只想看保存结果时可省略`testThermalPressureSurfaceFlow`。以`outputs/thermal_view_new`为目录运行本地静态服务，打开`/workbench/index.html`。当前实际验收页面为`http://127.0.0.1:8777/view05/workbench/index.html`。

```text
node tests/testPressureSurfaceFlowWorkbench.js outputs/thermal_view_new/workbench outputs/thermal_view_node_new.json
```

## 下一步及当前物理状态

固定p0/T0补片不能携带不同入口标签的熵差。下一步沿入口标签β显式携带`H0(β)`、`S0(β)`，通过同步射线族及横向尺度恢复正确的物理横向梯度，再接真实附着区域和边界层入口。候选关系为`∇s S = S0′ b/h`、`∇s H0 = H0′ b/h`、`∇s T = T(∇s S + R∇s p/p)/cp`，其中`b=n×t`、h为标签横向尺度；这些关系尚待实现和独立验证，不能把固定坐标下的辅助导数冒称物理横向梯度。

随后继续内部有来源的压力/横流—收缩—后体衔接、三角面元和边界层入口；分离再附及底面采用各自适用模型。不能用流线数量或几何闭合替代全区域物理解。正式合格热流面积仍为66.5626327424%，本轮新增0；Q41仍19次接受/82次主NS且未稳态。单时刻全面积和科研UI仍为当前优先目标，全周期目标保持active。
