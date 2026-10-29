# v96：非均匀入口H/S与物理横向梯度

本版已完成实际MATLAB计算、独立数值校核、三维工作台及真实下载验证，随源码交付。新增逐入口标签守恒的总焓/熵与物理横向梯度，为真实附着区外流和边界层入口衔接提供接口。本版仍使用明确给定的验证压力，**没有新增整箭外流或合格热流面积**。

## 模块和接口

`regionalFlow/surface/marchLabelledPressureSurfaceFlow.m`接受公共光滑曲面/连续压力采样器、入口曲线和方向，以及`entryTotalState(beta)`。每个入口标签独立冻结一次总状态，使用原v95热物性和曲面Euler内核积分；旧v95数值叶函数没有修改。

每条射线的固定标签辅助延拓可得到正确的沿线状态与转向，但其横向热导数不能作为真实物理场。本函数把它们明确标为辅助诊断，只在`center.samples`及`physicalGradients`中输出相邻标签重建的物理两向导数。不能从`rays`提取辅助导数冒充物理梯度。

`tests/labelledSurfaceBenchmarkProblem.m`提供非均匀剪切/转向及内部斜激波锚点验证；`tests/testLabelledPressureSurfaceFlow.m`执行八组实际核验、均匀极限和拒绝检查。`exportLabelledPressureSurfaceFlowStudy.m`只读取已通过的保存结果，不为看图重积分。

MATLAB文件采用中文模块注释，完整文件在项目中直接更新。X为主轴、鼻尖+X、零攻角来流−X、SI及原网格算法保持。

## 横向梯度的构造

在同一弧长站点，设相邻标签的位移导数为`J=∂r/∂β`，单位流向为t，法向为n，`b=n×t`。真正的横向尺度是`h=J·b`，而不是`|J|`。入口斜置时`J·t`可以非零，标签导数不能直接当作物理横向导数。

沿射线H0、S0守恒，则表面梯度为：

```text
grad(H0) = H0_beta * b/h
grad(S0) = S0_beta * b/h
grad(T)  = T/cp * (grad(S0) + R*grad(p)/p)
grad(U)  = (grad(H0) - cp*grad(T))/U
grad(rho)= rho * (grad(p)/p - grad(T)/T)
```

五条射线位于`β+[-1,-1/2,0,1/2,1]δβ`，内侧与外侧两组中心差分分别构造细/粗标签导数。中心物理梯度重新回代原EOS、沿程动量、总焓和熵合同；没有重置总压、裁温、重归一积分状态或放宽原判据。聚焦、标签不可分辨或相邻射线停止后，流管失效前缀之外保持NaN。

## 入口分辨率与微小残差

实际斜激波锚点的总焓理论上几乎相同，数值差分仍可能留下微小残差。本版不把它自动解释为已分辨的物理变化，也不通过补零隐去原数值。

新选项`maximumRelativeEntryDerivativeDifference`默认0.02。两组标签导数的差值指标为`|D_coarse-D_fine|/3`，在积分前分别核对H/S。允许的源状态一致性尺度为：

```text
kappa_H = 4*relativeStateTolerance*max(1, max(abs(Hj)))/delta_beta
kappa_S = 4*relativeStateTolerance*max(R, max(abs(Sj)), max(cpj))/delta_beta
indicator <= relativeLimit*abs(D_fine) + kappa
```

这里使用原状态一致性阈值，**不是物性精度、统计置信区间或严格误差界**。导数低于相应尺度时标为变化未分辨；探针和CSV同时保留原差分值、指标、尺度及分辨标记。刻意输入的未分辨熵跳跃会以`SurfaceFlow:EntryDerivativeResolution`拒绝，调用者应加密标签或拆分接口。

## 验证结果

| 检查 | 已取得证据 |
| --- | --- |
| MATLAB八组 | 定比热/NASA7精确剪切；三档标签间距的非均匀入口转向；斜坐标/刚体旋转；两个圆柱入口标签组 |
| 真实内部源 | 实际调用NASA7二维附体弱斜激波；M=5、T∞=220K、p∞=1200Pa，转角为`20+8 sin(beta)`度；十个锚点独立求解 |
| 独立复算 | Python独立Brent热恢复、DOP853轨迹、密度比能量根与弱激波角根；最大激波状态相对差约`4.56e-11` |
| 横向导数 | 温度/速度/密度链式差分误差在标签半距后约缩小4倍，六个比值约`4.0000` |
| Crocco校核 | 用保存的邻线速度和物理位移独立恢复法向涡量，最大相对差约`2.04e-6`，符合有限标签间距校核范围 |
| 分辨率门槛 | 新增拒绝入口跳跃；修订前后18378个公共保存对象逐值一致，仅排除两项新增元数据；均匀入口退化到v95逐值相同 |
| 工作台映射 | 最终三模型/四流管/84站点/80线段与保存MAT逐值对应；完整锥柱拓扑、定向、正体积及默认可见性通过 |
| Node | 15步骤、15过时/无关消息拒绝、63行物理梯度CSV、9次新增曲线选择通过；未知输入不变零 |
| 旧数据兼容 | 当前程序直接核对v95的93行和v94的164行CSV、原五步骤与探针，通过 |
| 实际浏览器 | 最终view03同页五步骤、三种新增曲线、圆柱第2管第11站、21行真实下载通过；48列全部对应保存数据 |

上述为解析校核和同模型数值复核，**不是实验验证**。圆柱案例是共压轴向熵剪切层；独立二维激波锚点不能拼成同一三维激波面，也不是已经解决真实20°箭体入口。完整锥柱仅为底图，锥段/底面没有新增结果。

第一轮数学核验已通过，随后补充入口分辨率门槛；最终`labelledResolution96`为terminal/exit0、源码稳定。三个界面导出任务都只读取保存MAT，均terminal/exit0、源码稳定，没有重积分。较早界面检查中的来源词组、Te/Tr/Tw完整说明及符号零比较问题已修正，原文与阶段说明保留。

CSV使用标准数字字符串时`-0`会写为`0`，Node对这些零值明确计数；所有非零数值没有容差放宽。实际最终下载的32个数值列和16个文本/布尔列逐项相等。入口导数差值及尺度的单位分别为J/(kg·标签单位)、J/(kg·K·标签单位)，由`entry_label_unit`说明标签单位。

最终浏览器保留一条`MutationObserver.observe`参数非Node的诊断，与v94控制页同签名；记录没有项目URL/完整来源栈，因此不对本次记录作完全归因，也不宣称控制台全无错误。实际点选、流程、曲线和下载验证分列。详见[验证索引](validation/v96/README.md)。

## 从仓库复现

最终283395字节MAT和分辨率门槛前的280111字节MAT位于`tests/data/labelled_surface_flow/`，用于独立重放及精确比较；没有大型矩阵、工程ZIP或隔离工程副本。Python需NumPy、SciPy，报告路径必须尚不存在：

```text
python tests/auditLabelledPressureSurfaceFlow.py tests/data/labelled_surface_flow docs/validation/v96/raw/labelledResolution96_process.json --report outputs/labelled_frozen_audit_new.json
python tests/auditLabelledResolutionEquivalence.py tests/data/labelled_surface_flow/before_resolution/labelled_pressure_surface_flow_evidence.mat tests/data/labelled_surface_flow/labelled_pressure_surface_flow_evidence.mat --report outputs/labelled_equivalence_new.json
python tests/auditLabelledSurfaceFlowCSV.py docs/validation/v96/raw/final_surface_flow_payload.json docs/validation/v96/labelled_cylinder_tube2_final.csv --report outputs/labelled_download_audit_new.json
```

MATLAB重算测试，或直接从冻结输入重导出界面：

```matlab
setupSurfacePressureFlow; addpath('tests');
testLabelledPressureSurfaceFlow('outputs/labelled_check_new');
exportLabelledPressureSurfaceFlowStudy( ...
    'tests/data/labelled_surface_flow/labelled_pressure_surface_flow_evidence.mat', ...
    'outputs/labelled_view_new');
```

只查看保存结果时可省略第一项测试；输出目录必须尚不存在。以`outputs/labelled_view_new`为静态服务目录，打开`/workbench/index.html`。当前实际核验页面为`http://127.0.0.1:8778/view03/workbench/index.html`。可在生成后运行：

```text
node tests/testPressureSurfaceFlowWorkbench.js outputs/labelled_view_new/workbench outputs/labelled_node_new.json
```

## 下一步及全面积目标

数学载体已支持入口标签间熵/总焓差，接下来应使用真实箭体、实际内部入口来源和连续压力补片，核对入口与曲面接口、横向动量及边界层所需条件，再映射到真实三角面元。不得继续用新增制造场、锥流案例或单纯覆盖审计替代这一主线；零攻角或独立激波锚点不能直接当成真实20°全箭外缘。

本函数不解未知压力、完整外流连续方程、激波层、热流、壁温、分离再附或底面。每个真正完成的区域仍须立即接回工作台，按变量保留算法与资格。正式合格热流面积仍66.5626327424%，本轮新增0；Q41仍19接受/82主NS且未稳态。当前单时刻全面积与科研UI优先，总目标保持active。
