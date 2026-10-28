# v82真实20°箭体精确通量研究包

本阶段保留MATLAB计算后端和离线三维工作台，解决真实收缩段旧Davis波界HLLC近似星压受限的一个数值问题。新增显式精确理想气体Godunov内部通量，以及可选的primitive切向非线性候选路径。原HLLC和保守候选仍为默认；给定墙温、未知外缘/恢复/求解壁温分别保留。

**本包不是稳态或实验合格的全面积气动热环境。** 常gamma、冻结mu层流研究的两条三次修正均未稳态。保守路径最大残差下降约61.3%，primitive路径RMS略低而最大残差更高，不能宣称整体加速。正式合格热流面积仍66.5626327424%，本阶段新增合格面积0。

## 查看和重建

- 直接打开`example_exact_v82/conservative/index.html`查看保守路径对照；`example_exact_v82/primitive/index.html`查看primitive路径对照。
- 在MATLAB工程根目录调用`demoRocketExactFluxStudy`，重建保守路径研究页；输出目录必须未存在。
- 选择另一端点：`demoRocketExactFluxStudy('outputs/primitive_demo','validation/v82','primitive_krylov01')`。
- 没有本机outputs输入时，入口直接读取包内validation保存的程序内部几何/检查点，不需要外部CFD数据。

两页均比较旧HLLC源状态、相同Q的新方程基线和新方程三次末态；历史从新方程零次开始，不把旧六次HLLC拼入新残差曲线。五步骤侧栏用于查看该冻结研究链，实际新计算仍由MATLAB执行。压力、入壁热流、邻壁气体温度、算子温度斜率均带来源和未收敛身份；Tw300K是输入，Te/Tr/求解Tw保持未知。

## 验收范围

当前工作区的数值/数据/事件测试已完成，详见`docs/validation/v82/progress.json`。本文件为归档前冻结说明；最终压缩包SHA和独立解包实际重算结果，以随包的`v82_release_verification.json`、`v82_isolated_verification.json`为准。只有独立重算真实终态通过后才登记本地交付。

独立验收实际重算精确通量、完整NS集成、保守与primitive各三次修正、合法/拒绝候选和32CV默认源兼容；逐字段核对完整数值对象，仅排除机器耗时/PID。旧HLLC六次末态作为冻结内部初猜，不冒称再次求解。两组原面JSON、逐变量CSV、任意截线和同宿主事件分别核对。

独立科研PNG不是浏览器截图，没有新增浏览器/WebGL视觉验收。数值守恒、公开压力曲线复算和接口测试不代替真实空气/湍流、近壁/网格/外域敏感性、区域热流精度或实验验证。后续继续真实收缩冷区/近壁能量与稳态收敛，逐步补齐主目标的大面积区域；不以本包宣称单时刻全面积或全周期已完成。小版本不推GitHub。
