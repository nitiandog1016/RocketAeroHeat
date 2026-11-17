# 当前交接：唯一恢复入口

更新时间：2026-11-17。最近交付v259：**固定lmr块推进执行分支未实现；旧Eilmer两片推进也未闭合当前区域反馈，具名候选拒绝。** 无新构建/CFD/热面积，完整总目标active。

<!-- records:last_completed=v259 -->
<!-- records:active=none -->
<!-- records:base_commit=fdc26459232763e01346d7524614db9522d9b652 -->

基线为本轮开始前v258；本轮本地提交另看Git。仅本地，不推送。v259只读收集监督终态已知/输入稳定/组空，CLI关闭；没有运行中数值任务。

## 1. 下一项

按[v259 NEXT_DECISION](validation/v259/NEXT_DECISION.md)实现并独立核对**同一相态标准Air的稳定分段量热E/H反算**。准确已链接起点为v258依赖中的`outputs/v220_validation/source/SU2_CFD/include/fluid/PhaseAwareAir.hpp`，SHA256 `7b3301b0e1b499c0398611333911965d113001868da65c85b9b76886de151366`。先核对维护版与真实源；不原位改写历史/共享硬链接。

Cp分段线性，h/e各段为严格单调二次式；现有E/H/S统一TOMS748/48位括区。新E/H用避免小根相减的稳定二次反算，S仍保持原有界方法。物性表、60—2000K及相态联合域、能量参考、输运和拒绝语义不变。复用v258实际微扰及少量节点/端点，用高精度独立参考直接核对精度与实际收益；准确v260预算先登记，当前未准入新helper或CFD。

此处尚未证明量热反算是v258失败根因，不把稳定公式或helper通过直接授予区域资格。有可证同模型精度收益才登记全链集成；不重跑未改的旧组合，不扫参数/物性域，不再追查已拒块推进或轮换品牌。总体路线仍低成本全局外流＋工程边界层＋局部黏性CFD，入口/反馈及全表面缺口仍需解决。

## 2. v259已核对的拒绝依据

- [专题报告](EILMER_MARCH_REVIEW_V259.md)、[计划](validation/v259/PLAN.md)、[最终证据](validation/v259/evidence/final_report.json)。固定GDTk提交f53f4609a0331d48efee69a4e4f3c3598378cc03的20份源/例/文档从已有git对象读取；与官方Git blob逐字匹配，7处官方行段交叉对照通过，GPL-3.0-or-later原文保留。
- lmr文档有block_marching选项，但`src/lmr/commands/runsim.d:373–375`为NOT IMPLEMENTED并exit1。没有启动必然失败任务；同仓库旧e4shared实现不能转授当前lmr/v210相态helper主程序。
- 旧`march_over_blocks`只保留两片活动，推进时停用旧上游片、临时下游改零阶外推；两个活动片间有连接，不是完全独立，但没有已证的全局反向反馈闭合。不能绕过原1984界面/827入射责任。
- 单进程在march之前为全部本地块初始化workspace、网格/状态及读取解；active仅表示参与积分，不能按两片估驻留。规则块拓扑/编号和原面映射未证；USG可选west-to-east播种函数未实现，但不一概宣称所有USG块推进都被代码禁止。
- 官方三维黏性Mach4平板、理想空气、Tw300K例确实启用块推进；CLBL是书中数值对照，不是α20°箭体实验/当前SA验证。作者旧15942秒不是本机实测。文件名带march的2D USG例没有开启该配置。

候选在实现/物理边界门槛已被拒，不做无意义峰值细算；未测新内存峰、未发生或声称新OOM。收集4.453秒，Windows采样峰73764864B；Linux git CPU0/每命令512MiB地址空间/10秒CPU/15秒timeout，Windows采样不冒充Linux RSS。额外定向搜索和网络/开发阅读另计；无求解、网格、MATLAB、物性查询或共享UI修改。

旧E盘wsl.exe路径缺失的直接发现命令未启动，后用系统现有wsl.exe读源；WSL的ERROR_INVALID_DATA及代理/NAT警告未阻止离线git读取，未改系统配置。工具语法错误未执行任何计算，原说明保留。

## 3. 最近数值状态与科研出口

v258的4次完整F初态/sigma/half/恢复已终态：原335726点/705490单元、9878壁节点、1974844自由量；接受更新0。固定方向sigma/half缩放相对差0.90646956，未达1e−3诊断阈值；half不是精确导数，A与half差0.99948081不能当A对精确Jacobian误差。一次纯空间ILU近似矩阵残差比26.5474573；约束行全0、方向坐标往返≤8.89e−16，未发现具体接口错误，唯一根因未知。

原六式残差未改善，3904重构警告保留；原−8、内部线性目标、守恒/热稳定/网格/外域/近壁及实验验证未过。v257固定KINSOL/SPGMR20/无重启/块ILU组合及v258方向诊断关闭，不能追加步长、Krylov、重启或盲续旧25/108/200。

科研页`outputs/v258_validation/workbench01/index.html`：3个原面方向诊断、59256值回读通过，明确非精确导数/非壁热流；原p/q/T仍引用v257未收敛输入。给定Tw300 K、壁节点流体T263.5938—361.0898 K、未知固体Tw分开。v259无新数值，不复制/重绘旧页。v255球帽1306面工程候选保留；v256已拒shock-shape质量表直接拼完整剖面，不能靠鼻帽细化替代真实供体。

## 4. 完整责任和保护

M5/30 km/α20°/β0°/Tw给定300 K，X主轴鼻尖+X；原19752面/1.9646292565715553m²、32层2µm/历史19层尾缘及完整下游保持。当前新合格面积0、profileReadyForCFD=false，旧工程q有值66.5626%不是合格覆盖；未知不补零，候选不累计资格。

v252旧三链父监督终态仍未知；v253拼接拒绝、v244反馈、v242/v243外缘未合格及各历史失败保持，不改写成已证唯一原因。v210/v211/v216的接口/资源/七域通过项仅对原范围有效，不转授新运行。

29项`outputs/v176_validation/protected_initial.json`、engineering54/候选、共享`web/section.js`、`web/viewer.js`、`tests/verifyWebSectionExport.py`和OpenCode保持。小证据在`docs/validation/v259/evidence`；旧完整场/矩阵/执行文件和成熟库仍E盘，Git不表示大场已备份。不推送/整包/独有删除/新代理/computer-use。

按[记录维护](RECORDS_GUIDE.md)和[索引](EXECUTION_INDEX.md)恢复；长期目标未改，不向AGENTS/PROJECT_GOALS追加版本流水账。
