# v259：固定Eilmer空间／块推进候选拒绝

2026-11-17，基线fdc26459232763e01346d7524614db9522d9b652。本轮完成v258指定的唯一H1候选核对：**当前固定lmr没有块推进执行实现；同仓库旧Eilmer的实现也不能直接闭合原箭体所需反馈，故不采用、不构建、不试跑。** 该结果不否定Eilmer的所有功能，也不表示所有区域方案无效；完整气动热总目标仍active，新增合格面积0。

## 实际源码与官方文档的区别

固定提交`f53f4609a0331d48efee69a4e4f3c3598378cc03`的20份相关源码/文档/算例来自已有WSL仓库的git对象，而非可能修改过的工作树。逐份Git blob SHA-1与官方GitHub完整树匹配；另有7处官方连接器原文行段比较通过，SHA256另存。[源清单](validation/v259/evidence/source01/report.json)保留准确路径、字节数及对照范围。

最关键的是两套入口不能混用：

- `doc/lmr/config-options-user.adoc:204–239`描述`config.block_marching`、每次两片活动块和主超声速i方向；配置枚举也接受这个名字。
- **`src/lmr/commands/runsim.d:373–375`实际输出`NOT IMPLEMENTED: block marching is not available right now.`，随后`exit(1)`。** 这证明固定lmr该执行分支不可用。没有启动一次必然失败的计算来重复证明它。
- `src/eilmer/main.d:763–765`会调用`march_over_blocks()`；实现位于旧`src/eilmer/simcore.d:801–926`。同一仓库中的该实现、旧e4shared算例及文档，不自动成为已安装lmr或v210隔离相态GasModel主程序的能力。

## 旧Eilmer具体做法及本项目拒绝理由

| 要求 | 固定代码证据 | 对当前原箭体的决定 |
|---|---|---|
| 真实下游反馈 | 首先只保留两片活动块；下游边界调用`pushExtrapolateCopyAction`；前进时将旧上游片`active=false`，没有返回全域反馈循环 | 不能据此闭合v244所列1984界面及827入射面的责任；不以全局M5取消近壁/分离入射信息 |
| 临时边界含义 | `boundary_condition.d:258–272`保存并清空原对流/黏性边界动作，只加入零阶`GhostCellExtrapolateCopy` | 这是人为的临时外推；不是相态特征闭合，也不是实际下游求解反馈 |
| 低内存顺序运行 | 非MPI入口将所有块纳入`localFluidBlocks`；在march前为所有本地块初始化workspace、网格/流动数组并读取解 | 两片“活动”不等于只驻留两片；不能按活动片比例宣称解决8GiB问题 |
| 复杂原网格和周向连接 | 要求单一规则i/j/k块阵列、固定编号关系和east连接；原始母网格不是现成该布局 | 没有原19752面、32层2µm及尾缘资产映射证明，不能换成平板或简化箭体后称原模型完成 |
| 可选剖面播种 | `propagate_inflow_data`将上游ghost状态铺入新块；`UFluidBlock.propagate_inflow_data_west_to_east`明确未实现 | 该选项不能提供原非结构区域的现成真实入口；复制初值也不是已解完整剖面 |

两个活动片之间仍有原生连接交换，**不是每块完全独立**；拒绝针对冻结旧上游和临时下游边界不能证明全局所需反馈闭合，不能夸大成代码完全没有任何反馈。固定代码也没有直接检查“必须SFluidBlock”，不能一概断言所有非结构网格都不能参与；已明确的限制是规则块拓扑和非结构可选传播函数未实现。

`src/eilmer/block.d:20–22`将active定义为“参加时间积分”；`simcore.d:182–186,221–225,245–256,299–308`保留全本地块的数据/初始化。因为实现与物理入口门槛已经拒绝，本轮没有再做精细峰内存预测，没有启动进程制造OOM，也没有把v211复数NK必要载荷套给旧实数时间推进器。

## 官方黏性算例究竟证明什么

同提交`examples/eilmer/3D/flat-plate-laminar/lam_flat_plate_march.lua`确实声明三维、黏性、理想空气、Mach约4、Tw300K、平板、22×4×1块，并开启块推进。入口`T∞=300K、p∞=1013Pa、u∞=1390m/s`；壁面固定温度，侧向slip，无真实箭体驻点/α20°/曲率/收缩/底面，也未在此例指定当前SA-neg/WITHFT2。

`post.sh`提取x=1m的边界层剖面，与Schetz书中CLBL数值解比较。该参考是书中数值对照，不是当前M5/30km实测箭体热验证。`readme.txt`所写15942秒是2020年作者在其机器上的记录，不能当本机预算/峰值或本轮已运行成绩。文件名含march的2D USG算例虽然转换了结构块，但其文件没有开启`config.block_marching`；不能从名字推出它已验证USG块推进。

v210气体实/复数helper通过、v211的必要资源账和v216七域索引均仍有效于各自范围；本轮未重做物性、ABI、网格或流动，不把这些旧通过项转授本候选。v219转向SU2的原理由及其后真实失败保持。

## 收口与剩余工作

本轮拒绝`fixed-lmr-block-marching-provider`，也拒绝直接移用旧`e4shared-two-slice-simple-outflow`来填当前原箭体入口。缺项仍为：能在本机预算内发展真实三维热壁剖面、保留必要上游/下游反馈，并有可复现数值和验证证据的区域方案。没有新p/q/T值可接入UI，因此不复制/重绘旧工作台，不改已过UI。

下一项见[唯一下一决定](validation/v259/NEXT_DECISION.md)：不再换品牌或追查这个已拒入口，转向现有相态量热反算的稳定实现。这不是声称它已被证明是v258方向失败根因；需要同一数据/物理域下直接给出精度证据后，才考虑再次接流动算子。

## 程序、来源和资源

本次仅新增具名只读源收集脚本，20份源共约0.6MB，整个收集输出743436B（后续许可/归档另计），低于32MiB。监督收集4.453秒、Windows采样峰73764864B；每个Linux git读取由taskset绑定CPU0，prlimit限制512MiB地址空间/10秒CPU并用15秒timeout，Windows采样不冒充Linux峰RSS。没有CFD、网格、MATLAB、库构建或物性原生查询。

直接发现阶段另有一次旧E盘wsl.exe路径不存在，命令未启动；随后改用系统现有`C:/Windows/System32/wsl.exe`，两次定向git搜索各约3.734/1.750秒，均返回源码。WSL曾输出ERROR_INVALID_DATA及localhost代理/NAT警告，未阻止离线git读取，不作系统修复、不解读为求解器失败。一次工具JS多余括号及一次PowerShell花括号路径解析失败均未执行脚本；按正确语法读取原文件，无数值重跑。网络读取、开发阅读与归档时间分列，不包含在4.453秒内。

固定上游采用GPL-3.0-or-later，COPYING和完整gpl.txt均保留；后者Git blob也与官方一致。只保存定向源码与小证据，不复制工具链/整仓库。29项、engineering54/候选、原完整壁/层/分母及并发Web/OpenCode保持；仅本地，不推送，不删除独有结果，不新代理/computer-use。
