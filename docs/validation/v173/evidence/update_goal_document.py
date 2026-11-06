"""本轮交付：只更新目标文档的当前入口，历史事实保留。"""
from pathlib import Path

p = Path('E:/RocketAeroHeat/docs/PROJECT_GOALS.md')
lines = p.read_text(encoding='utf-8').splitlines()

def replace(prefix, text):
    found = [i for i, line in enumerate(lines) if line.startswith(prefix)]
    assert len(found) == 1, (prefix, found)
    lines[found[0]] = text

replace('2026-11-06更新至v172',
'''2026-11-06更新至v173层流空间壁温与同守恒方向限步。**总目标进行中。独立组合已测、实际4次节点试探恢复；唯一实际候选91.281秒exit134，最终原因未知，完整local0/series1502、最后history56/series1558无场。17项筛选4失败/4待验，不采纳/不续算/不改默认，新增物理资格0。关闭自适应+节点后处理限步支线，下一项先有限评估官方IDEAL_GAS低焓NACA独立基线的物性与能量/空间Tw合同，不以定比热冒充原变比热成功。单工况全面积与科研UI最高，累计热量后置。** 已接受内嵌成熟引擎，无需外部CFD。[当前执行起点](validation/v173/NEXT_DECISION.md) [本轮实现与决定](NACA_SPATIAL_STEP_V173.md) [下一路线依据](validation/v173/REFERENCE_ROUTE_REVIEW.md)''')

idx = next(i for i, line in enumerate(lines) if line.startswith('| v172 |'))
lines.insert(idx + 1, '| v173 | 准确空间Tw+层流限步独立包、真实类/小组合及首帧严格热审核通过；实际4次试探恢复 | 唯一实际91.281秒exit134，最终原因未知、history56无场；17筛选4失败/4待验，支线关闭。小例逐字比较及容量元数据缺失失败保留，新资格0 |')

replace('### B. v172', '### B. v173之后的公开热验证与完整区域推进')
replace('原真实整箭数值运行记录停在v167，',
'''原真实整箭数值运行记录停在v167，空间Tw、冷重启及节点保护组合更新至v173，主默认不变。v164回退基线保留，v165/v167/v172/v173候选均未满足采纳条件，不自动接续。原19752面/面积1.9646292566m²/335726点/705490单元保持。v167完整0/100/199对应累计201/301/400；冷重启不称连续等价。v173真实4次节点恢复后仍异常退出，最终原因未知，残差无收益；关闭该后处理补丁支线。下一项评估官方IDEAL_GAS低焓NACA独立基线，先量化物性/能量/归约差和空间Tw合同，不直接运行或用其解释本次失败。[整箭数值记录](SST_ENERGY_CONSTRAINT_V167.md) [本轮决定](NACA_SPATIAL_STEP_V173.md)''')

replace('**v172当前具体待办：**',
'''**v172历史待办，已由v173执行并关闭候选：** 准确覆盖、真实类小检查、空间Tw节点保护独立包及唯一实际候选均已完成。v172原物性异常/未知调用者、首帧与末历史身份和残差失败保持，不重跑该待办。[历史记录](NACA_ADAPTIVE_DEVELOPMENT_V172.md) [历史下一决策](validation/v172/NEXT_DECISION.md)

**v173当前具体起点：** 原生91.281秒exit134、85输入稳定，完整local0/series1502、末history56/series1558无场。4次node.trial均theta=.5、旧态合法，独立EOS/共同方向/全局坐标通过；首帧严格能量2.004e−13/EOS2.22e−16/热梯度8.132e−16、3226实际Tw差0。17筛选9通过/4失败/4待验；末20rho/E平均log恶化.205700/.183384，原−8未过，不采用、不改默认、不接续，新资格0。完整日志无终止消息，最终原因/类型/调用者未知，不能移用v172。小例PHASE遗漏、原生容量未启用与逐字全场比较失败均保留；实际大任务已修正容量/清点。关闭adaptive+节点后处理限步支线，不再扫CFL/半步/初猜/通量或只为补栈重跑。准确NEWTON_KRYLOV的残差回溯仍TODO，不把它当现成全域保护，不持续自研事务回滚内核。[完整记录](NACA_SPATIAL_STEP_V173.md)

**下一交付限定为官方IDEAL_GAS低焓NACA独立基线的有限适用性检查。** 先≤120秒核对显式gamma/R、实验M/Re_D/T0/D/Tw、相容静态条件、同参考Δh、cp/声速/导热和Nu/√Re归约差，以及准确原生空间Tw/能量合同。量化模型误差并事前明确适用范围，旧非收敛场温区不是实验真实温区。不得沿用OTH的rhoE冷重启或把定cp当原变cp修复；不先构建/大CFD。检查支持后另登记必要最小接口和唯一实际候选的总预算，不支持则关闭该候选。保留高温变比热路线、原热流/−8/y+/网格/外域及实验门槛，低焓层流不授M5/20°、分离、底面和突起资格。[当前下一决策](validation/v173/NEXT_DECISION.md) [准确路线依据](validation/v173/REFERENCE_ROUTE_REVIEW.md)''')

replace('已过物性、原几何、初场、阶段保存、主UI及v169入口',
'''已过物性、原几何、初场、阶段保存、主UI及v169入口/v170网格/v171空间Tw/v172冷重启/v173保护小组合按准确证据复用；失败不追改。近期从B4的官方低焓独立基线适用性检查接续，不重复已过源/原列/边界或节点补丁。该基准只用于指定压缩性层流热交换验证，不授予高超声速湍流、真实20°分离/底面/突起资格；有效改进须返回完整域。停止限步、舍入因子、CFL、通量、初猜及网格器扫描；不得借新版本重开旧失败或将改物理模型包装成同方案修复，每项诊断必须形成有限决策。''')

replace('- **当前下一起点：** 按B4和[下一决策](validation/v172/',
'''- **v172历史下一起点，已由v173完成：** 准确覆盖、独立类、空间Tw限步组合及唯一实际候选已执行并关闭；原失败不追改。[历史下一决策](validation/v172/NEXT_DECISION.md)''')

lines.extend('''
- **v173已测试，随本次提交交付。** 层流空间Tw保护独立包、真实类/小组合通过，真实4次恢复；唯一190334点/179568单元候选91.281秒exit134、父108.640秒、85输入稳定。首帧严格热审通过，最后history56无场；最终原因未知。17筛选9通过/4失败/4待验，候选关闭、默认不变、新资格0。小例元数据缺失、逐字比较及筛选序列化失败保留；无新整箭/MATLAB/主UI，历史科学图已查看。收尾计算38.159秒/180秒，局部路线读源约3秒，开发协调另计。[记录](NACA_SPATIAL_STEP_V173.md)
- **当前下一起点：** 按B4及[下一决策](validation/v173/NEXT_DECISION.md)，关闭当前后处理限步支线；≤120秒评估官方显式IDEAL_GAS低焓NACA独立路线的物性/能量/归约及空间Tw合同。先明确实验锚点和适用域，不用OTH原rhoE、不先构建/大CFD，支持后才登记一次总预算。不是高温变cp失败修复或M5/20°资格；原几何/热流/近壁/−8及失败保持。有效能力回到真实全面积和科研UI，累计热量后置，总目标进行中。
- **v173存储。** 四个逐字完整帧别名硬链接去重122854214字节，8条原路径/SHA及回执保持，观测空闲增加122843136字节（含并发影响）；独有删除0，大场仅本地，无工程包，保护OpenCode/共享历史。[回执](validation/v173/storage_receipt.json)
'''.splitlines())

with p.open('w', encoding='utf-8', newline='\n') as stream:
    stream.write('\n'.join(lines) + '\n')
print('Updated PROJECT_GOALS.md current entry; historical evidence preserved.')
