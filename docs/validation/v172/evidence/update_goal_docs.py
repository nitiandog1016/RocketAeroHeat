"""模块用途：只更新目标中的当前摘要和执行入口，保留历史失败及门槛。"""
from pathlib import Path

root=Path(__file__).resolve().parents[2]
p=root/'docs/PROJECT_GOALS.md'
s=p.read_text(encoding='utf-8')
lines=s.splitlines()
assert lines[2].startswith('2026-11-06更新至v171')
lines[2]='2026-11-06更新至v172同场冷重启与自适应发展评估。**总目标进行中。接口及真实反馈已测，唯一实际候选79.782秒exit134，完整local0/series1502、最后history45/series1547无场；正常与残差收益失败，不采纳/不续算/不改默认，新增物理资格0。下一项有限核对原始态、重构/边界物性查询保护，评估层流+空间Tw同守恒方向限步/拒绝；不扫描CFL上限，不以未知异常调用者作因果结论。单工况全面积与科研UI最高，累计热量后置。** 已接受内嵌成熟引擎，无需外部CFD。[当前执行起点](validation/v172/NEXT_DECISION.md) [本轮实现与决定](NACA_ADAPTIVE_DEVELOPMENT_V172.md) [路线复审](reviews/TECHNICAL_ROUTE_REVIEW_V168.md)'
s='\n'.join(lines)+'\n'
old='### B. v171之后的公开热验证与完整区域推进'
assert s.count(old)==1
s=s.replace(old,'### B. v172之后的公开热验证与完整区域推进')
prefix='原真实整箭数值运行记录停在v167，空间Tw运行与实际三维热发展更新至v171，'
line=next(x for x in s.splitlines() if x.startswith(prefix))
s=s.replace(line,'原真实整箭数值运行记录停在v167，空间Tw及同场冷重启/实际发展更新至v172，主默认不变。v164回退基线保留，v165/v167及本次自适应候选均未满足采纳条件，不自动接续。原19752面/面积1.9646292566m²/335726点/705490单元保持。v167完整0/100/199对应累计201/301/400；冷重启不称连续等价。v171空间Tw和原段功率已通过，v172实际控制器增长/回退已执行但仍发生物性温区异常。当前先核清保护覆盖，形成层流+空间Tw可接受域限步/拒绝的小型组合合同，不能假定本次具体调用者或仅扫描CFL。[整箭数值记录](SST_ENERGY_CONSTRAINT_V167.md) [本轮决定](NACA_ADAPTIVE_DEVELOPMENT_V172.md)')
row=next(x for x in s.splitlines() if x.startswith('| v171 |'))
s=s.replace(row,row+'\n| v172 | 具名冷重启/自适应实际反馈、小例与来源拒绝、原列诊断通过；仅首完整场严格热审核通过 | 唯一实际79.782秒exit134，最后history45无场；正常/收益失败，候选关闭、默认不改；末场壁流体温度/y+/热稳定未知，新资格0 |')
oldline=next(x for x in s.splitlines() if x.startswith('**v171当前具体待办：**'))
s=s.replace(oldline,oldline.replace('**v171当前具体待办：**','**v171历史待办，已由v172完成有限核验和唯一候选：**').replace('[当前下一决策]','[历史下一决策]')+'\n\n**v172当前具体待办：** 唯一候选物性反解温区异常exit134，完整local0/series1502及最后history45/series1547分列；首帧热审核、3226施加Tw通过，末场壁流体温度及对给定Tw偏差/y+/p-q稳定未知。真实反馈已执行，9/46行线性未达1e−5，末20rho/E log分别恶化.176709/.155319，正常/收益失败，不采纳、不续算。残缺栈未定位查询模式/节点/越界方向，不能移用v163。下一项≤120秒有限准确源覆盖表及独立算子，评估同守恒方向可接受域限步/拒绝的层流+空间Tw组合，保留v165旧失败，不直接移植称根因修复。最小拒绝上下文随实际改动保存，不仅为补栈重跑；组合合同/覆盖明确才另行登记唯一候选和预算，否则停止支线。原−8/物性温区/热流/近壁/几何及失败保持，有效能力返回真实20°完整域和科研UI。[完整记录](NACA_ADAPTIVE_DEVELOPMENT_V172.md) [当前下一决策](validation/v172/NEXT_DECISION.md)')
oldline=next(x for x in s.splitlines() if x.startswith('已过物性、原几何、初场、阶段保存、主UI、公开输入及v169'))
s=s.replace(oldline,'已过物性、原几何、初场、阶段保存、主UI及v169入口/v170网格/v171空间Tw/v172冷重启与反馈按既有证据复用。近期从B4的保护覆盖与具名层流组合接续，不重复已过源/原列/边界或小步诊断。该基准只用于指定压缩性层流热交换验证，不授予高超声速湍流、真实20°分离/底面/突起资格；有效改进须返回完整域。停止同类限步参数、舍入因子、CFL上下界、通量、初猜及网格器扫描；新组合须明确覆盖差异和收益，不借新版本重开旧失败候选，每项诊断必须形成有限决策。')
oldline=next(x for x in s.splitlines() if x.startswith('- **当前下一起点：** 按B4和[下一决策](validation/v171'))
s=s.replace(oldline,oldline.replace('**当前下一起点：**','**v171历史下一起点，已由v172执行并关闭候选：**'))
s+='\n- **v172已测试，随本次提交交付。** 具名冷重启/1正向18拒绝/旧23回归、唯一8更新小例及准确CFL反馈通过；原列诊断不是网格资格。唯一实际79.782秒exit134、父98.390秒、74输入稳定；完整local0/series1502，最后history45/series1547无场。首帧热审核通过，正常/残差收益失败，终态壁流体温度/y+/p-q稳定未知，不采纳/不改默认/不接续，新增物理资格0。收尾计算29.719秒/180秒；开发协调另计。准确准备/验证失败保留，无新整箭/MATLAB/主UI。[记录](NACA_ADAPTIVE_DEVELOPMENT_V172.md)\n'
s+='- **当前下一起点：** 按B4和[下一决策](validation/v172/NEXT_DECISION.md)，先≤120秒有限准确源覆盖及独立算子，列清原始态、重构/边界反解保护，评估层流+空间Tw同守恒方向可接受域限步/拒绝。v165旧失败与本次未知调用者保持；最小拒绝上下文随实际改动，不仅补栈重跑。组合合同和保护范围明确才另定唯一候选总预算，否则停止支线；不扫CFL/通量/初猜或全编/无限诊断。有效能力回到原19752面真实20°全面积/分离再附/底面/突起与科研UI，累计热量后置，总目标继续。\n'
s+='- **v172存储。** 仅两份逐字末帧别名原子硬链接去重，4条路径/SHA保持，重复逻辑122482732字节，观测空闲增加122486784字节（含系统并发影响）；独有删除0，大场本地，无工程包，保护OpenCode/共享历史。[回执](validation/v172/storage_receipt.json)\n'
with p.open('w',encoding='utf-8',newline='\n') as stream:
    stream.write(s)

# 模块二：避免将给定Tw与异常时未知的壁流体温度混为一谈。
for name in ['docs/NACA_ADAPTIVE_DEVELOPMENT_V172.md','docs/validation/v172/README.md','docs/validation/v172/NEXT_DECISION.md']:
    q=root/name;t=q.read_text(encoding='utf-8')
    t=t.replace('末场Tw、y+','末场壁流体温度及对给定Tw偏差、y+')
    t=t.replace('终态Tw误差','终态壁流体温度对给定Tw偏差')
    t=t.replace('终态Tw、y+','终态壁流体温度及对给定Tw偏差、y+')
    t=t.replace('首帧Tw/y+不能替代终态','首帧壁流体温度偏差/y+不能替代终态')
    with q.open('w',encoding='utf-8',newline='\n') as stream:
        stream.write(t)
print('Updated current goal entries; historical evidence retained.')
