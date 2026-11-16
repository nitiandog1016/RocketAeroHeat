"""模块职责：合并本轮执行点，保留并发公开基准补充且不将其纳入GitHub交付。"""
from pathlib import Path
import subprocess,json,hashlib
ROOT=Path('E:/RocketAeroHeat');BASE=ROOT/'outputs/v244_validation'
private={};owned={}
for name in ('AGENTS.md','docs/PROJECT_GOALS.md'):
 p=ROOT/name;current=p.read_text(encoding='utf8')
 original=subprocess.check_output(['git','show','HEAD:'+name],cwd=str(ROOT)).decode('utf8').replace('\r\n','\n')
 start=current.index('\n\n')+2;stop=current.index('**2026-11-16 v243当前起点')
 prefix=current[start:stop];assert prefix.startswith('**2026-11-16') and '用户补充' in prefix
 assert current[:start]+current[stop:]==original,('UNEXPECTED_CONCURRENT_EDIT',name)
 private[name]=dict(prefix=prefix,initialSHA256=hashlib.sha256(p.read_bytes()).hexdigest())
 owned[name]=original

entry='''**2026-11-16 v244当前起点：真实黏性区域包与827入射界面已核对，直接切口生产拒绝。** 保留原160304母单元（141888棱柱/18416tet）、76550点、4434壁/2186外场/1984接口；完整边界/文本回读及原v236300状态2602700标量逐位同。独立Cv积分/Brent与法向对照max5.6704e−14，0/1/4/5入射分别1157/793/26/8。此为研究态原多边形特征，非原生双面通量/分离证明/黏性扩散闭合。见docs/VISCOUS_REGIONAL_BOUNDARY_V244.md。

真实壁接触46面，原报告按严格零壁距误计14保留；最高2.4533e−18m舍入量，修正函数与独立拓扑实际通过，不裁壁距/状态。首声速字段失败2.734秒、恢复3.578；首审核重复解压母面表120.188秒超时exit15，缓存恢复2.422秒，未重跑场。6监督终态稳定组空，离线7.093/审核出口124.141秒，峰398389248B，归档另计。完整箭体/接口科研图已看，无新CFD/构建/查询/MATLAB/共享WebGL/浏览器/热资格。

采用保留真实下游的同一物理域内部耦合，1984切面只作观察索引；原705490单元/19752壁全部保留，旧2.153GB仅资源参考。相态15静压配置仍禁；成熟Riemann有secondary不等于有真实外压，SA外推不等于双向黏性供体，MZ17组合仍关。原生通信/全局壁距ADT来源已绑定，未安装MPI、不当收敛修复。当前几何包无可盲启物理配置，不拿旧300/研究压力作边界。

**下一次按docs/validation/v244/NEXT_DECISION.md直接具名成熟FAS整体相容控制。** 已绑定原生限制/强迫/递归与延拓源，原细网格作为生产层保持；先一次≤120秒确认相态粗层联合域、无滑移投影、SA/secondary与能量零点及观察器层身份，再另预算一个固定整体控制。可用匹配物性/壁温/实测热的公开标准模型，复用TMR/Run14最新记录，不自动重启或新建无界支线；真实20°箭体仍最终全表面验收。FAS未实测/未授生产，不称低频误差唯一根因；通过才唯一具名真实发展，不盲续旧300/800/400/199/92，不扫层数/CFL/通量/日程或写通用内核。当前四CLI关闭，无新构建/查询/CFD准入。原面/近壁/全分母/门槛及所有区域责任保持；工程/UI独立最高、累计后置，保护29/engineering54/候选/并发工作，小证据GitHub大场E盘，无包/独有删除，不新代理/computer-use，总目标active，下方历史。

'''
s=owned['AGENTS.md'];i=s.index('\n\n')+2;owned['AGENTS.md']=s[:i]+entry+s[i:]
s=owned['docs/PROJECT_GOALS.md'];i=s.index('**2026-11-16 v243当前起点**');j=s.index('\n\n',i)
s=s[:i]+'''**2026-11-16 v244当前起点**：真实前体160304母单元/76550点/4434壁区域已导出，2602700原状态标量逐位保持；1984接口中827个需要下游信息，不能直接设出口。原壁接触统计14→46已独立纠正，失败保留。采用保留下游的同一物理域内部耦合，当前切口仅观察索引；无新CFD/热资格。按[v244执行点](validation/v244/NEXT_DECISION.md)定向检验成熟FAS与相态Air/SA/热壁/能量零点的整体相容性，尚未授生产。可用有匹配公开实测数据的标准模型诊断，真实20°箭体保持阶段集成和最终全面积验收；工程/科研UI独立最高，总目标active。'''+s[j:]
i=s.index('当前下一项为**');j=s.index('\n\n',i)
s=s[:i]+'''当前下一项为**保留真实反馈的成熟求解组合及工程入口闭环**。v244实际表明当前黏性切口有827个入射界面，不能直接用Euler全出射、p∞、零梯度、研究压力或给定SA填边界。原下游在同一物理域内部求解，前体/下游索引保留；MZ17失败组合仍关闭。下一具名FAS聚合多重网格候选先做完整层间相容控制，保留原细网格和近壁资产；不将FAS/MPI当已证收敛修复或另起自研通量Jacobian。可优先匹配公开标准算例，但不重启旧失败或削减真实箭体最终范围。收缩旧5483面、后体/底面等全部责任保留，工程/UI不因单项CFD门槛停工，不等待32GB主机。'''+s[j:]
anchor='| 成熟引擎接口 |';assert s.count(anchor)==1
s=s.replace(anchor,'| 真实黏性区域与反馈 | v244原160304母单元/76550点/4434壁导出，1984接口、2602700状态标量及独立Brent/法向通过，科研图已核对 | 827界面需下游输入；直接切口生产拒绝，旧研究场不作外压；采用单物理域内部反馈，FAS层间相容性待整体控制，新热资格0 |\n'+anchor)
s=s.replace('[v243执行点](validation/v243/NEXT_DECISION.md)','[v244执行点](validation/v244/NEXT_DECISION.md)')
s=s.replace('[v243决定](validation/v243/NEXT_DECISION.md)','[v244决定](validation/v244/NEXT_DECISION.md)')
i=s.index('每次先读AGENTS.md与本文件第4节H1—H5');j=s.index('\n\n',i)
s=s[:i]+'''每次先读AGENTS.md与本文件第4节H1—H5，以[v244执行点](validation/v244/NEXT_DECISION.md)覆盖历史下一步。真实区域包、入射需求及单物理域反馈决定已完成；直接一次具名FAS完整层间相容控制，先准确来源/方案再另预算必要构建与原生整体调用。公开匹配基准可以承担日常诊断，真实箭体仍为阶段集成及最终全面积验收。当前四固定CLI关闭、无新CFD/构建/查询准入；不盲续旧场、不扫参数或写通用内核，不以研究场当accepted-edge。新实际原面结果接科研工作台，工程/全面积/UI独立最高，累计后置。'''+s[j:]
owned['docs/PROJECT_GOALS.md']=s
for name,s in owned.items():
 p=BASE/'owned_source'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(s.encode('utf8'))
 i=s.index('\n\n')+2;actual=s[:i]+private[name]['prefix']+s[i:]
 (ROOT/name).write_bytes(actual.encode('utf8'))
(BASE/'concurrent_prefixes.json').write_text(json.dumps(private,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
p=ROOT/'README.md';s=p.read_text(encoding='utf8');anchor='## 最新状态与交付（2026-11-16）\n\n'
addition='''最新[v244真实黏性区域与边界](docs/VISCOUS_REGIONAL_BOUNDARY_V244.md)：已导出160304原单元/76550点/4434原壁的区域几何与状态映射。1984接口中827个需要下游信息，当前直接切口生产拒绝；保留真实下游的单物理域内部反馈。没有新CFD或新增合格热面积。

[完整箭体与接口科研图](docs/validation/v244/evidence/visual01/viscous_region_boundary.png)和独立原面/状态核对通过，壁面接触统计修正14→46，原失败保留。[下一执行点](docs/validation/v244/NEXT_DECISION.md)为成熟FAS与相态Air/热壁/SA的整体相容控制，原细网格保持；匹配公开基准可作日常诊断，全面积/UI仍独立最高。

'''
assert s.count(anchor)==1;s=s.replace(anchor,anchor+addition).replace('最新[v243前体空间对照]','旧[v243前体空间对照]',1);p.write_bytes(s.encode('utf8'))
p=ROOT/'.gitattributes';s=p.read_text(encoding='utf8');assert 'docs/validation/v244/evidence/' not in s
p.write_bytes((s.rstrip()+'\ndocs/validation/v244/evidence/** -text -eol -whitespace\n').encode('utf8'))
for name in ('prepareViscousRegionV244','auditViscousRegionV244','reviewViscousCouplingV244','renderViscousRegionV244'):
 p=ROOT/'scripts'/(name+'.py');s=p.read_text(encoding='utf8');assert s.count('CLOSED=False')==1
 p.write_bytes(s.replace('CLOSED=False','CLOSED=True').encode('utf8'))
print('Updated v244 execution points; preserved concurrent local-only prefixes; four fixed CLIs closed.')
