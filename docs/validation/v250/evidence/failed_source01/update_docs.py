"""模块职责：发布v250真实结果与下一具名决定，保留本地补充，不混入并发修改。"""
from pathlib import Path
import hashlib,json,re,subprocess
R=Path('E:/RocketAeroHeat');B=R/'outputs/v250_validation'
private=json.loads((R/'outputs/v249_validation/concurrent_prefixes.json').read_text(encoding='utf8'))
agent='''**2026-11-16 v250当前起点：原生NK真实108更新正常时间停止，末场与阶段独立通过；线性精度缺口已定位。** 原335726点/705490单元/19752壁/1.9646292565715553m²及全部近壁/下游保持；同v248200初态/CFL重置0.1，首0/1与v249逐位。仅编一个必要只读Newton对象，其他成熟对象复用，构建10.859秒；唯一1698.469秒exit0/峰2491424768B/1CPU，1536MiB及终态容量通过。见docs/NATIVE_NEWTON_DEVELOPMENT_V250.md。

1001完整Jv、0乘积恢复；717实际更新后合法恢复、方向与返回旧态逐位核对，714低于能量下界/3违反p-T联合域。首审核误当全部能量越界1.578秒exit1保持，完整联合域恢复8.281秒通过，未重跑CFD。0/1/50/100/108五场1678630状态、原生末CSV/零时间空间评估七量不变通过。末六式均未达−8，净质量/E残差比0.203454/0.288906未过，非完整边界账；50→100 p/q变化6.46731%/27.19021%，100→108短8步不授稳定。200双窗口未到，不称收益规则停止，新热资格0。

76更新达到原生外层10上限仍未达报告1e−5，首22、最差log−2.162102377、末−4.893137145；这不是独立Jv或唯一根因证明。五阶段22字段/395040原面值、两科研图和前端纯函数通过，实际完成/保存帧/停止原因分开；无MATLAB/共享WebGL/浏览器。6监督终态输入稳定组空，离线3.172/构建10.859/CFD1698.469/审核出口含失败13.109秒，开发归档另计。

**下一次按docs/validation/v250/NEXT_DECISION.md直接具名成熟线性精度闭环。** 优先一次核对并采用现有FGMRES＋块ILU(0)完整组合，准确工厂/内外层精度/flow-SA实际影响与存储先绑定，≤120秒定向准备；当前仅定位候选，未授兼容或资源。必要真实任务加入更新前实际线性精度接受检查，不能把到迭代上限当充分求解；上限须与成熟预条件/预算一起决定，不扫参数或只延长本108。可行后另合同唯一真实分支、复用当前阶段入口，不重复一步/五步/初始化/全表/冻结全域，不冒称精确续算；组合不成立关闭回混合区域。原生报告不授准确Jacobian，不长期自研内核/线搜索或扫CFL通量差分阻尼。当前三CLI关闭，无新构建/查询/CFD准入。工程外缘/流管热入口与科研UI独立最高，原3128前体/5483收缩/后体底面及全分母/近壁/所有科学门槛保持，累计后置、不等大主机。保护29/engineering54/候选/并发本地补充，小证据GitHub大场E盘，无包/独有删除，不新代理/computer-use，总目标active，下方历史。

'''
goal='''**2026-11-16 v250当前起点**：原生NK真实108次更新按时间正常停止，1001次乘积/717次实际恢复、五完整场/末CSV及零时间空间评估独立核对；首步与v249逐位。76次达到线性10上限仍未达报告精度，尚非唯一根因。200正式双窗口未到，六式/净残差/稳定未过，新增合格热面积0。五阶段22字段/395040原面研究值与两科研图已接。按[v250执行点](validation/v250/NEXT_DECISION.md)一次闭合成熟线性精度/块预条件候选，再另合同唯一真实分支；不盲续108，工程外缘/热入口与科研UI独立最高，总目标active。

'''
for name,prefix in [('AGENTS.md',agent),('docs/PROJECT_GOALS.md',goal)]:
    old=subprocess.check_output(['git','show','HEAD:'+name],cwd=str(R)).decode('utf8')
    live=(R/name).read_text(encoding='utf8');local=private[name]['prefix']
    assert live.count(local)==1 and live.replace(local,'',1)==old
    title,body=old.split('\n\n',1);owned=title+'\n\n'+prefix+body
    if name.endswith('PROJECT_GOALS.md'):
        owned,n=re.subn(r'^当前下一项为.*$', '当前下一项为**成熟原生线性精度/块预条件的完整采用决定与可执行入口，以及工程外缘/热入口推进**。v250原生NK已完成108次真实更新、准确阶段/末状态审核和科研出口；按时间停止，正式200双窗口未到，当前场不授物理资格。76次达到线性10上限仍未达报告目标，优先一次核对现有FGMRES＋块ILU(0)的真实flow-SA影响、内外层精度及内存，候选尚未授兼容。按[v250执行点](validation/v250/NEXT_DECISION.md)把更新前精度接受、预算和成熟组合一起决定，不只加时间/上限、不扫参数或重复小控制；必要真实任务另合同。原工程和科研UI独立最高，全部原面/近壁/收缩后体底面责任保持，不等待32GB主机。',owned,count=1,flags=re.M);assert n==1
        needle='| 成熟引擎接口 |';row='| 原生NK真实有限发展 | v250原完整域108更新正常时间停止，1001完整乘积/717实际恢复、首步逐位/末七量不变、5场22字段科研出口通过 | 76次线性上限未到报告容差；六式/净残差/稳定未过，200双窗口未到；首次审核误当全部能量越界失败保持，完整p-T联合域恢复通过，新热资格0 |\n'
        assert needle in owned;owned=owned.replace(needle,row+needle,1)
        oldh='按[v249执行点](validation/v249/NEXT_DECISION.md)与H3协同闭合成熟区域真实边界和反馈；v248未收敛末场和v249一步控制均不作accepted-edge，原3128面工程目标保持。'
        newh='按[v250执行点](validation/v250/NEXT_DECISION.md)与H3协同闭合成熟区域真实边界和反馈；v248未收敛末场、v249一步和v250未收敛108场均不作accepted-edge，原3128面工程目标保持。'
        assert oldh in owned;owned=owned.replace(oldh,newh,1)
        pattern=r'^v234具体多区组合保持关闭，历次源/失败/中断证据保留。.*$'
        replacement='v234具体多区组合保持关闭，历次源/失败/中断证据保留。v244切口827入射界面仍需真实下游反馈；v248 FAS200后窗恶化停止。v250原生NK108正常时间停止、1001乘积/717恢复及末空间独立核对，仍未收敛；76次原生线性上限未达报告精度为新缺口，非唯一根因。下一按[v250决定](validation/v250/NEXT_DECISION.md)直接一次闭合现有FGMRES＋块ILU(0)候选的准确工厂、内外层精度和完整flow-SA资源，成熟组合与更新前精度接受一起落实；可行后另预算唯一真实分支，不重做小控制或盲续108。原−8/y+/热稳定/收支/空间外域/公开验证及全部区域责任保持，匹配公开基准不转授真实面积。'
        owned,n=re.subn(pattern,replacement,owned,count=1,flags=re.M);assert n==1
    dest=B/'owned_source'/name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(owned,encoding='utf8')
    (R/name).write_text(owned.replace(title+'\n\n',title+'\n\n'+local,1),encoding='utf8')
(B/'concurrent_prefixes.json').write_text(json.dumps(private,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
p=R/'README.md';s=p.read_text(encoding='utf8');marker='## 最新状态与交付（2026-11-16）\n\n';assert s.count(marker)==1
s=s.replace('最新[v249','前版[v249',1)
intro='''最新[v250原生Newton-Krylov真实有限发展](docs/NATIVE_NEWTON_DEVELOPMENT_V250.md)：原完整箭体正常完成108次更新、按时间预算停止；1001完整乘积/717实际恢复与末场独立核对。**未收敛、无新增合格热面积，200双窗口未执行。** 76次线性求解达到上限仍未达报告目标，下一闭合成熟线性精度/预条件，不能盲续108。本机约28.3分钟/峰2.32GiB；五阶段22字段已接[首末云图](docs/validation/v250/evidence/workbench01/native_newton_development.png)和[基准趋势](docs/validation/v250/evidence/workbench01/nk_residual_trends.png)，本地入口`outputs/v250_validation/workbench01/index.html`。

[下一执行点](docs/validation/v250/NEXT_DECISION.md)：现有FGMRES＋块ILU(0)候选的完整采用决定、真实精度接受与资源合同；候选尚未授兼容，工程外缘/热入口和科研UI独立最高，全部责任/门槛不减，总目标active。

'''
p.write_text(s.replace(marker,marker+intro,1),encoding='utf8')
p=R/'.gitattributes';s=p.read_text(encoding='utf8');assert 'docs/validation/v250/evidence/**' not in s
p.write_text(s+'\ndocs/validation/v250/evidence/** -text -eol -whitespace\n',encoding='utf8')
for name in ('runNewtonDevelopmentV250','auditNewtonDevelopmentV250','exportNewtonDevelopmentV250'):
    p=R/'scripts'/(name+'.py');s=p.read_text(encoding='utf8');assert s.count('CLOSED=False')==1;p.write_text(s.replace('CLOSED=False','CLOSED=True'),encoding='utf8')
# 仅绑定下一候选已经定位的准确本地文本；不声明运行、内存或方法采用通过。
sources=[]
for rel,words in [('Common/include/linear_algebra/CPreconditioner.hpp',('class CILUPreconditioner','case ILU:', 'BuildILUPreconditioner')),
                  ('Common/src/CConfig.cpp',('LINEAR_SOLVER_ILU_FILL_IN',))]:
    p=R/'outputs/v220_validation/source'/rel;raw=p.read_bytes();lines=raw.decode('utf8').splitlines()
    hits=[dict(line=i+1,text=line.strip()) for i,line in enumerate(lines) if any(w in line for w in words)];assert hits
    sources.append(dict(path=str(p),sha256=hashlib.sha256(raw).hexdigest(),matches=hits))
(B/'next_method_sources.json').write_text(json.dumps(dict(candidate='native FGMRES with block ILU(0)',adopted=False,executed=False,sources=sources),indent=2)+'\n',encoding='utf8')
print('v250 execution points updated; three task CLIs closed; private supplements preserved.')
