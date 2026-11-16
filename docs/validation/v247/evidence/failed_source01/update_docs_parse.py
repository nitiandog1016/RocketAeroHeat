"""模块职责：更新当前执行起点，分别保存公开文档和含用户本地补充的工作副本。"""
from pathlib import Path
import json,re,subprocess
R=Path('E:/RocketAeroHeat');B=R/'outputs/v247_validation'
private=json.loads((R/'outputs/v246_validation/concurrent_prefixes.json').read_text(encoding='utf8'))
agent='''**2026-11-16 v247当前起点：真实FAS完成60更新，61合法原生恢复后被观察器中止。** 原335726/44484点、705490单元、19752壁/1.9646292565715553m²及全部近壁/下游保持。唯一270.468秒exit86/峰1968996352B/1CPU，无资源触发。960完整阶段＋61的10部分；完整场0/1，壁0/1/50，无完整60/61/原生末CSV。见docs/FAS_BOUNDED_DEVELOPMENT_V247.md。

细全局329559/内部42372、壁距33.5183µm，prolong_out之后试探物理e−257956.386低于60K域下界1953.749J/kg；旧态74.1798K合法，returned五量逐位old。原生恢复成功，但观察器按合同抛错；不是原生不可恢复崩溃/200收益停止。old不是延拓前态，真实方向未保存，不认唯一根因。671452保存状态/26Brent、原面237024值/14字段及科研图通过，热资格0。初0十一列逐位；初1最大缩放差7.575e−16、逐位false，机制未唯一定位。终态空间入口已编译未到达，不授末残差或保存资格。

准备chrono绑定、初1严格逐位审核、壁表通配符及外部源额外换行失败全部保留；恢复审核通过，无CFD重跑。只薄构建7.234秒；11监督终态输入稳定组空，离线3.672/CFD270.468/审核出口14.265秒，归档与定向源获取另记。PETSc FAS粗修正线搜索/SECANT/BT三个公开blob核对，仅具名备选，未安装或实测。

**下一次按docs/validation/v247/NEXT_DECISION.md直接保留成熟合法旧态恢复的完整真实FAS发展。** 只在recovered=true回调记录而不主动抛错，非法旧态仍拒绝；不裁态/放域/改CFL通量阻尼或新增自研线搜索。必要真实任务一起保存延拓前态/方向、实际返回后flow-SA及体积质量/能量扰动；不拿trial−old充粗修正，不将恢复授守恒。先另登记均匀新初场200/1800秒/1CPU/4GiB/1536MiB，保存0/1/每50与正常终点，原收益规则/终态原生空间评估和所有科学门槛保持。当前四CLI关闭，无新构建/CFD准入。不续不存在60/旧5、199、300、800、92，不再五步/全表/旧故障重跑或冻结整域。工程/科研UI独立最高、原所有区域与近壁责任不减，累计后置；保护29/engineering54/候选/并发补充，小证据GitHub大场E盘，无包/独有删除，不新代理/computer-use，总目标active，下方历史。

'''
goal='''**2026-11-16 v247当前起点**：原完整真实域FAS完成60更新，第61次试探越域后原生恢复合法旧态，观察器按合同中止；不是不可恢复崩溃或收益停止。原面0/1/50三阶段14字段/237024值已接工作台，新增合格热面积0。按[v247执行点](validation/v247/NEXT_DECISION.md)保留成熟恢复、完整记录接受与状态扰动账，并在唯一新有限发展中补齐正常终态残差/输出；不等待大主机、不重复小壳/五步或扫描参数。工程/科研UI独立最高，总目标active。

'''
for name,prefix in [('AGENTS.md',agent),('docs/PROJECT_GOALS.md',goal)]:
 old=subprocess.check_output(['git','show','HEAD:'+name],cwd=str(R)).decode('utf8')
 live=(R/name).read_text(encoding='utf8');local=private[name]['prefix']
 assert live.count(local)==1 and live.replace(local,'',1)==old
 title,body=old.split('\n\n',1);owned=title+'\n\n'+prefix+body
 if name.endswith('PROJECT_GOALS.md'):
  owned=re.sub(r'^当前下一项为.*$', '当前下一项为**保留原生合法恢复的固定真实FAS完整发展及工程入口闭环**。v247原生已恢复合法旧态，额外观察器中止使200收益与终态空间评估未到达；应先使观察器被动记录已成功恢复，不再把这一事件误判为不可恢复崩溃。实际延拓前态/方向及完整flow-SA接受必须同次记录，保持原有科学门槛；旧态恢复可能带来守恒扰动，不预授热资格。PETSc粗修正线搜索只作为已核对源码的备选，未安装或验证。收缩5483面、后体/底面和全部区域责任保持，工程/UI不等待单项CFD门槛或32GB主机。',owned,count=1,flags=re.M)
 out=B/'owned_source'/name;out.parent.mkdir(parents=True,exist_ok=True);out.write_text(owned,encoding='utf8')
 (R/name).write_text(owned.replace(title+'\n\n',title+'\n\n'+local,1),encoding='utf8')
(B/'concurrent_prefixes.json').write_text(json.dumps(private,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
p=R/'README.md';text=p.read_text(encoding='utf8');marker='## 最新状态与交付（2026-11-16）\n\n'
assert text.count(marker)==1;text=text.replace('最新[v246','前版[v246',1)
intro='''最新[v247真实FAS有限发展](docs/FAS_BOUNDED_DEVELOPMENT_V247.md)：原完整网格完成60更新，第61次原生已合法恢复旧态，但观察器按合同中止；**未收敛、无新增合格热面积**。不是资源不足或原生不可恢复崩溃。三阶段0/1/50、14字段原面研究已接[科研图](docs/validation/v247/evidence/workbench01/rocket_fas_control.png)，本地入口`outputs/v247_validation/workbench01/index.html`。

[下一执行点](docs/validation/v247/NEXT_DECISION.md)：保留成熟原生合法恢复，记录完整接受/扰动账并完成有界发展与正常终态原生残差；不扫参数或重跑旧故障。原面/近壁/全部区域与科学门槛不变，工程/UI最高，总目标active。

'''
p.write_text(text.replace(marker,marker+intro,1),encoding='utf8')
p=R/'.gitattributes';text=p.read_text(encoding='utf8');assert 'docs/validation/v247/evidence/**' not in text;p.write_text(text+'\ndocs/validation/v247/evidence/** -text -eol -whitespace\n',encoding='utf8')
for name in ('runFASDevelopmentV247','auditFASDevelopmentV247','exportFASDevelopmentV247','reviewFASRecoveryV247'):
 p=R/'scripts'/(name+'.py');text=p.read_text(encoding='utf8');assert text.count('CLOSED=False')==1;p.write_text(text.replace('CLOSED=False','CLOSED=True'),encoding='utf8')
