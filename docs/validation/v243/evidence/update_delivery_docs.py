"""模块职责：按已保存实际证据更新当前执行点，保留历史目标与并发修改。"""
from pathlib import Path
import json
ROOT=Path('E:/RocketAeroHeat')
entry='''**2026-11-16 v243当前起点：供体/插值已分解，一次前体细化真实800更新完成但未改善外缘。** 四管210经过单元/6停止面表明节点H已有主要漂移；具名scikit-fem RT/Piola仅调研，未安装运行，拒绝质量投影直接供热。一次Gmsh细化原结果8单元SICN<0.1保留；另登记4母单元中心16tet有限重接，其余162200tet/全部坐标边界保持，独立最低SICN0.102563。34956点/162216tet/17736计算壁面按面积回4434原面，原全分母/近壁保持。见docs/FOREBODY_SPATIAL_COMPARISON_V243.md。

同相态Euler16/Roe/原二阶/CFL，新相容均匀初态；成熟对象零重编，薄构建5.875秒。唯一800更新287.359秒exit0/峰264503296B/1CPU，正常预算终点非收益停止；0/1/每50/800及原生末CSV逐位同。3022出口点/2417600调用/12088000分量同，最低Mn2.11143。末预处理拒绝最大519/末3，准确ID与全部子阶段未知。五式末logRMS−5.721/−2.588/−3.505/−3.125/0.747均未达−8，净质量/E1.305e−7/4.586e−7未达1e−8；压力稳定通过，H均差1.6823%/最大5.2177%。粗细原面压力差13.2149%，两支均未稳态，不能分离空间与迭代误差或授网格无关；新热面积/资格0。

3阶段21字段/79812原面研究值及前端纯函数通过，原面聚合/共同色标科研图已看；q与固体Tw未知，给定Tw单列，无MATLAB/共享WebGL/浏览器。首源依赖查找0.968秒与首Mach父面聚合1.859秒失败保留，恢复通过未重跑CFD；原Gmsh科学拒绝另存。14监督终态输入稳定组空，离线5.968/网格1.579/构建5.875/CFD287.359/审核出口21.219秒，归档另计；早期256MiB/CFD前新增384MiB合同分开。

**下一次按docs/validation/v243/NEXT_DECISION.md直接真实前体单区成熟黏性边界与反馈方案。** 复用v233前160304母单元/4434壁/1984接口索引、原鼻部发展入口和已有单区相态Air/SA-neg-WITHFT2；索引非已可运行域，MZ17失败组合仍关。一次≤120秒定向明确局部法向特征、真实外压/下游供体、黏性/SA/全壁距与必要反馈，成立即区域包/原生整体控制/另预算唯一发展。无滑移层内亚声速，不能照搬Euler全出射或人工p∞/零梯度/研究压力充边界；缺供体则具名包含下游/反馈，不扫截断/参数或长期自研通量Jacobian。当前八固定CLI关闭，无新构建/查询/CFD准入；不盲续本800/旧400、200、199、92或叠加细化/物性小壳/冻结整域。3128目标、收缩5483/后体/底面全部责任和科学门槛保持，工程/UI独立最高、累计后置；保护29/engineering54/候选/并发web与文档，小证据GitHub大场E盘，无包/独有删除，不新代理/computer-use，总目标active，下方历史。

'''
p=ROOT/'AGENTS.md';s=p.read_text(encoding='utf8');assert 'v243当前起点' not in s
head='# RocketAeroHeat 项目执行约定\n\n';assert s.startswith(head);p.write_text(head+entry+s[len(head):],encoding='utf8')

p=ROOT/'docs/PROJECT_GOALS.md';s=p.read_text(encoding='utf8');start=s.index('**2026-11-16 v242当前起点**');end=s.index('\n\n',start)
s=s[:start]+'''**2026-11-16 v243当前起点**：已分解真实四管供体H与插值偏差，拒绝质量投影直接供热；一次前体细化经有限局部拓扑修正通过几何，162216tet/34956点，真实800更新正常完成（287.359秒/约265MB）。五式残差和净收支未过，H偏差未改善，粗细原面压力差13.2149%；两支均未完整稳态，不授网格无关或可信外缘。3阶段21字段/79812原面研究值与科研对照图已接，新热面积/资格0。按[v243执行点](validation/v243/NEXT_DECISION.md)直接闭合真实前体单区成熟黏性边界及下游反馈，不盲续800或继续Euler细化。全面积/科研UI独立最高，总目标active。'''+s[end:]
start=s.index('当前下一项为**');end=s.index('\n\n',start)
s=s[:start]+'''当前下一项为**H1/H2/H3真实区域入口与成熟局部黏性反馈闭环**。复用v233前体160304母单元/4434壁面/1984接口索引及已有单区相态Air/SA-neg-WITHFT2，先明确真实入口、截断特征、黏性/湍流条件、全壁距和下游入射信息。鼻部发展真实剖面；无滑移层内亚声速，不能把Euler全出射或人工p∞、零梯度、未合格研究压力当作黏性出口依据。一次具名成熟边界/反馈决定后直接形成区域输入及整体原生控制；缺供体则具名包含下游或反馈，不长期自研通用内核。v234 MZ17失败组合仍关闭。不得继续细化扫描、盲续本800/旧400、移种子或冻结H/S补解；−8及空间/外域/热资格门槛保持。收缩旧5483面、后体/底面等全部责任保留，工程/UI不因单项CFD门槛停工，也不等待32GB主机。'''+s[end:]
anchor='| 成熟引擎接口 |'
assert s.count(anchor)==1
s=s.replace(anchor,'| 前体供体分解与空间对照 | v243四管210实际单元分解；一次细化及4母单元有限重接，800真实更新、18完整帧/末CSV、原父面聚合和3阶段21字段通过 | 五式/净收支失败、H未改善、粗细压力差13.2149%；未分离空间与迭代误差，RT直接供热拒绝；新热面积/资格0，下一局部黏性真实边界/反馈 |\n'+anchor)
p.write_text(s,encoding='utf8')

p=ROOT/'README.md';s=p.read_text(encoding='utf8');anchor='## 最新状态与交付（2026-11-16）\n\n';assert s.count(anchor)==1
addition='''最新[v243前体空间对照](docs/FOREBODY_SPATIAL_COMPARISON_V243.md)：完成真实供体/插值分析与一次细网格800更新（287.359秒/约265MB）。几何及原面聚合通过，压力稳定，但残差、净收支及总焓偏差未改善；粗细原面压力差13.2149%，**不能据此认定网格无关或可信外缘，没有新增合格壁热面积**。

本地`outputs/v243_validation/workbench01/index.html`提供0/400/800三阶段21字段，79812研究数值及前端纯函数通过。[共同色标科研图](docs/validation/v243/evidence/visual01/spatial_comparison.png)同时显示压力、H偏差、残差和原面差值；未做浏览器交互验收。[下一执行点](docs/validation/v243/NEXT_DECISION.md)是有真实入口/出口及下游反馈依据的成熟局部黏性区域，工程/UI独立推进，不盲续800或叠加Euler细化。

'''
s=s.replace(anchor,anchor+addition).replace('最新[v242保持原壁','旧[v242保持原壁',1);p.write_text(s,encoding='utf8')

p=ROOT/'docs/validation/v243/research/sources.json';r=json.loads(p.read_text(encoding='utf8'))
assert not any(x['file']=='element_hdiv.py' for x in r['sources'])
r['sources'].append(dict(file='element_hdiv.py',url='https://raw.githubusercontent.com/kinnala/scikit-fem/f598e787feadcc63f077e41bf95d93a784f6c172/skfem/element/element_hdiv.py'))
p.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf8')

p=ROOT/'.gitattributes';s=p.read_text(encoding='utf8');assert 'docs/validation/v243/evidence/' not in s
p.write_text(s.rstrip()+'\ndocs/validation/v243/evidence/** -text -eol -whitespace\n',encoding='utf8')
names=('inspectForebodyTransportV243','refineForebodyV243','reviewRefinementDiagonalV243','reconnectRefinementV243','runForebodySpatialV243','auditForebodySpatialV243','exportForebodySpatialV243','renderForebodySpatialV243')
for name in names:
 p=ROOT/'scripts'/(name+'.py');s=p.read_text(encoding='utf8');assert s.count('CLOSED=False')==1
 p.write_text(s.replace('CLOSED=False','CLOSED=True'),encoding='utf8')
print('Updated current goals, history, sources and eight fixed CLI guards; no numerical work executed.')
