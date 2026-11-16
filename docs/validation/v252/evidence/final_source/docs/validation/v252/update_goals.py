"""本轮一次性文档更新：发布本轮内容，保留已有本地不发布的公开基准补充。"""
from pathlib import Path
import subprocess,json
ROOT=Path(__file__).resolve().parents[3];BASE=ROOT/'outputs/v252_validation'

def main():
    additions={
    'AGENTS.md':'''**2026-11-16 v252当前起点：同相态非相似剖面三分支与独立积分账通过，原面/剖面科研页已接。** 采用TP-3271层流Levy–Lees/沿xi后退差分＋成熟bvp4c，旧v61保持冻结。复用已有驻点，31/61/31站共120个新沿程站；独立物理T复算30站最大1.968e−6，步长减半q差9.803e−5，40°同终点Tw历史差2.70158%。原始积分偏差最大.234767%/扣已知后退截断后2.508e−5，均过1%。见docs/PHASE_NONSIMILAR_PROFILES_V252.md。

原19752面/1.9646292565715553m²保持，只原1306球帽面/.1790893%候选，基准功率862.014876W；新增合格面积0。34字段/671568值、三链剖面滑块/CSV与纯前端通过，两科研图已看；首3D刻度重叠保持、正投影恢复。无新CFD/构建/共享WebGL/浏览器。仅常H/同熵/β0—1附着层流，无真实三维横流/分离/实验资格，profileReadyForCFD=false。

原MATLAB三链报告完整，但父监督终态缺失、准确进程已无；原退出码/耗时/峰值/全程稳定未知。恢复观察输入SHA同/八结果稳定只属当前事实，原300秒保守扣足不补跑。另一次具名软件修复回归32.047秒exit0/峰1014579200B：全部积分后才授available，映射源元数据核对；第2站六量差0/原面与5拒绝通过。JSON重新序列化身份另列，不追改原摘要。监督准备.265秒/审核出口13.486秒、各已有回执稳定组空；原缺口和locale警告保留。

**下一次按docs/validation/v252/NEXT_DECISION.md直接鼻部之后3128目标面的真实外缘—熵层—横向流管入口。** 定向复用CR-112180/v96/v241/v243/v244，一次采用或拒绝具名成熟闭合并落可执行生产/区域入口与有区分力验证；不能又整版鼻部/平板/物性/BVP，不能固定驻点熵铺全后体或删零攻角守卫。研究源资格传播，5483收缩/后体底面仍需真实成熟黏性反馈；无人工p∞/SA常数，不续旧整箭/ILU参数支线，不等大主机。当前四CLI关闭、新任务另预算；原面/分母/近壁/科学门槛保持，全面积/UI最高、累计后置。保护29/engineering54/候选/并发，小证据GitHub大场E盘，无包/独有删除，不新代理/computer-use，总目标active，下方历史。

''',
    'docs/PROJECT_GOALS.md':'''**2026-11-16 v252当前起点**：同相态非相似剖面三链完成，独立30站复算和原始/截断积分账通过；34字段原面/剖面滑块及CSV已接。仍只1306球帽面/.1790893%候选，无新增合格面积；常H/同熵/β0—1范围不授真实三维横流。原MATLAB父监督终态缺失保留，单站软件修复回归32.047秒/六量差0/5拒绝通过，不补跑三链。按[v252执行点](validation/v252/NEXT_DECISION.md)直接鼻部之后真实外缘、熵层与横向流管入口，不重复鼻部或整箭参数支线；全面积/UI最高，总目标active。

'''}
    private={}
    for name,addition in additions.items():
        p=ROOT/name;live=p.read_text(encoding='utf8');head=subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT).decode('utf8')
        title,rest=head.split('\n\n',1);assert live.startswith(title+'\n\n')
        pos=live.index(rest[:100]);prefix=live[len(title)+2:pos];assert live.replace(prefix,'',1)==head
        private[name]={'prefix':prefix}
        own=title+'\n\n'+addition+rest
        if name=='docs/PROJECT_GOALS.md':
            start=own.index('当前下一项为**H2');end=own.index('\n\n',start)
            own=own[:start]+'''当前下一项为**鼻部之后3128目标面的真实外缘—熵层—横向流管入口**。v252已实现常H/同熵范围的实际非相似层流剖面、独立沿程复算、积分账和原面科研出口；不能再将“缺少实际剖面”作为同一区域重复开发的理由。按[v252执行点](validation/v252/NEXT_DECISION.md)复用CR-112180/v96/v241/v243/v244，一次形成成熟闭合采用决定并落地生产/区域接口与有区分力验证；不以球帽近似覆盖全后体。研究源资格逐级传播，强三维收缩/后体/底面仍需真实局部黏性反馈。原MATLAB父监督终态缺口与全部历史失败保留，不补跑旧三链、不追加ILU或整箭参数支线；全面积/UI独立最高，不等待32GB主机。'''+own[end:]
            marker='| 相态鼻部热剖面 |'
            start=own.index(marker);end=own.index('\n',start)
            row='\n| 同相态非相似剖面 | v252沿程后退差分/耦合bvp4c、三链120新站；独立30站/积分账/34字段与剖面CSV通过，单站修复回归通过 | 仍限原1306球帽面；无真实横流/变熵层/实验/CFD入口资格；原三链父监督终态缺失，完整结果保留 |'
            own=own[:end]+row+own[end:]
            start=own.index('每次先读AGENTS.md与本文件第4节H1—H5');end=own.index('\n\n',start)
            own=own[:start]+'''每次先读AGENTS.md与本文件第4节H1—H5，以[v252执行点](validation/v252/NEXT_DECISION.md)覆盖历史下一步。实际非相似剖面已完成受限验证，下一直接鼻部之后真实外缘/熵层/横向流管入口；先一次具名成熟方法决定，再执行生产/接口与具区分力验证。公开匹配基准可以承担诊断，真实箭体保持阶段集成和全面积责任。不重跑旧三链/物性/平板、不盲续整箭/扫参数，不把研究源升格accepted-edge；新任务另预算，所有原面/近壁/科学门槛保持。工程/全面积/UI独立最高，累计后置。'''+own[end:]
            # H2/H3正文只更新当前执行入口；前部历史记录不追改。
            start=own.index('### H2');end=own.index('### H4',start)
            region=own[start:end].replace('[v251执行点](validation/v251/NEXT_DECISION.md)','[v252执行点](validation/v252/NEXT_DECISION.md)').replace('[v251决定](validation/v251/NEXT_DECISION.md)','[v252决定](validation/v252/NEXT_DECISION.md)')
            region=region.replace('以非相似工程剖面生产补齐H2，并与H1/H3闭合真实外缘和反馈','复用v252非相似剖面，直接与H1/H3闭合真实外缘、熵层、横向流管和反馈')
            own=own[:start]+region+own[end:]
        dest=BASE/'owned_source'/name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(own,encoding='utf8')
        p.write_text(own.replace(title+'\n\n',title+'\n\n'+prefix,1),encoding='utf8')
    (BASE/'concurrent_prefixes.json').write_text(json.dumps(private,ensure_ascii=False,indent=2),encoding='utf8')
    p=ROOT/'README.md';s=p.read_text(encoding='utf8');marker='## 最新状态与交付（2026-11-16）\n\n'
    addition='''最新[v252非相似边界层剖面](docs/PHASE_NONSIMILAR_PROFILES_V252.md)：同相态Levy–Lees/沿程后退差分＋MATLAB bvp4c，三链及独立复算/积分账通过；上游壁温历史在40°处造成2.70%热流差。**仍只原1306球帽面候选，无新增合格热面积。** 原三链父监督终态缺失如实保留，单站修复回归已通过。34字段接原面WebGL，新增沿程滑块、实际速度/温度剖面与CSV；本地`outputs/v252_validation/workbench02/index.html`，图见[原面热流](docs/validation/v252/evidence/workbench02/rocket_nonsimilar.png)和[沿程/剖面](docs/validation/v252/evidence/workbench02/profiles_and_history.png)。

[下一执行点](docs/validation/v252/NEXT_DECISION.md)：鼻部之后的真实外缘—熵层—横向流管入口，不重做鼻部链或整箭参数支线；全面积/UI最高，总目标active。

'''
    assert marker in s;p.write_text(s.replace(marker,marker+addition,1).replace('最新[v251','前版[v251',1),encoding='utf8')

if __name__=='__main__':main()
