# -*- coding: utf-8 -*-
"""模块职责：归档后更新当前目标与关闭一次性入口；不执行数值任务。"""
from pathlib import Path
root=Path(__file__).resolve().parents[2]
p=root/'AGENTS.md';s=p.read_text(encoding='utf8');header='# RocketAeroHeat 项目执行约定\n';assert s.startswith(header)
new='''
**2026-11-16 v242当前起点：保持原壁的成熟PLC网格通过，真实前体400更新正常预算结束。** 13坏固定面均来自人工裁切（10远场/3内部），采用单壳零孔；原4434面/2241壁点/46环/0.3214887716893075m²与截断/外域保持。一次Gmsh4.15.2表面＋TetGen1.5.1生成20277tet/5354点，独立最低SICN0.133/无<0.1，Netgen未运行。原壁/导入/闭合/连通/文本往返通过；只编两个薄入口，成熟数值对象复用。见docs/FOREBODY_FIXED_WALL_PLC_V242.md。

同相态Euler16/Roe/原二阶/CFL和case.cfg逐字节保持，新均匀初态非旧200续算。唯一400更新29.422秒exit0/峰66039808B/1CPU，200/300/400收益规则均未停；不称收益停止。0/1/每50/400及强制原生末CSV逐位同；783出口点/313200调用/1566000五量分量逐位同，最低Mn2.246966。末预处理拒绝最大0不是全部内部事件证明。四式logRMS到−8，物理能量−7.84267未过；净质量/能量1.3023e−16/6.9565e−16及压力变化1.0282e−13通过，H均差0.8610%/最大4.7061%，外缘/新热资格0。

原物理九点入口/方向/标签/步长不动，四管共同136/135/137/56共464；5303站独立核对，36交点28壁/8出口，6496数CSV/科研图通过。H沿程漂移1.358—3.571%、局部非积分质量比0.832—2.093，仍不具外缘资格，碰壁非分离证明。10阶段63字段/266040原面值及前端纯函数通过，未知灰/q与固体Tw未知，图已看，无MATLAB/共享WebGL/浏览器。13监督任务全部终态稳定组空/无数值失败；离线15.500秒、审核出口14.156、两薄构建3.781/8.391、表面.938/体1.578、CFD29.422，归档另计；开发定向读取另述。文档批量命令首次stdin编码解析失败未执行，改UTF8文件后恢复，独立保留说明。

**下一次按docs/validation/v242/NEXT_DECISION.md直接具名成熟外流重构/真实流管入口采用与实施。** 复用新保存场，在原四管经过区域区分源场H/S空间误差、P1与原生双控制体通量/滑移壁不相容、壁邻→外缘/熵层及横向标签。一次≤120秒定向方法/源决定后做具体生产器或具名拒绝；必要未保存原生面通量/梯度另预算，不用净残差反推唯一通量、不重跑补日志。无成熟支持则具名局部黏性真实入口/反馈，不长期自研通用重构/通量Jacobian。不能仅因能量−7.84盲续400/200/800/199/92、放−8或称已证舍入极限；不重复旧网格/物性/小壳/冻结整域/扫参数，不移种子/投影/冻结H/S补解。当前固定CLI关闭、无新CFD/构建/查询准入。原面/全分母/近壁/所有科学门槛及3128目标、收缩5483/后体/底面责任保持，全面积/UI最高、累计后置；保护29/engineering54/候选/并发web，小证据GitHub大场E盘，无包/独有删除，不新代理/computer-use，总目标active，下方历史。
'''
p.write_text(header+new+s[len(header):],encoding='utf8')
p=root/'docs/PROJECT_GOALS.md';lines=p.read_text(encoding='utf8').splitlines()
for i,line in enumerate(lines):
    if line.startswith('**2026-11-16 v241当前起点**'):
        lines[i]='**2026-11-16 v242当前起点**：保持原4434壁面的单壳零孔成熟网格实际通过，20277tet/5354点、独立最低SICN0.133；唯一真实前体400更新正常完成（29.422秒/约66MB），净质量/能量与压力稳定通过，物理能量RMS−7.84267仍未达−8。原四管同物理入口的H/S与横向流率偏差仍在，外缘/新热资格0。10阶段63字段/266040原面值已接工作台，464共同流管站及6496数CSV/对照图通过。按[v242执行点](validation/v242/NEXT_DECISION.md)直接处理成熟外流重构和真实工程入口，不盲续400或回到大内存整域前置。全面积/科研UI最高，总目标active。'
    elif line.startswith('当前下一项仍为**H1/H2/H3'):
        lines[i]='当前下一项为**H1/H2/H3区域外流—工程边界层—必要黏性反馈的真实入口闭环**。v242已解决人工坏固定面并获得显著下降的区域残差与净收支，但同物理起点的四管H偏差、局部质量流率偏差和壁面停止未解决。一次具名成熟外流重构/流管方案采用，核对源场空间误差、原生双控制体通量/滑移壁与P1插值相容、壁邻→外缘/熵层及横向尺度；必要原生未保存量另预算，不能反推唯一面通量或重跑补日志。没有成熟支持则转具名局部黏性真实入口/反馈，不长期自研通用内核。不得移种子、投影失败流线或冻结H/S补解；−8及空间/外域/热资格门槛保持。收缩旧5483面、后体/底面等全部责任保留，工程/UI不因单项CFD门槛停工。'
    elif line.startswith('v224同模型工程热链、v225波系'):
        lines[i]='v224同模型工程热链、v225波系、v226鼻部外缘、v227局部相似热剖面/v228修复、v240区域生产与v241九线入口均复用。v242新网格/400发展及原四管比较已实际完成，但尚无可信新外缘。下一按[v242执行点](validation/v242/NEXT_DECISION.md)直接成熟外流重构/真实工程入口的具名采用和实施，原3128面目标保持。研究p/T/H/S/U和流管可查看；壁邻或稳定H/S不是完整三维BL外缘/局部CFD入口。'
    elif line.startswith('v234具体多区组合保持关闭'):
        lines[i]='v234具体多区组合保持关闭，v235—v242的原生入口、实际发展、接受/终点证据及所有历史失败/中断保留。新v242400按预算结束不盲续。按[v242决定](validation/v242/NEXT_DECISION.md)与H2协同建立有真实入口和下游影响依据的区域合同；近壁P1四管不能直接作黏性入口。局部化不等于任意截断。新控制/预算前不运行CFD；原−8/y+/热稳定/收支/空间外域与公开验证保持。'
    elif line.startswith('每次先读AGENTS.md与本文件第4节'):
        lines[i]='每次先读AGENTS.md与本文件第4节H1—H5，以[v242执行点](validation/v242/NEXT_DECISION.md)覆盖历史下一步。新合格网格、400真实发展、准确出口/末CSV及原四管已完成，直接成熟外流重构/真实流管入口采用与实施，复用已保存源；不重复全表/冻结/物性/小壳或盲续400/旧200/800/199/92。必要未保存原生量/新区域另预算，当前固定CLI关闭、无新CFD/构建/查询准入。不放宽科学门槛、不以研究壁邻当accepted-edge。新原面/流管数据接科研工作台，全面积/UI最高且独立，累计后置。'
idx=next(i for i,line in enumerate(lines) if line.startswith('| 内部三维体流管 |'))+1
lines.insert(idx,'| 保持原壁的成熟PLC与区域发展 | v242一次Gmsh/TetGen，独立SICN≥0.133；真实400更新完整保存/末CSV，净收支/压力稳定及原面63字段通过；同种子5303流管站审核 | 能量RMS−7.84267尚未过−8，H/S/质量流管与壁邻→外缘/空间外域未过；物理资格和新热面积0，不盲续400 |')
p.write_text('\n'.join(lines)+'\n',encoding='utf8')
p=root/'README.md';s=p.read_text(encoding='utf8');anchor='## 最新状态与交付（2026-11-16）\n'
new='''
最新[v242保持原壁的前体重建](docs/FOREBODY_FIXED_WALL_PLC_V242.md)：成熟Gmsh/TetGen网格通过原质量门槛，20277tet/5354点；内部前体实际400更新正常完成，用时29.422秒/采样峰约66MB。净质量/能量与压力稳定通过，**物理能量RMS及真实外缘尚未合格，没有新增合格壁热面积**。无需等待更大内存主机才能推进混合路线。

本地`outputs/v242_validation/workbench01/index.html`或`forebody_research.json`提供10阶段63字段，266040个原面研究值核对通过。[科研图](docs/validation/v242/evidence/visual01/forebody_research.png)保留全箭灰色未知区和前体放大；[同入口流管对照](docs/validation/v242/evidence/tube_audit01/tube_comparison.png)表明H/S与局部流管误差仍需处理。无共享WebGL修改/浏览器验收。[当前下一执行点](docs/validation/v242/NEXT_DECISION.md)是成熟外流重构与真实工程入口，不盲续400。
'''
assert anchor in s;s=s.replace(anchor,anchor+new,1).replace('最新[v241内部三维流管]','旧[v241内部三维流管]',1);p.write_text(s,encoding='utf8')
p=root/'.gitattributes';s=p.read_text(encoding='utf8');s+='\n# v242 成熟PLC、实际区域发展和准确执行源保留原字节。\ndocs/validation/v242/evidence/** -text -eol -whitespace\n';p.write_text(s,encoding='utf8')
for name in ('prepareForebodyPLCV242','surfaceForebodyPLCV242','buildForebodyPLCV242','meshForebodyPLCV242','auditForebodyPLCV242'):
    p=root/'scripts'/(name+'.py');s=p.read_text(encoding='utf8');start=s.index("if __name__=='__main__':")
    p.write_text(s[:start]+"if __name__=='__main__':raise SystemExit('v242 fixed task closed; use NEXT_DECISION')\n",encoding='utf8')
for name in ('runForebodyPLCV242','runForebodyEulerV242','auditForebodyEulerV242','exportForebodyEulerV242','renderForebodyEulerV242','compareForebodyTubesV242'):
    p=root/'scripts'/(name+'.py');s=p.read_text(encoding='utf8');assert s.count('CLOSED=False')==1;p.write_text(s.replace('CLOSED=False','CLOSED=True'),encoding='utf8')
