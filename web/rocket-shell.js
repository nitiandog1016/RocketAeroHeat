/* 模块用途：真实20度原面数值研究的同窗五步骤与三组初态/启动状态。 */
'use strict';
(function(){
 const d=window.ROCKET_SHELL_STUDY,$=id=>document.getElementById(id),buttons=[...document.querySelectorAll('[data-step]')];
 const number=x=>Number.isFinite(x)?Number(x.toPrecision(6)).toString():'未知';let selected=2,step=4;
 d.cases.forEach((c,i)=>{const o=document.createElement('option');o.value=String(i);o.textContent=c.label;$('caseChoice').append(o);});$('caseChoice').value=String(selected);
 function panel(){
  const c=d.cases[selected],g=d.gas,r=d.advanceReport;
  const panels=[
   ['输入模型 · 原真实箭体','原19752三角面直接建立独立三维计算体域，未改动既有表面重网格。',[['原三角面',d.sourceFaces],['坐标','X鼻尖正向；SI / m'],['模型身份','真实源模型；不是周向复制的零攻角解']]],
   ['生成网格 · 真实三维层域','检查原墙、近壁层和外边界真实节点。外边界不是激波。',[['控制体 / 层数',d.cellCount+' / '+d.layerCount],['第一层真实壁距',d.wallNormalDistanceRangeM.map(x=>number(x*1e6)).join(' ～ ')+' μm'],['体网格资格','正体积与几何账已校核；尚未做网格独立性']]],
   ['添加边界条件 · 完整矢量来流','M5、30 km、攻角20°；完整三分量远场特征储库，给定等温无滑移墙。',[['给定壁温','300 K；不是求解壁温'],['模型','γ='+number(g.gamma)+'，μ='+number(g.mu)+' Pa·s，Pr='+number(g.Pr)],['适用身份','定物性层流数值研究，无湍流/真实空气资格']]],
   ['计算 · 有界三步启动','真实共享dt SSP-RK2；原正性、温域及守恒护栏保留。步数上限不代表收敛。',[['接受 / 拒绝重试',r.acceptedSteps+' / '+r.rejectedAttempts],['启动积分时间',number(r.startupIntegrationTimeSeconds)+' s'],['t / (L/U∞)',number(r.fractionOfReferenceConvectiveTime)],['收敛状态','未收敛；不能用于正式飞行气动热判断']]],
   ['查看结果 · '+c.label,'同一原三角面比较墙压、入壁热流、邻壁气体温度及几何来流投影；支持探针、任意X截线和CSV。',[['有符号总入壁热功率',number(c.signedHeatPowerW/1000)+' kW（研究场）'],['当前数值积分时间',number(c.startupIntegrationTimeSeconds)+' s，非飞行时间'],['公共色标','压力 / 热流 / 邻壁气体温度共用三组范围'],['生产合格热流面积新增','0 m²']]]
  ];
  const p=panels[step];$('title').textContent=p[0];$('description').textContent=p[1];$('details').replaceChildren();
  p[2].forEach(([key,value])=>{const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=key;dd.textContent=String(value);$('details').append(dt,dd);});
  buttons.forEach((b,i)=>b.setAttribute('aria-selected',String(i===step)));
  const url=step===4?c.resultsURL:step===1?d.domainURL:d.geometryURL;
  if($('viewer').getAttribute('src')!==url)$('viewer').setAttribute('src',url);
  $('heatNote').textContent=selected===0?'旧初态保留184个错误热向面作为失败对照，不能作为可用热环境。':selected===1?'修正初态全体Tw高于气体T，热量由壁向气体；符号未裁剪。':'三步后出现的正负热流是当前未收敛状态派生，尚未形成可信飞行热流分布。';
 }
 d.history.forEach(row=>{const tr=document.createElement('tr');[row.acceptedStep,number(row.startupIntegrationTimeSeconds),number(row.dt),row.retries,number(row.residualPerSecond),number(row.globalBalanceRelative)].forEach(x=>{const td=document.createElement('td');td.textContent=String(x);tr.append(td);});$('historyRows').append(tr);});
 $('costNote').textContent='三步实际耗时 '+number(d.advanceReport.elapsedSeconds)+' s，t/(L/U∞)='+number(d.advanceReport.fractionOfReferenceConvectiveTime)+'。显式微步成本过高；下一步评估稳态/半隐式区域求解，不盲目延长显式运行。';
 buttons.forEach((b,i)=>b.addEventListener('click',()=>{step=i;panel();}));$('caseChoice').addEventListener('change',()=>{selected=Number($('caseChoice').value);panel();});panel();
})();
