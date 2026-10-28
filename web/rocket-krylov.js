/* 模块用途：同原面五组功策略/实际修正对照，保留非线性/物理时间和算法身份。 */
'use strict';
(function(){
 const d=window.ROCKET_KRYLOV_STUDY,$=id=>document.getElementById(id),buttons=[...document.querySelectorAll('[data-step]')];
 const num=x=>Number.isFinite(x)?Number(x.toPrecision(6)).toString():'未知';let selected=d.cases.length-1,step=4;
 d.cases.forEach((c,i)=>{const o=document.createElement('option');o.value=String(i);o.textContent=c.label;$('caseChoice').append(o);});$('caseChoice').value=String(selected);
 function panel(){
  const c=d.cases[selected],s=c.summary;
  const panels=[
   ['输入模型 · 原真实箭体','同一原19752面，保留X鼻尖正向及20°完整矢量来流。',[['三角面',d.sourceFaces],['场量身份','原面常值；没有周向复制零攻角结果']]],
   ['生成网格 · 三维法向线','原217272CV/11层；每原面55未知量的5×5块法向线仅作求解器预条件。',[['控制体 / 层数',d.cellCount+' / '+d.layerCount],['原空间方程','共享面HLLC + WLS应力/导热；未用线近似替代'],['网格资格','几何已验，网格/外域及物理精度尚未合格']]],
   ['添加边界条件 · 给定量','M5、30km、α20°、三分量特征储库；Tw给定300K，无滑移静止墙。',[['内面功速度',c.arithmeticWork?'显式一阶算术平均速度代理':'旧WLS重构速度'],['气体','常γ='+d.gas.gamma+'，冻结μ，Pr='+d.gas.Pr],['壁温身份','300K输入；外缘/恢复/求解壁温未知']]],
   ['计算 · '+c.label,'原空间五量残差重新验收；伪步和Krylov迭代均不算飞行或BE时间。',[['本段接受 / 原NS调用',s.acceptedNonlinearIterations+' / '+s.residualEvaluations],['累计接受 / 原NS调用',(s.totalAcceptedNonlinearIterations??s.acceptedNonlinearIterations)+' / '+(s.totalResidualEvaluations??s.residualEvaluations)],['原最大残差',num(s.maximumResidual)+' s⁻¹'],['RMS下降身份','仅线搜索规则；不替代最大残差稳态门槛'],['稳态','未达到原2e−6 / 定常边界5e−10要求']]],
   ['查看结果 · '+c.label,'比较原面墙压、热流、邻壁气体温度和算子温度变化率；共同色标、探针、任意截线及CSV。',[['邻壁气体温度',c.adjacentTemperatureRangeK.map(num).join(' ～ ')+' K'],['总有符号入壁热功率',num(c.signedHeatPowerW/1000)+' kW（研究场）'],['时间 / 资格','无物理时间推进；未收敛/未取得物理资格'],['新增合格面积','0 m²']]]
  ];
  const p=panels[step];$('title').textContent=p[0];$('description').textContent=p[1];$('details').replaceChildren();
  p[2].forEach(([key,value])=>{const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=key;dd.textContent=String(value);$('details').append(dt,dd);});
  buttons.forEach((b,i)=>b.setAttribute('aria-selected',String(i===step)));
  const url=step===4?c.resultsURL:step===1?d.domainURL:d.geometryURL;if($('viewer').getAttribute('src')!==url)$('viewer').setAttribute('src',url);
  $('historyRows').replaceChildren();const h=c.history;
  if(h&&Array.isArray(h.iteration))h.iteration.forEach((it,j)=>{const tr=document.createElement('tr');[it,num(h.pseudoDt[j]),num(h.alpha[j]),num(h.maximumResidualPerSecond[j]),num(h.rmsMerit[j]),num(h.spatialConservationRelative[j]),h.residualEvaluations[j]].forEach(x=>{const td=document.createElement('td');td.textContent=String(x);tr.append(td);});$('historyRows').append(tr);});
  $('workNote').textContent='同态受控检查：旧内面重构功速度最大 '+num(d.workReport.maximumOldWorkSpeed)+' m/s，算术代理最大 '+num(d.workReport.maximumNewWorkSpeed)+' m/s；只改变内部能量通量，壁热/牵引等保持。负配对功仍有 '+d.workReport.negativePairCount+' 面，不能宣称一般正性。';
  $('stateNote').textContent='当前 '+c.label+'；温度变化率来自空间算子，不是壁温历程。RMS下降与原最大残差分别列示；各离散方法的迭代成本不能当作同一物理误差的加速比。';
 }
 buttons.forEach((b,i)=>b.addEventListener('click',()=>{step=i;panel();}));$('caseChoice').addEventListener('change',()=>{selected=Number($('caseChoice').value);panel();});panel();
})();
