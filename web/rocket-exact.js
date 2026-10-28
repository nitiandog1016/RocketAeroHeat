/* 模块用途：三组真实原面通量研究，保留方法切换、新计数和逐变量未收敛身份。 */
'use strict';
(function(){
 const d=window.ROCKET_EXACT_STUDY,$=id=>document.getElementById(id),buttons=[...document.querySelectorAll('[data-step]')];
 const num=x=>Number.isFinite(x)?Number(x.toPrecision(6)).toString():'未知';let selected=d.cases.length-1,step=4;
 d.cases.forEach((c,i)=>{const o=document.createElement('option');o.value=String(i);o.textContent=c.label;$('caseChoice').append(o);});$('caseChoice').value=String(selected);
 function panel(){
  const c=d.cases[selected],s=c.summary,exact=c.interiorEulerFlux==='exact-godunov';
  const scheme=exact?'精确理想气体Godunov，xi=0波系采样':'Davis波界HLLC';
  const panels=[
   ['输入模型 · 原真实箭体','保留同一原19752面及20°完整三维来流。',[['三角面',d.sourceFaces],['坐标','SI / X鼻尖正向 / 零攻角来流-X']]],
   ['生成网格 · 原三维层域','真实原三角面层域；法向线仅作右预条件，不替代完整五量残差。',[['控制体 / 层数',d.cellCount+' / '+d.layerCount],['当前内部对流',scheme],['验证范围','几何与数值集成；网格/外域及物理精度未合格']]],
   ['添加边界条件 · 给定量','M5、30km、α20°；相同三维特征储库，静止无滑移墙。',[['壁温','300K输入，非求解Tw'],['气体','常γ='+d.gas.gamma+' / 冻结μ / Pr='+d.gas.Pr],['内面黏性功','算术平均速度代理；原应力与Fourier保留']]],
   ['计算 · '+c.label,'不同内部对流方程单独计数；旧HLLC六次仅为初猜来源。',[['内部对流',scheme],['本方法修正 / NS调用',s.acceptedNonlinearIterations+' / '+s.residualEvaluations],['本方程最大残差',num(s.maximumResidual)+' s⁻¹'],['RMS残差',num(s.rmsMerit)+' s⁻¹'],['稳态','未达到原2e−6 / 边界5e−10门槛']]],
   ['查看结果 · '+c.label,'原面云图、探针和任意截线共用当前状态与色标；每个场量保留方法。',[['邻壁气体温度',c.adjacentTemperatureRangeK.map(num).join(' ～ ')+' K'],['有符号入壁热功率',num(c.signedHeatPowerW/1000)+' kW（研究值）'],['身份','未收敛；Te/Tr/求解Tw未知'],['新增合格面积','0 m²']]]
  ];
  const p=panels[step];$('title').textContent=p[0];$('description').textContent=p[1];$('details').replaceChildren();
  p[2].forEach(([key,value])=>{const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=key;dd.textContent=String(value);$('details').append(dt,dd);});
  buttons.forEach((b,i)=>b.setAttribute('aria-selected',String(i===step)));
  const url=step===4?c.resultsURL:step===1?d.domainURL:d.geometryURL;if($('viewer').getAttribute('src')!==url)$('viewer').setAttribute('src',url);
  $('historyRows').replaceChildren();const h=c.history;
  if(h&&Array.isArray(h.iteration))h.iteration.forEach((it,j)=>{const tr=document.createElement('tr');[it,num(h.pseudoDt[j]),num(h.alpha[j]),num(h.maximumResidualPerSecond[j]),num(h.rmsMerit[j]),num(h.spatialConservationRelative[j]),h.residualEvaluations[j]].forEach(x=>{const td=document.createElement('td');td.textContent=String(x);tr.append(td);});$('historyRows').append(tr);});
  $('fluxNote').textContent='同Q基线对照，收缩内面70995：旧近似星压约0.003960Pa，精确Riemann星压 '+num(d.integrationReport.actualFacePressureStarPa)+' Pa，xi=0采样温度 '+num(d.integrationReport.actualFaceSampleTemperatureK)+' K。星压与界面采样不是壁压或热流解。';
  $('stateNote').textContent='当前 '+c.label+'；气体温度变化率来自当前空间算子，不是壁温历程。两套方程的残差历史分别记录；研究有限值不计入正式合格面积。';
 }
 buttons.forEach((b,i)=>b.addEventListener('click',()=>{step=i;panel();}));$('caseChoice').addEventListener('change',()=>{selected=Number($('caseChoice').value);panel();});panel();
})();
