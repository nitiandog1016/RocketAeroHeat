/* 模块用途：同窗五步骤、三组名义模型、原图符号及条件压力候选的科研对照。 */
'use strict';
(function(){
 const d=window.NASA_BOATTAIL_STUDY,$=id=>document.getElementById(id),buttons=[...document.querySelectorAll('[data-step]')];
 const number=x=>Number.isFinite(x)?Number(x.toPrecision(6)).toString():'未知';let selected=0,step=4;
 d.cases.forEach((c,i)=>{const o=document.createElement('option');o.value=String(i);o.textContent=c.label;$('caseChoice').append(o);});
 const colors=['#246b94','#a06028','#367657'];
 function panel(){
  const c=d.cases[selected],g=c.geometry;
  const panels=[
   ['输入模型 · 名义几何','按原图尺寸生成15°锥鼻、柱体及直尾锥；未生成支杆和砂粒带。',[['鼻底至尾底',number(g.noseBaseToBaseM)+' m（14 in）'],['柱径',number(g.diameterM)+' m'],['尾锥半角',c.angleDeg+'°'],['几何身份','名义闭合模型；非实验安装全几何']]],
   ['生成网格 · 几何检查','闭合、有向、连通的真实三角面；周向细化与解析体积/面积另有测试。',[['三角面',c.faceCount],['周向节点',g.resolution],['收缩轴长',number(g.boattailLengthM)+' m'],['本页网格用途','几何与关联映射；未建立NS/RANS体网格']]],
   ['添加边界条件 · 已知与未知','只按原报告限定零攻角和鼻型；绝对当地柱压与壁温未知。',[['攻角','0°'],['名义远场Mach / 当地柱面Mach','5.98 / 5.80，分别保留来源身份'],['压力分母','同工况当地柱面静压；不替换成远场压力'],['实验壁温 / 完整体入口','未知 / 未重建']]],
   ['计算 · 条件压力关联','实际复算声速流线膨胀压力比，再对照原图5个圆圈；不是完整NS/RANS求解。',[['算法','定比热空气γ=1.4的声速流线Prandtl-Meyer压力近似'],['声速Mach身份','边界层内流线近似，非外缘Mach'],['候选 p_bt/p_c',number(c.comparison.predictionRatio)],['实验比较身份','归一化压力关联对照；热流未验证']]],
   ['查看结果 · 候选压力比','同一工作台查看尾锥候选压力比、面元方法和任意截面曲线；其余区域与热场保持未知。',[['初始变量','候选压力比 p_bt/p_c / 1'],['共同色标',d.commonRatioRange.join(' ～ ')],['原图符号',String(c.comparison.points.length)+' 个'],['合格热流面积新增','0 m²']]]
  ];
  const p=panels[step];$('title').textContent=p[0];$('description').textContent=p[1];$('details').replaceChildren();
  p[2].forEach(([key,value])=>{const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=key;dd.textContent=String(value);$('details').append(dt,dd);});
  buttons.forEach((b,i)=>b.setAttribute('aria-selected',String(i===step)));
  const url=step===4?c.resultsURL:c.geometryURL;if($('viewer').getAttribute('src')!==url)$('viewer').setAttribute('src',url);
 }
 function draw(){
  const svg=$('pressureCurve');svg.replaceChildren();const ns='http://www.w3.org/2000/svg';
  const add=(tag,a={},text)=>{const e=document.createElementNS(ns,tag);Object.entries(a).forEach(([k,v])=>e.setAttribute(k,String(v)));if(text!==undefined)e.textContent=text;svg.append(e);return e;};
  const X=x=>65+x/.36*850,Y=y=>250-(y-.30)/.45*210;
  for(const x of [0,.06,.12,.18,.24,.30,.36]){add('line',{x1:X(x),x2:X(x),y1:40,y2:250,stroke:'#e0e6eb'});add('text',{x:X(x),y:269,'text-anchor':'middle'},number(x));}
  for(const y of [.3,.4,.5,.6,.7,.75]){add('line',{x1:65,x2:915,y1:Y(y),y2:Y(y),stroke:'#e0e6eb'});add('text',{x:55,y:Y(y)+4,'text-anchor':'end'},number(y));}
  add('text',{x:490,y:292,'text-anchor':'middle'},"原图 x′/D / 1（收缩肩部向下游）");add('text',{x:10,y:20},'p_bt / p_c / 1');
  d.cases.forEach((c,k)=>{
   if(!$('allCases').checked&&k!==selected)return;const points=c.comparison.points,ratio=c.comparison.predictionRatio,color=colors[k];
   add('line',{x1:65,x2:X(c.geometry.boattailLengthM/c.geometry.diameterM),y1:Y(ratio),y2:Y(ratio),stroke:color,'stroke-width':1.5,'data-role':'candidate'});
   points.forEach(p=>{const x=p.sourceFigureXprimeOverD,y=p.pressureRatio,h=p.ratioDigitizationHalfWidth,w=p.xDigitizationHalfWidth;
    add('line',{x1:X(x),x2:X(x),y1:Y(y-h),y2:Y(y+h),stroke:color,'data-role':'read-interval'});
    add('line',{x1:X(x-w),x2:X(x+w),y1:Y(y),y2:Y(y),stroke:color});
    const e=add('circle',{cx:X(x),cy:Y(y),r:4,fill:'#fff',stroke:color,'stroke-width':1.5,'data-id':p.id});
    e.addEventListener('click',()=>{$('pointDetails').textContent=p.id+'；原像素 ('+p.pixelXY.join(', ')+')；原图压力比 '+number(y)+' ± '+number(h)+'（仅读图）；未强行匹配测孔；绝对压力/Tw未知。';});
   });add('text',{x:670+k*78,y:22,fill:color},c.angleDeg+'°');
  });
 }
 function table(){
  $('comparisonRows').replaceChildren();d.cases.forEach((c,i)=>{const tr=document.createElement('tr');tr.setAttribute('aria-selected',String(i===selected));
   [c.label,c.comparison.points.length,number(c.comparison.predictionRatio),number(100*c.comparison.meanAbsoluteRelativeDifference)+'%',number(100*c.comparison.maximumAbsoluteRelativeDifference)+'%'].forEach(v=>{const td=document.createElement('td');td.textContent=String(v);tr.append(td);});$('comparisonRows').append(tr);});
 }
 buttons.forEach((b,i)=>b.addEventListener('click',()=>{step=i;panel();}));
 $('caseChoice').addEventListener('change',()=>{selected=Number($('caseChoice').value);panel();table();draw();});
 $('allCases').addEventListener('change',draw);panel();table();draw();
})();
