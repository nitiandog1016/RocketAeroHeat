/* 模块用途：已确认单参数工况在同一五步骤宿主切换；分布比较保留同一真实墙段和单位。 */
(function(){'use strict';
 const d=window.ROCKET_SENSITIVITY,$=id=>document.getElementById(id),buttons=[...document.querySelectorAll('[data-step]')];
 if(!d||!['complete-rocket-single-parameter-sensitivity-workbench-v1','complete-rocket-wall-resolution-workbench-v1','complete-rocket-contact-discretization-workbench-v1'].includes(d.kind)||d.gridIndependentClaimed!==false)throw Error('敏感性快照身份无效');
 const colors=['#34495e','#2879a5','#b17035','#52875c'],fmt=x=>Number(x).toPrecision(5),ns='http://www.w3.org/2000/svg';let chosen=0,stage=4;
 d.cases.forEach((c,i)=>{const o=document.createElement('option');o.value=String(i);o.textContent=c.label;$('caseChoice').append(o);});
 const titles=['输入模型','生成网格','添加边界条件','计算','查看结果'];
 function rowsTable(target,head,rows){const table=document.createElement('table'),tr=document.createElement('tr');head.forEach(x=>{const th=document.createElement('th');th.textContent=x;tr.append(th);});table.append(tr);rows.forEach(row=>{const tr=document.createElement('tr');row.forEach(x=>{const td=document.createElement('td');td.textContent=String(x);tr.append(td);});table.append(tr);});target.replaceChildren(table);}
 function select(){const c=d.cases[chosen],s=c.study,m=c.metrics;const flux=s.numericalStrategies?.interiorConvectiveFluxPolicy||'unknown';const fluxLabel=flux==='hllc'?'HLLC':flux==='hlle-einfeldt'?'HLLE-Einfeldt':flux;
   $('title').textContent=titles[stage]+' · '+c.label;
   $('description').textContent=['四档使用同一个完整闭合箭体，包含真实收缩、后体和底面。',d.meshDescription||'分别改变径向层数、近壁偏移或外域半径；未进行流向加密。','同一30 km / Mach 5 / 零攻角；远场、静止墙、Tw和SA输入保持。','本工况已按原门槛达到数值定常；Newton准备与BE确认分别列出。','三维云图保留面元探针、各变量方法与任意X截面周向曲线。压力、热流和气体温度采用四档统一色标。'][stage];
   if(d.gridChanges&&(stage===3||stage===4))$('description').textContent+=' HLLE两网格总热变化 '+fmt(100*d.gridChanges.hlle.signedHeatChangeFraction)+'%；HLLC '+fmt(100*d.gridChanges.hllc.signedHeatChangeFraction)+'%，HLLC面积加权热流绝对差 '+fmt(100*d.gridChanges.hllc.areaWeightedHeatAbsoluteChangeFraction)+'%。这不是实验验证或网格独立结论。';
   const details=[['源模型',s.sourceFaces+' 三角面；'+fmt(s.totalLengthM)+' m'],['流场网格',m.cellCount+' 格 / '+m.radialLayers+' 层'],['首节点法向偏移',fmt(m.firstNormalOffsetM*1e6)+' μm'],['实际体质心墙距',s.nearWallNormalDistanceM.map(x=>fmt(x*1e6)).join(' ～ ')+' μm'],['外域半径',fmt(m.outerRadiusM)+' m'],['原定常残差',fmt(m.originalResidual)+' s⁻¹；门槛2e-6'],['Newton / BE',m.NewtonIterations+' 次准备 / 3个独立BE端点'],['总有符号热功率',fmt(m.totalSignedHeatW/1000)+' kW'],['热流算法','SA-noft2 RANS / Fourier / Sutherland / 给定等温墙'],['压力算法','内部'+fluxLabel+'通量 / 五量RANS；墙owner压力代理'],['温度身份','邻壁单元气体代理；Tw给定，Te/Tr/求解Tw未知']];
   $('details').replaceChildren();details.forEach(([name,value])=>{const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=name;dd.textContent=value;$('details').append(dt,dd);});
   buttons.forEach((b,i)=>b.setAttribute('aria-selected',String(i===stage)));$('iterations').hidden=stage!==3;
   const url=stage===4?c.resultsURL:c.geometryURL;if($('viewer').getAttribute('src')!==url)$('viewer').setAttribute('src',url);
   rowsTable($('newtonRows'),['Newton迭代','原残差 / s⁻¹','边界/源相对账'],s.rootIterationHistory.map(r=>[r.iteration,fmt(r.residual),fmt(r.boundaryBalance)]));
   rowsTable($('beRows'),['BE步','原残差 / s⁻¹','实际伪步 / s','五量相对账'],s.steps.map(r=>[r.number,fmt(r.residual),fmt(r.dt),fmt(r.balance)]));
   [...$('comparisonRows').children].forEach((r,i)=>r.setAttribute('aria-selected',String(i===chosen)));draw();drawProfile();
 }
 function node(tag,attrs,text){const e=document.createElementNS(ns,tag);Object.entries(attrs).forEach(([k,v])=>e.setAttribute(k,String(v)));if(text!==undefined)e.textContent=String(text);return e;}
 function draw(){const svg=$('comparisonCurve');svg.replaceChildren();const key=$('curveField').value,scale=['signedHeatWm2','pressurePa'].includes(key)?.001:1;
   const xy=d.cases[0].metrics.fields.wallCentroidXR;let arc=[0];for(let i=1;i<xy.length;i++)arc.push(arc[i-1]+Math.hypot(xy[i][0]-xy[i-1][0],xy[i][1]-xy[i-1][1]));
   // 墙心到墙心弧长：首墙心定义为0；不把底面同X的点重合成一条流向曲线。
   const values=d.cases.flatMap(c=>c.metrics.fields[key].map(v=>v*scale));let lo=Math.min(...values),hi=Math.max(...values);if(lo===hi){lo-=1;hi+=1;}
   const pad=.06*(hi-lo);lo-=pad;hi+=pad;const x=v=>65+v/arc.at(-1)*860,y=v=>210-(v-lo)/(hi-lo)*160;
   for(let i=0;i<5;i++){const a=lo+(hi-lo)*i/4;svg.append(node('line',{x1:65,y1:y(a),x2:925,y2:y(a),stroke:'#e3e8ed'}));svg.append(node('text',{x:58,y:y(a)+4,'text-anchor':'end'},fmt(a)));const b=arc.at(-1)*i/4;svg.append(node('text',{x:x(b),y:232,'text-anchor':'middle'},fmt(b)));}
   d.cases.forEach((c,i)=>{const points=c.metrics.fields[key].map((v,k)=>x(arc[k])+','+y(v*scale)).join(' ');svg.append(node('polyline',{points,fill:'none',stroke:colors[i],'stroke-width':i===chosen?2.4:1.3}));svg.append(node('text',{x:75+215*i,y:22,fill:colors[i]},c.label));});
   svg.append(node('text',{x:495,y:252,'text-anchor':'middle'},'相邻墙段中心累积距离 / m（首墙心为0）'));
 }
 function drawProfile(){
   const available=d.cases.every(c=>Array.isArray(c.wallNormalProfiles));if(!available)return;
   const svg=$('profileCurve');svg.replaceChildren();const index=Number($('profileWall').value)-1,key=$('profileField').value;
   if(!Number.isInteger(index)||index<0||index>=d.cases[0].wallNormalProfiles.length)throw Error('未知墙段');
   const first=d.cases[0].wallNormalProfiles[index],limit=Number($('profileExtent').value)*.001;
   $('profilePosition').textContent='第'+(index+1)+'段 · '+first.region+' · X='+fmt(first.wallCentroidXR[0])+' m';
   $('profileBoundaryNote').textContent='墙边界：Tw='+first.wallTemperatureInputK+' K，速度0；均为给定，未混入内部计算点。';
   const face=['EulerEnergyFluxMinusArithmeticPhysicalWm2','viscousTotalEnergyFluxWm2'].includes(key);
   const fieldMethod=key==='temperatureK'?'由同一实际Q按定比热空气EOS恢复内部气体温度；':key==='muTOverMu'?'同一实际SA-noft2有效黏度除以Sutherland层流黏度；':face?'保存的内面实际能量通量；Euler项另减两侧算术物理参考，仅为离散诊断；':'由同一实际Q和真实墙切向/位置读取；';
   $('profileAlgorithmText').textContent=fieldMethod+(d.cases[chosen].profileSourceMethod||'实际端点未提供方法说明');
   if(!(Number.isFinite(limit)&&limit>0)){svg.append(node('text',{x:65,y:70},'请输入有限正的最大距离。'));return;}
   const scale=['pressurePa','EulerEnergyFluxMinusArithmeticPhysicalWm2','viscousTotalEnergyFluxWm2'].includes(key)?.001:1;
   const samples=d.cases.map(c=>{const p=c.wallNormalProfiles[index],dist=face?p.radialFaceNormalDistanceM:p.distanceM;return dist.map((x,i)=>[x,p[key][i]*scale]).filter(([x,y])=>x>0&&x<=limit&&Number.isFinite(y));});
   const all=samples.flat();if(!all.length){svg.append(node('text',{x:65,y:70},'当前范围没有实际控制体/内面采样；未补零。'));return;}
   const log=$('profileLog').checked,minX=log?Math.min(...all.map(p=>p[0])):0,maxX=limit;
   if(minX>=maxX){svg.append(node('text',{x:65,y:70},'请扩大法向距离范围。'));return;}
   const transform=x=>log?Math.log10(x):x,xlo=transform(minX),xhi=transform(maxX);
   const x=v=>70+(transform(v)-xlo)/(xhi-xlo)*855;let lo=Math.min(...all.map(p=>p[1])),hi=Math.max(...all.map(p=>p[1]));if(lo===hi){lo-=1;hi+=1;}const pad=.06*(hi-lo);lo-=pad;hi+=pad;const y=v=>225-(v-lo)/(hi-lo)*170;
   for(let i=0;i<5;i++){
    const value=lo+(hi-lo)*i/4;svg.append(node('line',{x1:70,y1:y(value),x2:925,y2:y(value),stroke:'#e3e8ed'}));svg.append(node('text',{x:62,y:y(value)+4,'text-anchor':'end'},fmt(value)));
    const coordinate=xlo+(xhi-xlo)*i/4,valueX=log?10**coordinate:coordinate;
    svg.append(node('text',{x:70+855*i/4,y:247,'text-anchor':'middle'},fmt(1e3*valueX)));
   }
   samples.forEach((points,i)=>{if(points.length)svg.append(node('polyline',{points:points.map(([a,b])=>x(a)+','+y(b)).join(' '),fill:'none',stroke:colors[i],'stroke-width':i===chosen?2.4:1.3}));svg.append(node('text',{x:75+215*i,y:22,fill:colors[i]},d.cases[i].label));});
   svg.append(node('text',{x:495,y:270,'text-anchor':'middle'},'实际墙法向投影距离 / mm'+(log?'（对数）':'')));
 }
 if(d.cases.every(c=>Array.isArray(c.wallNormalProfiles))){$('wallProfilePanel').hidden=false;$('profileWall').max=d.cases[0].wallNormalProfiles.length;['profileWall','profileField','profileExtent','profileLog'].forEach(id=>$(id).addEventListener(id==='profileWall'?'input':'change',drawProfile));}
 d.cases.forEach((c,i)=>{const tr=document.createElement('tr');[c.label,c.metrics.cellCount,fmt(c.metrics.originalResidual),fmt(c.metrics.totalSignedHeatW/1000),i===0?'基线':fmt(100*(c.metrics.totalSignedHeatW/d.cases[0].metrics.totalSignedHeatW-1))+'%',fmt(c.metrics.yPlusProxyRange[1])].forEach(v=>{const td=document.createElement('td');td.textContent=String(v);tr.append(td);});$('comparisonRows').append(tr);});
 if(Array.isArray(d.failedTrials)&&d.failedTrials.length){const panel=document.createElement('details'),summary=document.createElement('summary');summary.textContent='未收敛试验（不计作有效定常结果）';panel.append(summary);d.failedTrials.forEach(t=>{const p=document.createElement('p');p.textContent=t.kind+'：'+t.stopReason+'，原残差 '+fmt(t.originalResidual)+' s⁻¹；最低流体温度 '+fmt(t.minimumFluidTemperatureK)+' K；保护拒绝 '+JSON.stringify(t.guardRejections)+'。原失败记录保留；当前有效精网格使用粗HLLC根初猜。';panel.append(p);});$('iterations').append(panel);}
 $('caseChoice').addEventListener('change',()=>{chosen=Number($('caseChoice').value);if(!Number.isInteger(chosen)||chosen<0||chosen>=d.cases.length)throw Error('未知工况');select();});
 $('curveField').addEventListener('change',draw);buttons.forEach((b,i)=>b.addEventListener('click',()=>{stage=i;select();}));select();
})();
