/* 模块职责：稀疏原面结果接既有WebGL；46入口积分量、来源及CSV同页查看。 */
'use strict';
function applyCoupledSurface(model,data){
  if(model.faces.length!==data.originalFaces)throw new Error('ORIGINAL_FACE_COUNT');
  const n=data.originalFaces,known=new Set(data.faceIDs);
  for(const spec of data.surfaceFields){
    const values=Array(n).fill(null),status=Array(n).fill(0),reason=Array(n).fill('不属于本轮已求解球帽；保留未知');
    data.withheldFaceIDs.forEach(id=>{reason[id-1]='驻点起步区两档步长差超2%或整面跨入该区；原始值保留在诊断出口，不进入本轮默认统计';});
    data.faceIDs.forEach((id,j)=>{if(id<1||id>n||!Number.isFinite(spec.values[j]))throw new Error('INVALID_SURFACE_VALUE');values[id-1]=spec.values[j];status[id-1]=1;reason[id-1]='';});
    if(known.size!==data.faceIDs.length)throw new Error('DUPLICATE_FACE_ID');
    model.fields.push({...spec,values,status,missingReason:reason,methodIndex:Array(n).fill(0),methodNames:[data.method],
      computedStateLabel:'变熵积分工程候选',availability:'available',scope:data.scope,
      stateNote:'积分量已有值；实际剖面和三维/实验资格未知。给定Tw=300 K，求解固体Tw仍未知。'});
  }
  model.defaultField='entropyIntegralHeat';model.caseNote=data.caseNote;
  model.researchIdentity={...model.researchIdentity,coupledEntropy:data.identity};
  model.resultCondition={caption:'M5 · 30 km · α20° · 给定Tw300 K · 原19752面 · v254球帽积分工程候选'};
  model.comparisonColorRanges={...model.comparisonColorRanges,...data.colorRanges};
  return model;
}
function integralCSV(data,branch){
  const b=data.branches[branch];if(!b)throw new Error('UNKNOWN_BRANCH');
  const keys=['index','upstreamFaceID','downstreamFaceID','phiDeg','thetaDeg','pressurePa','heatFluxWm2','edgeTemperatureK','edgeSpeedMps','massThicknessM','betaPressureProxy','fittedMetricM','actualMetricM','givenWallTemperatureK','solvedSolidWallTemperatureK','profileReadyForCFD','physicalQualification'];
  return keys.concat(['sourceSHA256','method']).join(',')+'\n'+b.stations.map(s=>keys.map(k=>s[k]===null?'':s[k]).concat([b.sourceSHA256,data.method]).join(',')).join('\n')+'\n';
}
function integralSegments(x,y){const lines=[];let line=[];x.forEach((v,i)=>{if(Number.isFinite(v)&&Number.isFinite(y[i]))line.push([v,y[i]]);else if(line.length){lines.push(line);line=[];}});if(line.length)lines.push(line);return lines;}
if(typeof module!=='undefined')module.exports={applyCoupledSurface,integralCSV,integralSegments};
if(typeof window!=='undefined'){
 const data=window.COUPLED_ENTROPY;
 if(window.ROCKET_CASE)applyCoupledSurface(window.ROCKET_CASE,data);
 if(document.getElementById('integralBranch')){
  const selector=document.getElementById('integralBranch'),slider=document.getElementById('position');let hit=[];
  Object.entries(data.branches).forEach(([k,b])=>{const option=document.createElement('option');option.value=k;option.textContent=b.label;selector.appendChild(option);});
  const fmt=(v,d=6)=>Number.isFinite(v)?Number(v.toPrecision(d)).toString():'未知';
  function plot(id,x,y,marker){const canvas=document.getElementById(id),w=canvas.clientWidth,h=canvas.clientHeight,d=window.devicePixelRatio||1;canvas.width=w*d;canvas.height=h*d;const ctx=canvas.getContext('2d');ctx.scale(d,d);
   const seg=integralSegments(x,y),flat=seg.flat();if(!flat.length){ctx.fillText('未知',60,60);return;}
   let lo=Math.min(...flat.map(p=>p[0])),hi=Math.max(...flat.map(p=>p[0])),bot=Math.min(...flat.map(p=>p[1])),top=Math.max(...flat.map(p=>p[1]));if(lo===hi)hi++;if(bot===top){bot--;top++;}
   const px=v=>65+(w-85)*(v-lo)/(hi-lo),py=v=>h-30-(h-50)*(v-bot)/(top-bot);ctx.font='12px system-ui';
   for(let k=0;k<5;k++){const xx=lo+(hi-lo)*k/4,yy=bot+(top-bot)*k/4;ctx.strokeStyle='#e4e9ee';ctx.beginPath();ctx.moveTo(px(xx),20);ctx.lineTo(px(xx),h-30);ctx.moveTo(65,py(yy));ctx.lineTo(w-20,py(yy));ctx.stroke();ctx.fillStyle='#596775';ctx.textAlign='center';ctx.fillText(fmt(xx,4),px(xx),h-10);ctx.textAlign='right';ctx.fillText(fmt(yy,4),60,py(yy)+4);}
   ctx.strokeStyle='#246689';ctx.lineWidth=2;seg.forEach(line=>{ctx.beginPath();line.forEach((p,i)=>i?ctx.lineTo(px(p[0]),py(p[1])):ctx.moveTo(px(p[0]),py(p[1])));ctx.stroke();});
   if(marker!==undefined&&Number.isFinite(y[marker])){ctx.fillStyle='#ac5728';ctx.beginPath();ctx.arc(px(x[marker]),py(y[marker]),5,0,Math.PI*2);ctx.fill();}
  }
  function ring(stations,index){const canvas=document.getElementById('ring'),w=canvas.clientWidth,h=canvas.clientHeight,d=window.devicePixelRatio||1;canvas.width=w*d;canvas.height=h*d;const ctx=canvas.getContext('2d');ctx.scale(d,d);
   const span=Math.max(...stations.map(s=>Math.hypot(s.wallXYZ[1],s.wallXYZ[2])))*1.25,scale=Math.min(w-100,h-65)/(2*span),cx=w/2,cy=h/2;
   hit=stations.map(s=>[cx+s.wallXYZ[1]*scale,cy-s.wallXYZ[2]*scale]);ctx.font='12px system-ui';ctx.fillStyle='#596775';ctx.textAlign='center';
   [-.03,0,.03].forEach(v=>{ctx.strokeStyle='#e4e9ee';ctx.beginPath();ctx.moveTo(cx+v*scale,20);ctx.lineTo(cx+v*scale,h-30);ctx.moveTo(40,cy-v*scale);ctx.lineTo(w-15,cy-v*scale);ctx.stroke();ctx.fillText(v.toFixed(2),cx+v*scale,h-10);ctx.fillText(v.toFixed(2),24,cy-v*scale);});
   ctx.fillText('Y / m',w-28,h-10);ctx.fillText('Z / m',25,15);ctx.strokeStyle='#b7c2cc';ctx.beginPath();hit.forEach((p,i)=>i?ctx.lineTo(...p):ctx.moveTo(...p));ctx.closePath();ctx.stroke();
   hit.forEach((p,i)=>{ctx.fillStyle=i===index?'#ac5728':'#246689';ctx.beginPath();ctx.arc(...p,i===index?6:3.5,0,Math.PI*2);ctx.fill();});
  }
  function render(){const b=data.branches[selector.value],s=b.stations[+slider.value],i=+slider.value,phi=b.stations.map(s=>s.phiDeg);
   document.getElementById('positionLabel').textContent=`入口 ${i+1}/46 · φ ${fmt(s.phiDeg,5)}° · θ ${fmt(s.thetaDeg,5)}°`;
   document.getElementById('probe').textContent=`原上游面 ${s.upstreamFaceID} → 前体面 ${s.downstreamFaceID}；p=${fmt(s.pressurePa)} Pa；q=${fmt(s.heatFluxWm2/1000)} kW/m²；Te=${fmt(s.edgeTemperatureK)} K；Ue=${fmt(s.edgeSpeedMps)} m/s；质量厚度=${fmt(s.massThicknessM*1000)} mm。给定Tw=300 K；固体Tw未知。`;
   document.getElementById('source').textContent=`β=${fmt(s.betaPressureProxy)}；拟合hφ=${fmt(s.fittedMetricM)} m，原环hφ=${fmt(s.actualMetricM)} m。${s.oldProfileAvailable?'旧固定熵剖面有值，未与新外缘重新求解':'旧剖面缺失，本轮仅补出积分量'}。实际剖面/SA未知，CFD入口未就绪。源 ${b.sourceSHA256}。`;
   ring(b.stations,i);[['q','heatFluxWm2',.001],['thickness','massThicknessM',1000],['temperature','edgeTemperatureK',1],['speed','edgeSpeedMps',1],['beta','betaPressureProxy',1]].forEach(([id,key,c])=>plot(id,phi,b.stations.map(s=>s[key]*c),i));
   plot('streamwise',b.along.angleDeg,b.along.heatFluxWm2.map(v=>v/1000));
  }
  document.getElementById('ring').onclick=e=>{const r=e.currentTarget.getBoundingClientRect(),x=e.clientX-r.left,y=e.clientY-r.top;slider.value=hit.reduce((best,p,i)=>Math.hypot(p[0]-x,p[1]-y)<Math.hypot(hit[best][0]-x,hit[best][1]-y)?i:best,0);render();};
  document.getElementById('download').onclick=()=>{const url=URL.createObjectURL(new Blob([integralCSV(data,selector.value)],{type:'text/csv;charset=utf-8'})),a=document.createElement('a');a.href=url;a.download=`v254_${selector.value}_46_inlets.csv`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
  slider.oninput=render;selector.onchange=render;window.addEventListener('resize',render);render();
 }
}
