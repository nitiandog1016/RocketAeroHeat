/* 模块职责：真实保存剖面的沿程选择、带单位曲线和可追溯CSV；不计算或填补物理场。 */
'use strict';
function profileCSV(s,branch){
  const rows=['branch,angle_deg,y_m,u_over_Ue,fluid_T_K,given_Tw_K,q_into_wall_W_m2,tau_Pa,source_sha256'];
  for(let i=0;i<s.y.length;i++)rows.push([branch,s.angleDeg,s.y[i],s.speedRatio[i],s.temperature[i],s.wallTemperature,s.heatFlux,s.wallShearStress,s.sourceChainSHA256].join(','));
  return rows.join('\n')+'\n';
}
if(typeof module!=='undefined')module.exports={profileCSV};
if(typeof window!=='undefined'){
  const data=window.PHASE_PROFILES,branch=document.getElementById('branch'),slider=document.getElementById('station');
  for(const [key,v] of Object.entries(data)){const o=document.createElement('option');o.value=key;o.textContent=v.label;branch.appendChild(o);}
  function plot(id,x,y,color){
    const c=document.getElementById(id),dpr=window.devicePixelRatio||1,w=c.clientWidth,h=c.clientHeight;c.width=w*dpr;c.height=h*dpr;
    const ctx=c.getContext('2d');ctx.scale(dpr,dpr);ctx.clearRect(0,0,w,h);
    let x0=Math.min(...x),x1=Math.max(...x),y0=Math.min(...y),y1=Math.max(...y);if(x1===x0)x1=x0+1;if(y1===y0){y0-=1;y1+=1;}
    const px=v=>58+(w-76)*(v-x0)/(x1-x0),py=v=>h-32-(h-50)*(v-y0)/(y1-y0);
    ctx.font='12px system-ui';ctx.lineWidth=1;
    for(let k=0;k<5;k++){const xx=x0+(x1-x0)*k/4,yy=y0+(y1-y0)*k/4;ctx.strokeStyle='#e5eaef';ctx.beginPath();ctx.moveTo(px(xx),18);ctx.lineTo(px(xx),h-32);ctx.moveTo(58,py(yy));ctx.lineTo(w-18,py(yy));ctx.stroke();ctx.fillStyle='#4c5965';ctx.textAlign='center';ctx.fillText(Number(xx.toPrecision(4)).toString(),px(xx),h-12);ctx.textAlign='right';ctx.fillText(Number(yy.toPrecision(4)).toString(),51,py(yy)+4);}
    ctx.strokeStyle=color;ctx.lineWidth=2;ctx.beginPath();for(let i=0;i<x.length;i++){if(i)ctx.lineTo(px(x[i]),py(y[i]));else ctx.moveTo(px(x[i]),py(y[i]));}ctx.stroke();
  }
  function render(){const b=data[branch.value];slider.max=b.stations.length-1;const s=b.stations[Number(slider.value)];
    document.getElementById('angle').textContent=s.angleDeg.toFixed(1)+'°';
    document.getElementById('values').textContent=`入壁热流 ${(s.heatFlux/1e6).toFixed(5)} MW/m²　壁剪应力 ${s.wallShearStress.toFixed(3)} Pa　给定壁温 ${s.wallTemperature.toFixed(2)} K　外缘温度 ${s.edgeTemperature.toFixed(2)} K　外缘速度 ${s.edgeSpeed.toFixed(2)} m/s`;
    document.getElementById('method').textContent=b.method+'；'+b.scope;
    plot('heat',b.stations.map(s=>s.angleDeg),b.stations.map(s=>s.heatFlux/1e6),'#25698c');
    plot('wall',b.stations.map(s=>s.angleDeg),b.stations.map(s=>s.wallTemperature),'#9d6633');
    plot('velocity',s.speedRatio,s.y.map(y=>y*1000),'#25698c');plot('temperature',s.temperature,s.y.map(y=>y*1000),'#9d6633');
  }
  let lastBranch=branch.value;
  branch.onchange=()=>{const angle=data[lastBranch].stations[Number(slider.value)].angleDeg,stations=data[branch.value].stations;
    slider.max=stations.length-1;slider.value=stations.reduce((best,s,i)=>Math.abs(s.angleDeg-angle)<Math.abs(stations[best].angleDeg-angle)?i:best,0);lastBranch=branch.value;render();};slider.oninput=render;window.addEventListener('resize',render);
  document.getElementById('csv').onclick=()=>{const s=data[branch.value].stations[Number(slider.value)],url=URL.createObjectURL(new Blob([profileCSV(s,branch.value)],{type:'text/csv;charset=utf-8'}));const a=document.createElement('a');a.href=url;a.download=`profile_${branch.value}_${s.angleDeg.toFixed(1)}deg.csv`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};render();
}
