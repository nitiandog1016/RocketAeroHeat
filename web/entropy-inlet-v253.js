/* 模块职责：真实入口选择、缺失段显示及准确剖面CSV。浏览不触发计算。 */
'use strict';
function inletCSV(s){
  if(!s.available)throw new Error('INLET_PROFILE_UNAVAILABLE');
  const rows=['branch,edge_index,upstream_face_id,downstream_face_id,phi_deg,normal_y_m,x_m,y_m,z_m,u_m_s,v_m_s,w_m_s,T_K,rho_kg_m3,rho_Ephysical_J_m3,pressure_Pa,source_chain_sha256,profile_ready_for_CFD'];
  for(let j=0;j<s.profile.y.length;j++){
    const y=s.profile.y[j],speed=s.profile.speed[j],xyz=s.wallXYZ.map((v,k)=>v+y*s.normal[k]),vel=s.tangent.map(v=>v*speed);
    rows.push([s.branch,s.index,s.upstreamFaceID,s.downstreamFaceID,s.phiDeg,y,...xyz,...vel,s.profile.temperature[j],s.profile.density[j],s.profile.rhoE[j],s.pressurePa,s.sourceChainSHA256,false].join(','));
  }return rows.join('\n')+'\n';
}
function finiteSegments(x,y){
  const lines=[];let line=[];for(let i=0;i<x.length;i++){if(Number.isFinite(x[i])&&Number.isFinite(y[i]))line.push([x[i],y[i]]);else if(line.length){lines.push(line);line=[];}}if(line.length)lines.push(line);return lines;
}
if(typeof module!=='undefined')module.exports={inletCSV,finiteSegments};
if(typeof window!=='undefined'){
  const data=window.ENTROPY_INLET,select=document.getElementById('branch'),slider=document.getElementById('station');let hit=[];
  for(const [key,b]of Object.entries(data.branches)){const o=document.createElement('option');o.value=key;o.textContent=b.label;select.appendChild(o);}
  function plot(id,x,y){const canvas=document.getElementById(id),w=canvas.clientWidth,h=canvas.clientHeight,d=window.devicePixelRatio||1;canvas.width=w*d;canvas.height=h*d;const ctx=canvas.getContext('2d');ctx.scale(d,d);
    const seg=finiteSegments(x,y),flat=seg.flat();if(!flat.length){ctx.fillStyle='#596775';ctx.font='15px system-ui';ctx.fillText('缺少已求得的上游剖面',65,h/2);return;}
    let lo=Math.min(...flat.map(p=>p[0])),hi=Math.max(...flat.map(p=>p[0])),bot=Math.min(...flat.map(p=>p[1])),top=Math.max(...flat.map(p=>p[1]));if(lo===hi)hi++;if(bot===top){bot--;top++;}
    const px=v=>62+(w-82)*(v-lo)/(hi-lo),py=v=>h-30-(h-50)*(v-bot)/(top-bot);ctx.font='12px system-ui';
    for(let k=0;k<5;k++){const xx=lo+(hi-lo)*k/4,yy=bot+(top-bot)*k/4;ctx.strokeStyle='#e5eaef';ctx.beginPath();ctx.moveTo(px(xx),20);ctx.lineTo(px(xx),h-30);ctx.moveTo(62,py(yy));ctx.lineTo(w-20,py(yy));ctx.stroke();ctx.fillStyle='#596775';ctx.textAlign='center';ctx.fillText(Number(xx.toPrecision(4)).toString(),px(xx),h-10);ctx.textAlign='right';ctx.fillText(Number(yy.toPrecision(4)).toString(),56,py(yy)+4);}
    ctx.strokeStyle='#246689';ctx.lineWidth=2;for(const line of seg){ctx.beginPath();for(let i=0;i<line.length;i++){const p=line[i];if(i)ctx.lineTo(px(p[0]),py(p[1]));else ctx.moveTo(px(p[0]),py(p[1]));}ctx.stroke();}
  }
  function ring(records,index){const c=document.getElementById('ring'),w=c.clientWidth,h=c.clientHeight,d=window.devicePixelRatio||1;c.width=w*d;c.height=h*d;const ctx=c.getContext('2d');ctx.scale(d,d);
    const range=Math.max(...records.map(s=>Math.hypot(s.wallXYZ[1],s.wallXYZ[2])))*1.18,scale=Math.min(w-100,h-70)/(2*range),cx=w/2,cy=h/2;
    const project=s=>[cx+s.wallXYZ[1]*scale,cy-s.wallXYZ[2]*scale];hit=records.map(project);ctx.font='12px system-ui';ctx.fillStyle='#596775';ctx.textAlign='center';
    for(const v of [-.03,0,.03]){ctx.strokeStyle='#e5eaef';ctx.beginPath();ctx.moveTo(cx+v*scale,30);ctx.lineTo(cx+v*scale,h-30);ctx.moveTo(45,cy-v*scale);ctx.lineTo(w-20,cy-v*scale);ctx.stroke();ctx.fillText(v.toFixed(2),cx+v*scale,h-10);ctx.fillText(v.toFixed(2),25,cy-v*scale);}
    ctx.fillText('Y / m',w-38,h-10);ctx.fillText('Z / m',25,16);ctx.strokeStyle='#c7ced5';ctx.beginPath();hit.forEach((p,i)=>{if(i)ctx.lineTo(...p);else ctx.moveTo(...p);});ctx.closePath();ctx.stroke();
    hit.forEach((p,i)=>{ctx.fillStyle=i===index?'#bd632c':records[i].available?'#246689':'#b8bec5';ctx.beginPath();ctx.arc(...p,i===index?6:3.5,0,2*Math.PI);ctx.fill();});
  }
  function render(){const records=data.branches[select.value].stations,i=Number(slider.value),s=records[i];document.getElementById('position').textContent=`φ ${s.phiDeg.toFixed(2)}° · 入口 ${i+1}/46`;
    const button=document.getElementById('csv');button.disabled=!s.available;
    document.getElementById('status').textContent=s.available?`前体面 ${s.downstreamFaceID}；原外缘 Te=${s.sourceEdgeTemperatureK.toFixed(2)} K，Ue=${s.sourceEdgeSpeedMps.toFixed(2)} m/s；匹配后 Te=${s.matchedEdge.temperatureK.toFixed(2)} K，Ue=${s.matchedEdge.speedMps.toFixed(2)} m/s；给定 Tw=${s.wallTemperatureK.toFixed(2)} K。入口未就绪。`:`前体面 ${s.downstreamFaceID}；${s.reason}。入口未就绪。`;
    document.getElementById('source').textContent=`交界 X=${s.wallXYZ[0].toPrecision(8)} m；h=${s.metricM.toPrecision(6)} m；风轴方位标签Jacobian=${s.labelJacobian.toPrecision(6)}。原剖面保持，变熵修正未反算回剖面。`;
    const phi=records.map(s=>s.phiDeg);ring(records,i);plot('metric',phi,records.map(s=>s.metricM/Math.hypot(s.wallXYZ[1],s.wallXYZ[2])));
    plot('edgeT',phi,records.map(s=>s.available?100*s.edgeTemperatureRelativeCorrection:null));plot('edgeU',phi,records.map(s=>s.available?100*s.edgeSpeedRelativeCorrection:null));
    plot('temperature',s.available?s.profile.temperature:[],s.available?s.profile.y.map(v=>1000*v):[]);plot('velocity',s.available?s.profile.speed:[],s.available?s.profile.y.map(v=>1000*v):[]);
  }
  document.getElementById('ring').onclick=e=>{const box=e.currentTarget.getBoundingClientRect(),x=e.clientX-box.left,y=e.clientY-box.top;slider.value=hit.reduce((best,p,i)=>Math.hypot(p[0]-x,p[1]-y)<Math.hypot(hit[best][0]-x,hit[best][1]-y)?i:best,0);render();};
  document.getElementById('csv').onclick=()=>{const s=data.branches[select.value].stations[Number(slider.value)],url=URL.createObjectURL(new Blob([inletCSV(s)],{type:'text/csv;charset=utf-8'})),a=document.createElement('a');a.href=url;a.download=`inlet_${s.branch}_${s.index+1}.csv`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
  slider.oninput=render;select.onchange=render;window.addEventListener('resize',render);render();
}
