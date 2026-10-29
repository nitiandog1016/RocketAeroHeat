/* 模块用途：真实保存的流线/流管研究，共用五步骤、三维站点拾取和定量曲线。 */
'use strict';
(function(){
 const study=window.PRESSURE_SURFACE_FLOW_STUDY,$=id=>document.getElementById(id),buttons=[...document.querySelectorAll('[data-step]')];
 if(!study||study.kind!=='pressure-surface-flow-study-v1'||study.productionQualified!==false)throw Error('流线研究身份无效。');
 let selected=2,tube=0,station=0,step=4,plotBounds=null;
 const fmt=x=>Number.isFinite(x)?Number(x.toPrecision(6)).toString():'未求解';
 const labels=['输入模型','生成网格','添加边界条件','计算','查看结果'];
 const fields={metricRatio:{label:'条件数值 h/h₀',unit:'1',factor:1},metric:{label:'条件数值 h',unit:null,factor:1},pressureInputPa:{label:'输入压力',unit:'kPa',factor:.001},temperatureInputK:{label:'输入外缘温度',unit:'K',factor:1},speedInputMS:{label:'输入速度幅值',unit:'m/s',factor:1}};
 study.cases.forEach((c,i)=>{const option=document.createElement('option');option.value=String(i);option.textContent=c.caption;$('caseChoice').append(option);});
 $('caseChoice').value=String(selected);
 function appendRows(node,rows){node.replaceChildren();rows.forEach(([name,value])=>{const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=name;dd.textContent=String(value);node.append(dt,dd);});}
 function current(){return study.cases[selected];}
 function show(){
  const c=current(),summary=c.summary;$('title').textContent=labels[step]+' · '+c.caption;
  $('description').textContent=c.scope;
  const rows=[['输入压力/热状态','内部解析验证输入，未由流线模块求解'],['流线方向','曲面切向Euler压力转向 + RK4'],['横向度量','两组标签间距的横向投影差分'],['入口标签单位',c.labelUnit],['流管标签组',c.tubes.length],['实际有效站点',summary.validStations+' / '+summary.requestedStations],['几何范围',c.geometryScope],['热流与求解壁温','未求解'],['新增合格热流面积','0 m²']];
  if(step===1)rows.splice(0,0,['原三角面',c.geometry.faces],['几何闭合',c.geometry.closed?'闭合；每边2面且方向一致':'平面验证片，明确非闭合箭体']);
  if(step===2)rows.splice(0,0,['给定总温 / 总压',fmt(c.totalTemperatureK)+' K / '+fmt(c.totalPressurePa/1000)+' kPa'],['入口方向来源',c.sources.entryDirection]);
  if(step===3)rows.splice(0,0,['实际求解状态','已保存的条件数值积分；切换不重算'],['最大单位切向漂移',fmt(summary.maximumUnitTangentError)],['标签差值限制','分辨率指标≤'+fmt(c.maximumRelativeMetricDifference)+' h；不是严格误差上界']);
  appendRows($('details'),rows);$('scope').textContent=c.methodScope;$('sources').textContent=Object.entries(c.sources).map(([k,v])=>k+'：'+v).join('\n');
  buttons.forEach((b,i)=>b.setAttribute('aria-selected',String(i===step)));
  const url=study.geometryURL+(study.geometryURL.includes('?')?'&':'?')+'case='+selected+'&flow='+(step>=3?'1':'0')+'&wire='+(step===1?'1':'0');
  if($('viewer').getAttribute('src')!==url)$('viewer').setAttribute('src',url);
  if($('viewer').parentElement)$('viewer').parentElement.scrollTop=0;updateProbe();draw();
 }
 function updateProbe(){
  const c=current(),t=c.tubes[tube],n=t.arcLengthM.length;station=Math.max(0,Math.min(n-1,station));
  $('station').max=String(n-1);$('station').value=String(station);
  const selectedPressure=t.pressureInputPa[station];
  $('stationNote').textContent='流管 '+(tube+1)+' · s = '+fmt(t.arcLengthM[station])+' m · 数值h '+fmt(t.metric[station])+' '+c.metricUnit+' · 输入p '+fmt(Number.isFinite(selectedPressure)?selectedPressure*.001:NaN)+' kPa · 输入T '+fmt(t.temperatureInputK[station])+' K · q未求';
  const vector=x=>Array.isArray(x)?x.map(fmt).join(', '):'未求解';
  const pressure=t.pressureInputPa[station];
  appendRows($('probe'),[['实际积分站点',(station+1)+' / '+n],['XYZ / m',vector(t.position[station])],['数值单位流向',vector(t.tangent[station])],['数值横向h / '+c.metricUnit,fmt(t.metric[station])],['数值h/h₀',fmt(t.metricRatio[station])],['标签差值指标 / '+c.metricUnit,fmt(t.metricDifferenceIndicator[station])+'（非严格误差上界）'],['输入压力 / kPa',fmt(Number.isFinite(pressure)?pressure*.001:NaN)],['输入外缘温度 / K',fmt(t.temperatureInputK[station])],['输入密度 / kg/m³',fmt(t.densityInputKgM3[station])],['输入速度幅值 / m/s',fmt(t.speedInputMS[station])],['热流 / 求解壁温','未求解；不能由轨迹或输入温度推得'],['算法与状态',c.methodScope]]);
 }
 function draw(){
  const canvas=$('curve'),ctx=canvas.getContext('2d');if(!ctx)return;
  const dpr=Math.min(window.devicePixelRatio||1,2),w=Math.max(350,canvas.clientWidth),h=Math.max(120,canvas.clientHeight);canvas.width=w*dpr;canvas.height=h*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,w,h);
  const c=current(),t=c.tubes[tube],key=$('field').value,field=fields[key],unit=field.unit||c.metricUnit;
  const sameUnit=x=>key!=='metric'||x.metricUnit===c.metricUnit;
  const banks=$('sharedScale').checked?study.cases.filter(sameUnit):[c];
  const values=banks.flatMap(x=>x.tubes.flatMap(y=>y[key])).filter(Number.isFinite).map(x=>x*field.factor);
  let min=Math.min(...values),max=Math.max(...values);if(min===max){const pad=Math.max(Math.abs(min)*.05,.01);min-=pad;max+=pad;}
  const left=68,right=w-18,top=20,bottom=h-27,xmax=t.arcLengthM.at(-1),px=x=>left+x/xmax*(right-left),py=y=>bottom-(y-min)/(max-min)*(bottom-top);
  ctx.font='11px Segoe UI,Microsoft YaHei,sans-serif';ctx.fillStyle='#536879';ctx.fillText(field.label+' / '+unit,8,12);
  for(let i=0;i<3;i++){const y=min+(max-min)*i/2;ctx.strokeStyle='#e0e6eb';ctx.beginPath();ctx.moveTo(left,py(y));ctx.lineTo(right,py(y));ctx.stroke();ctx.fillText(fmt(y),5,py(y)+4);}
  for(let i=0;i<5;i++){const x=xmax*i/4;ctx.fillText(fmt(x),px(x)-8,h-8);}ctx.fillText('s / m',right-32,h-8);
  ctx.strokeStyle='#176b85';ctx.lineWidth=1.7;ctx.beginPath();let started=false;
  t[key].forEach((value,i)=>{if(!Number.isFinite(value)){started=false;return;}const x=px(t.arcLengthM[i]),y=py(value*field.factor);if(started)ctx.lineTo(x,y);else ctx.moveTo(x,y);started=true;});ctx.stroke();
  const value=t[key][station];if(Number.isFinite(value)){ctx.fillStyle='#176b85';ctx.beginPath();ctx.arc(px(t.arcLengthM[station]),py(value*field.factor),3,0,Math.PI*2);ctx.fill();}
  plotBounds={left,right,xmax};$('plotNote').textContent='三维点击可见流线站点，或点击曲线/拖动站点滑块；'+field.label+'，只显示实际有效积分前缀。横向h取n×t投影，不取两线总距离。';
 }
 buttons.forEach((b,i)=>b.addEventListener('click',()=>{step=i;show();}));
 $('caseChoice').addEventListener('change',()=>{selected=Number($('caseChoice').value);tube=0;station=0;show();});
 $('field').addEventListener('change',draw);$('sharedScale').addEventListener('change',draw);
 $('station').addEventListener('input',()=>{station=Number($('station').value);updateProbe();draw();});
 $('curve').addEventListener('click',e=>{if(!plotBounds)return;const r=$('curve').getBoundingClientRect(),x=(e.clientX-r.left-plotBounds.left)/(plotBounds.right-plotBounds.left)*plotBounds.xmax,t=current().tubes[tube];station=t.arcLengthM.reduce((best,v,i)=>Math.abs(v-x)<Math.abs(t.arcLengthM[best]-x)?i:best,0);updateProbe();draw();});
 window.addEventListener('message',e=>{
  if(e.origin!==location.origin||e.source!==$('viewer').contentWindow||!e.data||e.data.kind!=='pressure-surface-view-event'||e.data.caseIndex!==selected)return;
  const p=e.data.payload;if(!p||p.revision!==selected+1)return;
  if(p.kind==='geometry-error'){$('plotNote').textContent='三维显示失败：'+p.message;return;}
  if(p.kind!=='surface-flow-probe'||!Number.isInteger(p.curveIndex)||!Number.isInteger(p.stationIndex))return;
  const c=current(),curve=c.curves[p.curveIndex];if(!curve)return;const next=curve.sourceTubeIndexOneBased-1;
  if(next<0||next>=c.tubes.length||p.stationIndex<0||p.stationIndex>=c.tubes[next].arcLengthM.length)return;
  tube=next;station=p.stationIndex;updateProbe();draw();
 });
 $('saveCSV').addEventListener('click',()=>{
  const c=current(),t=c.tubes[tube],quote=x=>'"'+(x==null?'':String(x)).replace(/"/g,'""')+'"';
  const names=['s_m','x_m','y_m','z_m','direction_x','direction_y','direction_z','conditional_metric','metric_unit','conditional_metric_ratio','metric_difference_not_rigorous_bound','input_pressure_Pa','input_temperature_K','input_density_kg_m3','input_speed_m_s','heat_flux_status','wall_solved_status','method_scope'];
  const lines=[names.map(quote).join(',')];t.arcLengthM.forEach((s,i)=>lines.push([s,...t.position[i],...t.tangent[i],t.metric[i],c.metricUnit,t.metricRatio[i],t.metricDifferenceIndicator[i],t.pressureInputPa[i],t.temperatureInputK[i],t.densityInputKgM3[i],t.speedInputMS[i],'unknown','unknown',c.methodScope].map(quote).join(',')));
  const blob=new Blob(['\uFEFF'+lines.join('\n')+'\n'],{type:'text/csv;charset=utf-8'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='surface_flow_'+c.id+'_tube'+(tube+1)+'.csv';a.click();URL.revokeObjectURL(url);
 });
 if(typeof ResizeObserver==='function')new ResizeObserver(draw).observe($('curve'));show();
})();
