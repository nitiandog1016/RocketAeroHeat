/* 模块用途：同一五步骤工作台显示流线、横向尺度及逐变量的输入/条件反算身份。 */
'use strict';
(function(){
 const study=window.PRESSURE_SURFACE_FLOW_STUDY,$=id=>document.getElementById(id),buttons=[...document.querySelectorAll('[data-step]')];
 if(!study||study.kind!=='pressure-surface-flow-study-v1'||study.productionQualified!==false)throw Error('流线研究身份无效。');
 let selected=Math.min(2,study.cases.length-1),tube=0,station=0,step=4,plotBounds=null;
 const fmt=x=>Number.isFinite(x)?Number(x.toPrecision(6)).toString():'未求解';
 const scaled=(x,factor)=>Number.isFinite(x)?x*factor:NaN;
 const labels=['输入模型','生成网格','添加边界条件','计算','查看结果'];
 const fields={metricRatio:{label:'条件数值 h/h₀',unit:'1',factor:1},metric:{label:'条件数值 h',unit:null,factor:1},pressureInputPa:{label:'输入压力',unit:'kPa',factor:.001},temperatureInputK:{label:'输入外缘温度',unit:'K',factor:1},speedInputMS:{label:'输入速度幅值',unit:'m/s',factor:1}};
 const restored=c=>c.thermodynamicStateIsConditionallyRestored===true;
 const labelled=c=>c.inletVariesByLabel===true;
 const aliases={temperatureInputK:'temperatureRecoveredK',densityInputKgM3:'densityRecoveredKgM3',speedInputMS:'speedRecoveredMS'};
 const temperatureName=c=>restored(c)?'条件反算外缘温度Te':'输入外缘温度';
 const speedName=c=>restored(c)?'条件反算速度幅值':'输入速度幅值';
 function series(c,t,key){const a=t[restored(c)&&aliases[key]?aliases[key]:key];if(!Array.isArray(a))throw Error('缺少声明场量序列：'+key);return a;}
 const hasSeries=(c,t,key)=>Array.isArray(t[restored(c)&&aliases[key]?aliases[key]:key]);
 const extraOptions=[];
 for(const [key,field] of Object.entries({totalEnthalpyJPerKg:{label:'沿线总焓H₀',unit:'kJ/kg',factor:.001},entropyJPerKgK:{label:'沿线熵',unit:'J/(kg·K)',factor:1},entropyGradientNormJPerKgKM:{label:'物理熵梯度幅值',unit:'J/(kg·K·m)',factor:1}})){
  if(!study.cases.some(c=>c.tubes.every(t=>hasSeries(c,t,key))))continue;
  fields[key]=field;const option=document.createElement('option');option.value=key;option.textContent=field.label+' / '+field.unit;$('field').append(option);extraOptions.push([key,option]);
 }
 function fieldLabel(c,key){if(key==='temperatureInputK')return temperatureName(c);if(key==='speedInputMS')return speedName(c);return fields[key].label;}
 study.cases.forEach((c,i)=>{
  if(restored(c)&&c.thermodynamicStateIsInput!==false)throw Error('条件反算状态不能同时标为给定外缘。');
  const option=document.createElement('option');option.value=String(i);option.textContent=c.caption;$('caseChoice').append(option);
 });
 $('caseChoice').value=String(selected);
 function appendRows(node,rows){node.replaceChildren();rows.forEach(([name,value])=>{const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=name;dd.textContent=String(value);node.append(dt,dd);});}
 function current(){return study.cases[selected];}

 // 模块一：阶段切换保持同一实际数据；压力来源与外缘反算的身份分别呈现。
 function show(){
  const c=current(),summary=c.summary;$('title').textContent=labels[step]+' · '+c.caption;$('description').textContent=c.scope;
  const role=labelled(c)?'连续压力与逐入口标签H₀/S₀条件反算；不是独立完整外流解':restored(c)?'连续压力与固定总状态的同物性条件反算；不是独立完整外流解':'内部解析验证输入，未由流线模块求解';
  const rows=[['压力来源（流线输入）',c.pressureSourceRole||'内部解析验证压力'],['外缘热状态',role],['流线方向','曲面切向Euler压力转向 + RK4'],['横向度量','两组标签间距的横向投影差分'],['入口标签单位',c.labelUnit],['流管标签组',c.tubes.length],['实际有效站点',summary.validStations+' / '+summary.requestedStations],['几何范围',c.geometryScope],['热流与求解壁温','未求解'],['新增合格热流面积','0 m²']];
  if(c.thermalConventionCaption)rows.splice(2,0,['热物性约定',c.thermalConventionCaption]);
  if(step===1)rows.splice(0,0,['原三角面',c.geometry.faces],['几何闭合',c.geometry.closed?'闭合；每边2面且方向一致':'平面验证片，明确非闭合箭体']);
  if(labelled(c))rows.splice(3,0,['入口总状态','随标签变化；每条流线独立守恒，沿程不重置'],['物理横向梯度','相邻入口标签 + n×t横向尺度；辅助延拓导数不作物理场']);
  if(step===2){
   if(labelled(c))rows.splice(0,0,['入口总温范围 / K',c.entryTotalTemperatureRangeK.map(fmt).join(' ～ ')],['入口总压范围 / kPa',c.entryTotalPressureRangePa.map(x=>fmt(x/1000)).join(' ～ ')],['入口H₀/S₀来源',c.sources.entryProfile]);
   else rows.splice(0,0,['固定总温 / 总压',fmt(c.totalTemperatureK)+' K / '+fmt(c.totalPressurePa/1000)+' kPa']);
   rows.splice(0,0,['入口方向来源',c.sources.entryDirection]);
  }
  if(step===3)rows.splice(0,0,['实际求解状态','已保存的条件数值积分；切换不重算'],['最大单位切向漂移',fmt(summary.maximumUnitTangentError)],['标签差值限制','分辨率指标≤'+fmt(c.maximumRelativeMetricDifference)+' h；不是严格误差上界']);
  appendRows($('details'),rows);$('scope').textContent=c.methodScope;$('sources').textContent=Object.entries(c.sources).map(([k,v])=>k+'：'+v).join('\n');
  if($('roleNote'))$('roleNote').textContent=labelled(c)?'各入口标签分别携带总焓与熵；中心物理两向梯度由相邻标签重建，单射线辅助延拓导数不作物理场。这里的Te不是绝热恢复温度Tr，也不是壁温Tw；热流和求解壁温仍未知。':restored(c)?'压力来源单独列明；Te、密度和速度由同一热物性与固定总状态条件反算。这里的Te不是绝热恢复温度Tr，也不是壁温Tw；热流和求解壁温仍未知，不能据流线数量增加合格面积。':'压力、外缘温度、密度和速度幅值为解析验证输入。数值流线与横向尺度不构成完整外流、热流或壁温求解，不能按着色增加合格面积。';
  extraOptions.forEach(([key,option])=>{option.disabled=!c.tubes.every(t=>hasSeries(c,t,key));});
  if(!c.tubes.every(t=>hasSeries(c,t,$('field').value)))$('field').value='metricRatio';
  buttons.forEach((b,i)=>b.setAttribute('aria-selected',String(i===step)));
  const url=study.geometryURL+(study.geometryURL.includes('?')?'&':'?')+'case='+selected+'&flow='+(step>=3?'1':'0')+'&wire='+(step===1?'1':'0');
  if($('viewer').getAttribute('src')!==url)$('viewer').setAttribute('src',url);
  if($('viewer').parentElement)$('viewer').parentElement.scrollTop=0;updateProbe();draw();
 }

 // 模块二：探针保留原坐标/数值；Te、Tr、给定Tw与求解Tw分开。
 function updateProbe(){
  const c=current(),t=c.tubes[tube],n=t.arcLengthM.length;station=Math.max(0,Math.min(n-1,station));
  $('station').max=String(n-1);$('station').value=String(station);
  const pressure=t.pressureInputPa[station],temperature=series(c,t,'temperatureInputK')[station];
  $('stationNote').textContent='流管 '+(tube+1)+' · s = '+fmt(t.arcLengthM[station])+' m · 数值h '+fmt(t.metric[station])+' '+c.metricUnit+' · 输入p '+fmt(Number.isFinite(pressure)?pressure*.001:NaN)+' kPa · '+(restored(c)?'反算Te ':'输入T ')+fmt(temperature)+' K · q未求';
  const vector=x=>Array.isArray(x)?x.map(fmt).join(', '):'未求解';
  const rows=[['实际积分站点',(station+1)+' / '+n],['XYZ / m',vector(t.position[station])],['数值单位流向',vector(t.tangent[station])],['数值横向h / '+c.metricUnit,fmt(t.metric[station])],['数值h/h₀',fmt(t.metricRatio[station])],['标签差值指标 / '+c.metricUnit,fmt(t.metricDifferenceIndicator[station])+'（非严格误差上界）'],['输入压力 / kPa',fmt(Number.isFinite(pressure)?pressure*.001:NaN)],[temperatureName(c)+' / K',fmt(temperature)],[(restored(c)?'条件反算密度':'输入密度')+' / kg/m³',fmt(series(c,t,'densityInputKgM3')[station])],[speedName(c)+' / m/s',fmt(series(c,t,'speedInputMS')[station])],['恢复温度Tr / K','未计算'],['给定壁温Tw / K','本研究未设置'],['热流 / 求解壁温','未求解；不能由轨迹或反算Te推得'],['压力方法',c.sources.pressure],['外缘状态方法',c.sources.thermodynamics],['算法与状态',c.methodScope]];
  if(labelled(c))rows.splice(10,0,['入口标签 / '+c.labelUnit,fmt(t.labelCenter)],['该入口总温 / 总压',fmt(t.inletTotalTemperatureK)+' K / '+fmt(scaled(t.inletTotalPressurePa,.001))+' kPa'],['沿线总焓H₀ / kJ/kg',fmt(scaled(t.totalEnthalpyJPerKg[station],.001))],['沿线熵 / J/(kg·K)',fmt(t.entropyJPerKgK[station])],['物理gradH₀ / J/(kg·m)',vector(t.totalEnthalpyGradient[station])],['物理gradS₀ / J/(kg·K·m)',vector(t.entropyGradient[station])],['物理gradTe / K/m',vector(t.temperatureGradient[station])],['横向梯度方法',c.sources.transverseGradients]);
  if(labelled(c)){const r=t.inletDerivativeResolution;rows.splice(14,0,['入口导数分辨状态','总焓变化'+(r.enthalpyVariationResolved?'已分辨':'未分辨（低于源状态一致性尺度）')+'；熵变化'+(r.entropyVariationResolved?'已分辨':'未分辨（低于源状态一致性尺度）')+'；一致性尺度不是严格误差界']);}
  appendRows($('probe'),rows);
 }

 // 模块三：同量纲曲线可共用纵轴；断值不补零，条件反算场量使用正确标签。
 function draw(){
  const canvas=$('curve'),ctx=canvas.getContext('2d');if(!ctx)return;
  const dpr=Math.min(window.devicePixelRatio||1,2),w=Math.max(350,canvas.clientWidth),h=Math.max(120,canvas.clientHeight);canvas.width=w*dpr;canvas.height=h*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,w,h);
  const c=current(),t=c.tubes[tube],key=$('field').value,field=fields[key],unit=field.unit||c.metricUnit;
  const banks=($('sharedScale').checked?study.cases.filter(x=>key!=='metric'||x.metricUnit===c.metricUnit):[c]).filter(x=>x.tubes.every(y=>hasSeries(x,y,key)));
  const values=banks.flatMap(x=>x.tubes.flatMap(y=>series(x,y,key))).filter(Number.isFinite).map(x=>x*field.factor);
  if(!values.length){ctx.fillText('当前物理量未求解',8,24);plotBounds=null;return;}
  let min=Math.min(...values),max=Math.max(...values);if(max-min<=1e-10*Math.max(1,Math.abs(min),Math.abs(max))){const pad=Math.max(Math.abs(min)*.05,.01);min-=pad;max+=pad;}
  const left=68,right=w-18,top=20,bottom=h-27,xmax=t.arcLengthM.at(-1),px=x=>left+x/xmax*(right-left),py=y=>bottom-(y-min)/(max-min)*(bottom-top);
  ctx.font='11px Segoe UI,Microsoft YaHei,sans-serif';ctx.fillStyle='#536879';ctx.fillText(fieldLabel(c,key)+' / '+unit,8,12);
  for(let i=0;i<3;i++){const y=min+(max-min)*i/2;ctx.strokeStyle='#e0e6eb';ctx.beginPath();ctx.moveTo(left,py(y));ctx.lineTo(right,py(y));ctx.stroke();ctx.fillText(fmt(y),5,py(y)+4);}
  for(let i=0;i<5;i++){const x=xmax*i/4;ctx.fillText(fmt(x),px(x)-8,h-8);}ctx.fillText('s / m',right-32,h-8);
  ctx.strokeStyle='#176b85';ctx.lineWidth=1.7;ctx.beginPath();let started=false;
  const displayed=series(c,t,key);displayed.forEach((value,i)=>{if(!Number.isFinite(value)){started=false;return;}const x=px(t.arcLengthM[i]),y=py(value*field.factor);if(started)ctx.lineTo(x,y);else ctx.moveTo(x,y);started=true;});ctx.stroke();
  const value=displayed[station];if(Number.isFinite(value)){ctx.fillStyle='#176b85';ctx.beginPath();ctx.arc(px(t.arcLengthM[station]),py(value*field.factor),3,0,Math.PI*2);ctx.fill();}
  plotBounds={left,right,xmax};$('plotNote').textContent='三维点击可见流线站点，或点击曲线/拖动站点滑块；'+fieldLabel(c,key)+'，只显示实际有效积分前缀。横向h取n×t投影，不取两线总距离。';
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

 // 模块四：导出列名明确区分输入和条件反算；所有未求热量保持unknown。
 $('saveCSV').addEventListener('click',()=>{
  const c=current(),t=c.tubes[tube],quote=x=>'"'+(x==null?'':String(x)).replace(/"/g,'""')+'"';
  const names=['s_m','x_m','y_m','z_m','direction_x','direction_y','direction_z','conditional_metric','metric_unit','conditional_metric_ratio','metric_difference_not_rigorous_bound','input_pressure_Pa',restored(c)?'conditionally_recovered_edge_temperature_K':'input_temperature_K',restored(c)?'conditionally_recovered_density_kg_m3':'input_density_kg_m3',restored(c)?'conditionally_recovered_speed_m_s':'input_speed_m_s','heat_flux_status','wall_solved_status','thermal_state_role','thermodynamic_convention','pressure_method','edge_state_method','direction_metric_method','method_scope'];
  if(labelled(c))names.push('entry_label','entry_label_unit','inlet_total_pressure_Pa','inlet_total_temperature_K','total_enthalpy_J_kg','entropy_J_kg_K','gradH_x_J_kg_m','gradH_y_J_kg_m','gradH_z_J_kg_m','gradS_x_J_kg_K_m','gradS_y_J_kg_K_m','gradS_z_J_kg_K_m','gradTe_x_K_m','gradTe_y_K_m','gradTe_z_K_m','entry_state_method','physical_gradient_method','derivative_role');
  if(labelled(c))names.push('inlet_H_derivative_difference','inlet_S_derivative_difference','inlet_H_consistency_scale','inlet_S_consistency_scale','inlet_H_variation_resolved','inlet_S_variation_resolved','inlet_consistency_scale_is_rigorous_error_bound');
  const lines=[names.map(quote).join(',')];t.arcLengthM.forEach((s,i)=>{
   const values=[s,...t.position[i],...t.tangent[i],t.metric[i],c.metricUnit,t.metricRatio[i],t.metricDifferenceIndicator[i],t.pressureInputPa[i],series(c,t,'temperatureInputK')[i],series(c,t,'densityInputKgM3')[i],series(c,t,'speedInputMS')[i],'unknown','unknown',restored(c)?'conditional-recovery':'given-edge-input',c.thermodynamicConvention||'calorically-perfect-v94',c.sources.pressure,c.sources.thermodynamics,'曲面切向Euler压力转向 + RK4；五射线双间距横向差分',c.methodScope];
   if(labelled(c))values.push(t.labelCenter,c.labelUnit,t.inletTotalPressurePa,t.inletTotalTemperatureK,t.totalEnthalpyJPerKg[i],t.entropyJPerKgK[i],...t.totalEnthalpyGradient[i],...t.entropyGradient[i],...t.temperatureGradient[i],c.sources.entryProfile,c.sources.transverseGradients,'physical-surface-gradient-from-inlet-labels');
   if(labelled(c)){const r=t.inletDerivativeResolution;values.push(r.enthalpyDifferenceIndicator,r.entropyDifferenceIndicator,r.enthalpyConsistencyScale,r.entropyConsistencyScale,r.enthalpyVariationResolved,r.entropyVariationResolved,r.consistencyScalesAreRigorousErrorBounds);}
   lines.push(values.map(quote).join(','));
  });
  const blob=new Blob(['\uFEFF'+lines.join('\n')+'\n'],{type:'text/csv;charset=utf-8'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='surface_flow_'+c.id+'_tube'+(tube+1)+'.csv';a.click();URL.revokeObjectURL(url);
 });
 if(typeof ResizeObserver==='function')new ResizeObserver(draw).observe($('curve'));show();
})();
