/* 模块用途：同一侧栏查看保存的模型/网格/边界/计算/结果；逐面探针不补造未知值。 */
'use strict';
(function(){
 const data=window.ALE_STUDY,$=id=>document.getElementById(id),buttons=[...document.querySelectorAll('[data-step]')];
 if(!data||!['shock-fitted-kernel-study-v101','shock-fitted-dynamics-study-v102'].includes(data.kind)||data.productionQualified!==false||data.newProductionSolvedAreaM2!==0)throw Error('ALE研究身份无效。');
 const sectionEnabled=data.kind==='shock-fitted-dynamics-study-v102'&&Boolean(window.RocketSection);
 const fields={pressurePa:['流体压力','Pa',[0,40000]],temperatureK:['流体温度','K',[80,2000]],speedMS:['速度幅值','m/s',[0,2000]],wallPressurePa:['Riemann壁压','Pa',[0,40000]],heatFluxWm2:['热流','W/m²',[0,1e6]],wallTemperatureK:['求解壁温','K',[80,2000]],edgeTemperatureK:['外缘温度 Te','K',[80,2000]],recoveryTemperatureK:['恢复温度 Tr','K',[80,2000]],inputWallTemperatureK:['给定壁温','K',[80,2000]]},labels=['输入模型','生成网格','添加边界条件','计算','查看结果'];
 const fmt=v=>Number.isFinite(v)?Number(v.toPrecision(8)).toString():'未求解',kind=k=>['','储库侧面','移动/固定滑移壁面','移动激波 / 候选外壳','超声速全特征出流'][k]||'原箭体';
 let selected=Number.isSafeInteger(data.defaultCaseIndex)&&data.cases[data.defaultCaseIndex]?data.defaultCaseIndex:0,frame=1,step=4,view='all',field='pressurePa',face=0,revision=0,pickPosition=null;
 let sliceX=null,section=null,sectionSamples=[],plotSegments=[],lastOverlayKey='',meshCache=null;
 function rows(node,items){node.replaceChildren();items.forEach(([name,value])=>{const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=name;dd.textContent=String(value);node.append(dt,dd);});}
 function current(){return data.cases[selected];}
 function state(){const c=current(),f=c.frames[frame];if(view==='rocket')return{c,f,V:data.originalGeometry.vertices,F:data.originalGeometry.faces,map:data.originalGeometry.faces.map((_,i)=>i),full:true};const map=c.faces.map((_,i)=>i).filter(i=>view!=='wall'&&view!=='outer'||c.boundaryKind[i]===(view==='wall'?2:3));return{c,f,V:f.position,F:map.map(i=>c.faces[i]),map,full:false};}
 function value(s,key,index){return s.full||!Array.isArray(s.f[key])?null:s.f[key][s.map[index]];}
 function method(s,key){return s.full?'原箭体物理场未求':s.c.methods[key]||(key==='inputWallTemperatureK'?'本研究未设置壁温边界':'本研究未计算此场量');}
 function source(s){return s.c.source||data.source;}
 function range(s,key){return ['pressurePa','wallPressurePa'].includes(key)?s.c.pressureRangePa||data.pressureRangePa||fields[key][2]:fields[key][2];}
 function center(s,index){const ids=s.F[index];return [0,1,2].map(j=>ids.reduce((t,i)=>t+s.V[i][j],0)/3);}
 function show(){
  const c=current();if((view==='wall'&&!c.boundaryKind.includes(2))||(view==='outer'&&!c.boundaryKind.includes(3))){view='all';$('view').value=view;}
  const s=state();face=Math.max(0,Math.min(face,s.F.length-1));$('title').textContent=labels[step]+' · '+s.c.caption;$('scope').textContent=s.c.scope;$('gasScope').textContent=data.NASA7Scope;
  const detail=[['当前范围',s.F.length+' 个三角面'],['研究身份',s.c.stateRole],['数值核验',s.c.diagnostic],['正式新增热流面积','0 m²']];
  if(s.c.conditionLabel)detail.unshift(['工况',s.c.conditionLabel]);
  if(step===0)detail.unshift(['模型',s.c.modelRole||'保留原坐标与原面号；活塞为独立解析基准']);
  if(step===1)detail.unshift(['体网格','共享三角面的三维四面体，未复制轴对称剖面']);
  if(step===2)detail.unshift(['边界',s.c.boundaryDescription||'储库 / 滑移壁面 / 移动激波；壁压由内部状态相容求得'],['壁温边界','本研究未设置']);
  if(step===3)detail.unshift(['计算入口','查看已保存结果；此页不启动新求解'],['未完成项',s.c.unsolvedNote||'真实曲面激波与壁面同时闭合；黏性热流未求']);
  rows($('details'),detail);buttons.forEach((b,i)=>b.setAttribute('aria-selected',String(i===step)));
  $('frame').max=String(s.c.frames.length-1);$('frame').value=String(frame);$('time').textContent=(s.c.id==='piston'||s.c.timeKind==='transient'?'实际保存 t = ':'几何试验参数 t = ')+fmt(s.f.timeS)+' s · 仅实际保存帧';
  $('face').max=String(s.F.length-1);$('face').value=String(face);revision++;lastOverlayKey='';
  $('viewer').setAttribute('src','geometry/index.html?'+new URLSearchParams({case:String(selected),frame:String(frame),revision:String(revision),view,field,layer:step===4?'1':'0',wire:step===1||(step===4&&s.c.id==='piston')?'1':'0'}));
  $('renderState').textContent='加载三维网格';probe();if(sectionEnabled)updateSection();else draw();
 }
 function probe(){
  const s=state(),index=s.map[face],id=s.full?index+1:s.c.sourceFaceIDs[index],xyz=pickPosition||center(s,face);
  $('role').textContent=s.full?'原完整箭体几何；真实物理场尚未求解。':s.c.stateRole;
  rows($('probe'),[['当前面元',face+1+' / '+s.F.length],['坐标 / m',xyz.map(fmt).join(', ')],[s.full?'所属原箭面号':s.c.sourceFaceLabel||'所属原箭面号',id>0?id:'基准面元，无原箭体对应'],['边界身份',s.full?'原整箭未求解面':kind(s.c.boundaryKind[index])],...Object.entries(fields).map(([key,[name,unit]])=>[name+' / '+unit,key==='inputWallTemperatureK'&&!Number.isFinite(value(s,key,face))?'未设置':fmt(value(s,key,face))])]);
  rows($('methods'),[['数据SHA256',source(s).evidenceSHA256],['审计SHA256',source(s).auditSHA256||'未提供'],['适用范围',s.c.scope],...Object.keys(fields).map(key=>[fields[key][0]+'方法',method(s,key)]),['字段状态',Number.isFinite(value(s,field,face))?s.f.valueRole||s.c.fieldStateLabel||'已保存的限定基准数值':'未求解；不能用零替代'],['物理资格','研究验证；未授予真实箭体解资格']]);
 }
 function draw(){
  if(sectionEnabled){drawSection();return;}
  const canvas=$('curve'),s=state(),rect=canvas.getBoundingClientRect(),w=Math.max(300,rect.width),h=160,dpr=window.devicePixelRatio||1;canvas.width=w*dpr;canvas.height=h*dpr;const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);ctx.clearRect(0,0,w,h);ctx.fillStyle='#fff';ctx.fillRect(0,0,w,h);
  const [name,unit]=fields[field],limits=range(s,field),left=58,right=w-18,top=22,bottom=126;ctx.strokeStyle='#c6d2de';ctx.beginPath();ctx.moveTo(left,top);ctx.lineTo(left,bottom);ctx.lineTo(right,bottom);ctx.stroke();ctx.fillStyle='#506779';ctx.font='11px system-ui';ctx.fillText(name+' / '+unit,left,13);ctx.fillText('面元序号',right-55,151);ctx.fillText(fmt(limits[1]),2,top+4);ctx.fillText(fmt(limits[0]),2,bottom+4);
  let drawing=false,known=0;ctx.strokeStyle='#256686';ctx.lineWidth=1.4;ctx.beginPath();s.F.forEach((_,i)=>{const v=value(s,field,i);if(!Number.isFinite(v)){drawing=false;return;}known++;const x=left+(right-left)*i/Math.max(1,s.F.length-1),y=bottom-(bottom-top)*(v-limits[0])/(limits[1]-limits[0]);if(drawing)ctx.lineTo(x,y);else ctx.moveTo(x,y);drawing=true;});ctx.stroke();if(!known){ctx.fillStyle='#63768a';ctx.fillText('当前范围的此场量尚未求解',left+25,77);}
 }
 // 模块用途：复用固定版本的真实截面求交；先压缩未引用顶点，避免错误的截面范围/取值侧。
 function sectionMesh(s){
  const key=selected+'|'+frame+'|'+view;if(meshCache&&meshCache.key===key)return meshCache;
  const used=[...new Set(s.F.flat())].sort((a,b)=>a-b),ids=new Map(used.map((i,j)=>[i,j]));
  const V=used.map(i=>s.V[i]),F=s.F.map(f=>f.map(i=>ids.get(i))),xs=V.map(v=>v[0]);
  meshCache={key,V,F,min:Math.min(...xs),max:Math.max(...xs)};return meshCache;
 }
 function sectionField(s){
  const values=s.F.map((_,i)=>value(s,field,i));
  return {label:fields[field][0],unit:fields[field][1],values,status:values.map(v=>Number.isFinite(v)?1:0),
   methodNames:[method(s,field)],methodIndex:values.map(()=>0),scope:s.full?'原箭体物理场未求':s.c.scope,
   computedStateLabel:s.f.valueRole||s.c.fieldStateLabel||'限定基准数值；未授予工程区域资格'};
 }
 function updateSection(){
  const s=state(),mesh=sectionMesh(s);if(!Number.isFinite(sliceX))sliceX=(mesh.min+mesh.max)/2;
  sliceX=Math.max(mesh.min,Math.min(mesh.max,sliceX));
  section=window.RocketSection.geometry(mesh.V,mesh.F,sliceX);sectionSamples=window.RocketSection.sample(section,sectionField(s));
  $('sliceX').min=String(mesh.min);$('sliceX').max=String(mesh.max);$('sliceX').value=String(sliceX);
  $('slice').value=String(mesh.max>mesh.min?1000*(sliceX-mesh.min)/(mesh.max-mesh.min):0);$('slice').disabled=!(mesh.max>mesh.min);
  $('sectionSummary').textContent=fmt(sliceX)+' m · '+section.segments.length+' 条交线 / '+section.components.length+' 个轮廓 · '+sectionSamples.filter(s=>Number.isFinite(s.value)).length+' 项有效值';
  drawSection();sendSection();
 }
 function sendSection(){
  if(!sectionEnabled||!section||typeof $('viewer').contentWindow.postMessage!=='function')return;
  const key=revision+'|'+sliceX;if(key===lastOverlayKey)return;lastOverlayKey=key;
  $('viewer').contentWindow.postMessage({kind:'ale-section-overlay',caseIndex:selected,frameIndex:frame,revision,view,segments:section.segments.map(s=>[s.a,s.b])},location.origin);
 }
 // 模块用途：面元常值按实际周向角跨度绘制；0/360°拆段，未求和无效角度保持断开。
 function drawSection(){
  const s=state(),canvas=$('curve'),rect=canvas.getBoundingClientRect(),w=Math.max(300,rect.width),h=160,dpr=window.devicePixelRatio||1;
  canvas.width=w*dpr;canvas.height=h*dpr;const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);ctx.clearRect(0,0,w,h);ctx.fillStyle='#fff';ctx.fillRect(0,0,w,h);
  const limits=range(s,field),left=58,right=w-18,top=22,bottom=126;ctx.strokeStyle='#c6d2de';ctx.beginPath();ctx.moveTo(left,top);ctx.lineTo(left,bottom);ctx.lineTo(right,bottom);ctx.stroke();
  ctx.fillStyle='#506779';ctx.font='11px system-ui';ctx.fillText(fields[field][0]+' / '+fields[field][1],left,13);ctx.fillText('0°',left-5,148);ctx.fillText('180°',(left+right)/2-12,148);ctx.fillText('360°',right-23,148);ctx.fillText(fmt(limits[1]),2,top+4);ctx.fillText(fmt(limits[0]),2,bottom+4);
  plotSegments=[];const colors=['#256686','#876c48','#567c6d'];
  sectionSamples.forEach(segment=>{if(!Number.isFinite(segment.value))return;window.RocketSection.spans(segment).forEach(([a,b])=>{
   const x0=left+(right-left)*a/360,x1=left+(right-left)*b/360,y=bottom-(bottom-top)*(segment.value-limits[0])/(limits[1]-limits[0]);
   ctx.strokeStyle=colors[segment.component%colors.length];ctx.lineWidth=segment.face===face?2.8:1.5;ctx.beginPath();ctx.moveTo(x0,y);ctx.lineTo(x1,y);ctx.stroke();plotSegments.push({x0,x1,y,segment});
  });});
  if(!plotSegments.length){ctx.fillStyle='#63768a';ctx.fillText(section&&section.segments.length?'此截面没有可绘的有效角度/场量；几何交线仍保留':'此 X 平面没有非退化交线',left+15,77);}
 }
 function download(text,name){
  if(sectionEnabled){$('csvPreviewDetails').hidden=false;$('csvPreview').value=text.slice(0,65536);}
  const blob=new Blob([text],{type:'text/csv;charset=utf-8'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;a.click();URL.revokeObjectURL(url);
 }
 function caseDefaults(){const c=current();if(c.defaultView&&['all','wall','outer'].includes(c.defaultView)){view=c.defaultView;$('view').value=view;}if(c.defaultField&&fields[c.defaultField]){field=c.defaultField;$('field').value=field;}}
 data.cases.forEach((c,i)=>{const o=document.createElement('option');o.value=String(i);o.textContent=c.caption;$('caseChoice').append(o);});
 $('caseChoice').value=String(selected);caseDefaults();
 buttons.forEach((b,i)=>b.addEventListener('click',()=>{step=i;pickPosition=null;show();}));
 $('caseChoice').addEventListener('change',()=>{const i=Number($('caseChoice').value);if(!Number.isInteger(i)||!data.cases[i])return;selected=i;if(i!==1&&view==='rocket'){view='all';$('view').value=view;}frame=current().frames.length-1;face=0;pickPosition=null;sliceX=null;caseDefaults();show();});
 $('view').addEventListener('change',()=>{if(!['all','wall','outer','rocket'].includes($('view').value))return;view=$('view').value;if(view==='rocket'){selected=1;$('caseChoice').value='1';}face=0;pickPosition=null;sliceX=null;show();});
 $('field').addEventListener('change',()=>{if(!fields[$('field').value])return;field=$('field').value;show();});
 $('frame').addEventListener('input',()=>{const i=Number($('frame').value);if(!Number.isInteger(i)||!current().frames[i])return;frame=i;pickPosition=null;show();});
 $('face').addEventListener('input',()=>{const i=Number($('face').value);if(!Number.isInteger(i)||i<0||i>=state().F.length)return;face=i;pickPosition=null;probe();if(sectionEnabled){sliceX=center(state(),face)[0];updateSection();}});
 window.addEventListener('message',event=>{const m=event.data;if(event.origin!==location.origin||event.source!==$('viewer').contentWindow||!m||m.kind!=='ale-view-event'||m.caseIndex!==selected||m.frameIndex!==frame||!m.payload||m.payload.revision!==revision)return;const p=m.payload;
  if(p.kind==='geometry-rendered'){$('renderState').textContent=p.webglError===0&&p.visiblePixelSamples>0?'三维网格已显示 · 可点选面元':'三维显示需要检查';if(sectionEnabled&&!Object.prototype.hasOwnProperty.call(p,'sectionSegments')){lastOverlayKey='';sendSection();}return;}
  if(step!==4||p.kind!=='geometry-face-probe'||!Number.isInteger(p.faceIndex)||p.faceIndex<0||p.faceIndex>=state().F.length||!Array.isArray(p.barycentric)||p.barycentric.length!==3||p.barycentric.some(x=>!Number.isFinite(x)||x< -1e-9||x>1+1e-9)||Math.abs(p.barycentric.reduce((a,b)=>a+b,0)-1)>1e-9)return;
  face=p.faceIndex;const s=state();pickPosition=[0,1,2].map(j=>p.barycentric.reduce((v,t,i)=>v+t*s.V[s.F[face][i]][j],0));$('face').value=String(face);probe();if(sectionEnabled){sliceX=pickPosition[0];updateSection();}
 });
 $('csv').addEventListener('click',()=>{
  const s=state(),quote=v=>'"'+(v==null?'':String(v)).replace(/"/g,'""')+'"';
  const header=['view_face_id','original_face_id','x_m','y_m','z_m','saved_time_s','pressure_Pa','fluid_temperature_K','speed_m_s','Riemann_wall_pressure_Pa','heat_flux_W_m2','solved_wall_temperature_K','edge_temperature_K','recovery_temperature_K','input_wall_temperature_K','state_role','pressure_method','temperature_method','heat_flux_status','production_qualified'];
  if(sectionEnabled)header.push('case_id','source_evidence_sha256','source_audit_sha256','value_role','wall_pressure_method','edge_temperature_method','recovery_temperature_method','input_wall_temperature_method','solved_wall_temperature_method','heat_flux_method');
  const lines=[header.map(quote).join(',')];s.F.forEach((_,i)=>{
   const index=s.map[i],id=s.full?index+1:s.c.sourceFaceIDs[index];
   const values=[i+1,id>0?id:null,...center(s,i),s.f.timeS,...['pressurePa','temperatureK','speedMS','wallPressurePa','heatFluxWm2','wallTemperatureK'].map(key=>value(s,key,i)),...['edgeTemperatureK','recoveryTemperatureK','inputWallTemperatureK'].map(key=>value(s,key,i)),s.full?'原完整箭体几何；真实物理场未求':s.c.stateRole,method(s,'pressurePa'),method(s,'temperatureK'),'unknown',false];
   if(sectionEnabled)values.push(s.c.id,source(s).evidenceSHA256,source(s).auditSHA256,s.full?'原箭体未求':s.f.valueRole||s.c.fieldStateLabel,...['wallPressurePa','edgeTemperatureK','recoveryTemperatureK','inputWallTemperatureK','wallTemperatureK','heatFluxWm2'].map(key=>method(s,key)));
   lines.push(values.map(quote).join(','));
  });download('\uFEFF'+lines.join('\n')+'\n','ale_'+current().id+'_'+view+'_frame'+frame+'.csv');
 });
 if(sectionEnabled){
  $('sectionControls').hidden=false;$('curveTitle').textContent='真实截面周向分布';
  $('curveNote').textContent='点击三维面元选择真实 X，也可输入或拖动滑块。沿交线保留面元常值；+Y 为 0°，+Z 为 90°。开放/分支轮廓单列，未解值断开。';
  $('sliceX').addEventListener('change',()=>{const x=Number($('sliceX').value);if(!Number.isFinite(x))return;sliceX=x;updateSection();});
  $('slice').addEventListener('input',()=>{const t=Number($('slice').value);if(!Number.isFinite(t)||t<0||t>1000)return;const mesh=sectionMesh(state());sliceX=mesh.min+t/1000*(mesh.max-mesh.min);updateSection();});
  $('sectionCsv').addEventListener('click',()=>{const s=state();const text=window.RocketSection.csv(section,sectionField(s),'t='+s.f.timeS+' s; '+s.c.id+'; evidenceSHA256='+source(s).evidenceSHA256,s.c.conditionLabel||s.c.scope);download(text,'ale_'+s.c.id+'_section_'+sliceX+'.csv');});
  $('curve').addEventListener('click',event=>{const rect=$('curve').getBoundingClientRect(),x=event.clientX-rect.left,y=event.clientY-rect.top;let chosen=null,best=144;
   plotSegments.forEach(p=>{const dx=x-Math.max(p.x0,Math.min(p.x1,x)),dy=y-p.y,d=dx*dx+dy*dy;if(d<best){best=d;chosen=p.segment;}});
   if(!chosen)return;face=chosen.face;pickPosition=chosen.a.map((v,j)=>(v+chosen.b[j])/2);$('face').value=String(face);probe();drawSection();
  });
 }
 if(typeof ResizeObserver==='function')new ResizeObserver(draw).observe($('curve'));show();
})();
