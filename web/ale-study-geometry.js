/* 模块用途：向通用科研几何管线传递实际逐面值、固定色标和可追溯拾取事件。 */
'use strict';
(function(){
 const data=window.ALE_STUDY,q=new URLSearchParams(location.search),caseIndex=Number(q.get('case')),frameIndex=Number(q.get('frame')),revision=Number(q.get('revision')),view=q.get('view'),field=q.get('field');
 if(!data||!['shock-fitted-kernel-study-v101','shock-fitted-dynamics-study-v102','shock-fitted-interface-study-v103'].includes(data.kind)||!Number.isSafeInteger(caseIndex)||!data.cases[caseIndex]||!Number.isSafeInteger(frameIndex)||!data.cases[caseIndex].frames[frameIndex]||!Number.isSafeInteger(revision)||revision<1)throw Error('研究来源或帧身份无效。');
 const study=data.cases[caseIndex],frame=study.frames[frameIndex],pressureRange=study.pressureRangePa||data.pressureRangePa,fields={pressurePa:['流体压力','Pa',pressureRange],temperatureK:['流体温度（不是壁温）','K',data.temperatureRangeK],speedMS:['速度幅值','m/s',[0,2000]],wallPressurePa:['Riemann壁压','Pa',pressureRange],heatFluxWm2:['热流（未求）','W/m²',[0,1000000]],wallTemperatureK:['求解壁温（未求）','K',[80,2000]],edgeTemperatureK:['外缘Te（未求）','K',[80,2000]],recoveryTemperatureK:['恢复Tr（未求）','K',[80,2000]],inputWallTemperatureK:['给定壁温（未设置）','K',[80,2000]]};
 if(!fields[field]||!['all','wall','outer','rocket'].includes(view))throw Error('场量或范围无效。');
 const indices=study.faces.map((_,i)=>i).filter(i=>view==='wall'?study.boundaryKind[i]===2:view==='outer'?[3,5].includes(study.boundaryKind[i]):true);
 let vertices=frame.position,faces=indices.map(i=>study.faces[i]),values=indices.map(i=>q.get('layer')==='1'&&Array.isArray(frame[field])?frame[field][i]:null);
 if(view==='rocket'){vertices=data.originalGeometry.vertices;faces=data.originalGeometry.faces;values=faces.map(()=>null);}
 // 只向绘图器传递当前面集合真正使用的顶点；壁面范围不能被未显示的外壳节点扩大。
 const used=[...new Set(faces.flat())].sort((a,b)=>a-b),indexMap=new Map(used.map((id,i)=>[id,i]));
 vertices=used.map(id=>vertices[id]);faces=faces.map(f=>f.map(id=>indexMap.get(id)));
 const listeners=[];let value=null;const bridge={addEventListener:(name,fn)=>{if(name==='DataChanged')listeners.push(fn);}};
 Object.defineProperty(bridge,'Data',{get:()=>value,set:v=>{value=v;if(v&&['geometry-rendered','geometry-face-probe','geometry-error'].includes(v.kind))parent.postMessage({kind:'ale-view-event',caseIndex,frameIndex,payload:v},location.origin);}});
 window.setup(bridge);
 const scene={kind:'aerothermal-geometry',schemaVersion:1,revision,vertices,faces,caption:study.caption,showWire:q.get('wire')==='1',resetCamera:true,initialCameraAzimuthRad:-1.05,initialCameraElevationRad:.4,coordinateNote:study.coordinateNote||'SI / m · 已保存研究 · 原箭体物理场未求',faceField:{values,range:fields[field][2],name:fields[field][0],unit:fields[field][1],productionQualified:false}};
 function update(){value=scene;listeners.forEach(fn=>fn());const note=document.querySelector('footer span:last-child');if(note)note.textContent=view==='rocket'||study.id==='real-nose'?'真实箭体场量未求 · 灰色不代表零':study.geometryNote||'移动活塞解析基准 · 不是箭体解';}
 update();
 // 模块用途：截面改变只更新空间交线，不重置相机；外部窗口、旧帧和非法坐标不接受。
 if(typeof window.addEventListener==='function')window.addEventListener('message',event=>{
  const m=event.data;if(event.source!==parent||event.origin!==location.origin||!m||m.kind!=='ale-section-overlay'||m.caseIndex!==caseIndex||m.frameIndex!==frameIndex||m.revision!==revision||m.view!==view||!Array.isArray(m.segments)||m.segments.length>100000||m.segments.some(s=>!Array.isArray(s)||s.length!==2||s.some(p=>!Array.isArray(p)||p.length!==3||p.some(x=>!Number.isFinite(x)))))return;
  scene.sectionSegments=m.segments.map(s=>s.map(p=>p.slice()));scene.resetCamera=false;update();
 });
})();
