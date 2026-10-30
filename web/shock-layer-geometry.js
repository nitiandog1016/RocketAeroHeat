/* 模块用途：复用几何WebGL画线/拾取原语，明确叠加的是法向列而非表面流线。 */
'use strict';
(function(){
 const study=window.SHOCK_LAYER_STUDY,q=new URLSearchParams(location.search),index=Number(q.get('case')),revision=Number(q.get('revision'));
 if(!study||study.kind!=='shock-layer-study-v1'||study.productionQualified!==false||!Number.isInteger(index)||!study.cases[index]||!Number.isSafeInteger(revision)||revision<1)throw Error('激波层几何请求身份无效。');
 const c=study.cases[index],focus=q.get('view')==='nose',geometry=focus?study.geometry.focus:study.geometry;
 const listeners=[];let value=null;
 const bridge={addEventListener:(name,fn)=>{if(name==='DataChanged')listeners.push(fn);}};
 Object.defineProperty(bridge,'Data',{get:()=>value,set:v=>{value=v;if(v&&['geometry-rendered','surface-flow-probe','geometry-error'].includes(v.kind))parent.postMessage({kind:'shock-layer-view-event',caseIndex:index,payload:v},location.origin);}});
 window.setup(bridge);
 value={kind:'aerothermal-geometry',schemaVersion:1,revision,vertices:geometry.vertices,faces:geometry.faces,
  caption:focus?'原箭体鼻部开放片 · 层内法向列':'原完整箭体 · 层内法向列',resetCamera:true,showWire:q.get('wire')==='1',
  coordinateNote:'SI / m · 原X轴与鼻尖方向保持 · 层内候选值未赋给真实壁面',
  initialCameraAzimuthRad:Number.isFinite(study.initialCameraAzimuthRad)?study.initialCameraAzimuthRad:-.35,
  initialCameraElevationRad:Number.isFinite(study.initialCameraElevationRad)?study.initialCameraElevationRad:.3};
 if(q.get('layer')==='1')value.surfaceFlowLines={kind:'pressure-surface-flow-research',productionQualified:false,curves:c.columns.map(col=>({points:col.position,normals:col.pickNormals}))};
 listeners.forEach(fn=>fn());
 const note=document.querySelector('footer span:last-child');if(note)note.textContent=q.get('layer')==='1'?'激波内法向坐标列（不是流线） · 未通过体形匹配 · 热流未求':'原几何预览 · 无新增物理解面积';
})();
