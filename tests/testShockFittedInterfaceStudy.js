/* 模块用途：操作实际导出页面，核对当前帧壁压、五量诊断、类型5外壳及CSV时间身份。 */
'use strict';
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert/strict'),crypto=require('crypto');
const folder=path.resolve(process.argv[2]),target=path.resolve(process.argv[3]);
assert.ok(!fs.existsSync(target),'验证目录须为新目录');
const data=JSON.parse(fs.readFileSync(path.join(folder,'ale-data.json'),'utf8'));
assert.equal(data.kind,'shock-fitted-interface-study-v103');assert.equal(data.productionQualified,false);
const sectionCore=require(path.join(folder,'ale-section-core.js'));
const context=new Proxy({}, {get:(_,k)=>['strokeStyle','fillStyle','font','lineWidth'].includes(k)?undefined:(...args)=>assert.ok(args.filter(x=>typeof x==='number').every(Number.isFinite))});
class Element{
 constructor(id){this.id=id;this.events={};this.attrs={};this.children=[];this.value='';this.hidden=false;}
 append(...items){this.children.push(...items);}replaceChildren(...items){this.children=items;}
 addEventListener(k,fn){this.events[k]=fn;}setAttribute(k,v){this.attrs[k]=String(v);}
 getBoundingClientRect(){return{left:0,top:0,width:800,height:160};}getContext(){return context;}
}
const names=['caseChoice','title','scope','gasScope','details','frame','time','face','viewer','renderState','role','probe','methods','curve','view','field','csv','sectionControls','sliceX','slice','sectionCsv','sectionSummary','curveTitle','curveNote','csvPreviewDetails','csvPreview','interfaceDiagnostics','diagnosticTime','interfaceProbe'];
const ids=Object.fromEntries(names.map(k=>[k,new Element(k)])),buttons=Array.from({length:5},()=>new Element('step'));
const overlays=[],listeners={};ids.viewer.contentWindow={postMessage:m=>overlays.push(m)};
let blob,downloads=0;
const document={getElementById:id=>{assert.ok(ids[id],id);return ids[id];},querySelectorAll:()=>buttons,createElement:name=>{const e=new Element(name);if(name==='a')e.click=()=>downloads++;return e;}};
const location={origin:'http://127.0.0.1:8783'},window={ALE_STUDY:data,RocketSection:sectionCore,devicePixelRatio:1,addEventListener:(k,fn)=>listeners[k]=fn};
vm.runInNewContext(fs.readFileSync(path.join(folder,'ale-study.js'),'utf8'),{window,document,location,Blob,URLSearchParams,URL:{createObjectURL:b=>{blob=b;return'blob:test';},revokeObjectURL:()=>{}}});
const fmt=v=>Number.isFinite(v)?Number(v.toPrecision(8)).toString():'未求解';
const parse=line=>[...line.matchAll(/"((?:[^"]|"")*)"(?:,|$)/g)].map(x=>x[1].replace(/""/g,'"'));
function select(id,value,event='change'){ids[id].value=String(value);ids[id].events[event]();}
function row(id,key){const i=ids[id].children.findIndex(x=>x.textContent===key);assert.ok(i>=0,key);return ids[id].children[i+1].textContent;}
function query(){return new URL(ids.viewer.attrs.src,location.origin).searchParams;}
function mapFor(c,view){return c.faces.map((_,i)=>i).filter(i=>view==='wall'?c.boundaryKind[i]===2:view==='outer'?[3,5].includes(c.boundaryKind[i]):true);}
function message(k,f,i,extra={}){return{origin:location.origin,source:ids.viewer.contentWindow,data:{kind:'ale-view-event',caseIndex:k,frameIndex:f,payload:{kind:'geometry-face-probe',revision:Number(query().get('revision')),faceIndex:i,barycentric:[1/3,1/3,1/3],...extra}}};}
const exported=[];let probes=0,workflow=0,csvRows=0,geometryChecks=0;
(async()=>{
 assert.equal(ids.caseChoice.value,String(data.defaultCaseIndex));
 for(let k=0;k<data.cases.length;k++){
  const c=data.cases[k];select('caseChoice',k);
  for(let stage=0;stage<5;stage++){buttons[stage].events.click();assert.equal(query().get('layer'),stage===4?'1':'0');workflow++;}
  for(const view of ['all','wall','outer'])for(let f=0;f<2;f++){
   select('view',view);select('frame',f,'input');const frame=c.frames[f],map=mapFor(c,view);
   assert.equal(query().get('view'),view);assert.equal(frame.wallPressureTimeS,frame.timeS);
   for(let i=0;i<map.length;i++){
    listeners.message(message(k,f,i));const index=map[i];
    assert.equal(row('probe','流体压力 / Pa'),fmt(frame.pressurePa[index]));
    assert.equal(row('probe','Riemann壁压 / Pa'),fmt(frame.wallPressurePa[index]));
    assert.equal(row('probe','热流 / W/m²'),'未求解');assert.equal(row('probe','给定壁温 / K'),'未设置');
    assert.equal(row('methods','数据SHA256'),c.source.evidenceSHA256);
    if(frame.diagnostics.available){
     assert.ok(frame.diagnostics.inputTimeS<frame.timeS);
     assert.equal(row('interfaceProbe','原单元压力相容差 / 1'),fmt(frame.diagnostics.cellPressureMismatch[index]));
     for(const [j,label] of ['RH质量通量缺陷 / 1','RH动量X通量缺陷 / 1','RH动量Y通量缺陷 / 1','RH动量Z通量缺陷 / 1','RH能量通量缺陷 / 1'].entries())assert.equal(row('interfaceProbe',label),fmt(frame.diagnostics.RHFluxDefect5[index][j]));
     assert.ok(ids.diagnosticTime.textContent.includes(fmt(frame.diagnostics.inputTimeS)));
    }else assert.ok(ids.diagnosticTime.textContent.includes('没有末步'));
    if(c.boundaryKind[index]===5)assert.equal(row('probe','边界身份'),'实际移动面Riemann候选激波');
    probes++;
   }
   for(const key of Object.keys(c.methods))select('field',key);
   const field=view==='wall'?'wallPressurePa':'pressurePa';select('field',field);
   const used=[...new Set(map.flatMap(i=>c.faces[i]))],xs=used.map(i=>frame.position[i][0]);
   const x=(Math.min(...xs)+Math.max(...xs))/2;select('sliceX',x);
   const overlay=overlays.at(-1);assert.equal(overlay.view,view);assert.equal(overlay.frameIndex,f);
   overlay.segments.forEach(s=>s.forEach(p=>assert.ok(Math.abs(p[0]-x)<1e-12)));
   ids.sectionCsv.events.click();const text=(await blob.text()).replace(/^\uFEFF/,''),lines=text.trim().split(/\r?\n/);
   lines.slice(1).forEach(line=>{const r=parse(line),face=map[Number(r[2])-1];assert.equal(r.length,18);assert.equal(r[10],frame[field][face]===null?'':String(frame[field][face]));assert.equal(r[12],c.methods[field]);csvRows++;});
   exported.push({name:c.id+'_'+view+'_frame'+f+'.csv',caseIndex:k,frameIndex:f,view,field,x,rows:lines.length-1,text});
   ids.csv.events.click();const all=(await blob.text()).replace(/^\uFEFF/,'').trim().split(/\r?\n/),header=parse(all[0]);
   assert.equal(header.length,42);assert.equal(all.length,map.length+1);
   all.slice(1).forEach((line,i)=>{const v=parse(line),r=Object.fromEntries(header.map((h,j)=>[h,v[j]]));assert.equal(v.length,42);assert.equal(Number(r.current_frame_time_s),frame.timeS);assert.equal(Number(r.wall_pressure_time_s),frame.timeS);assert.equal(r.production_qualified,'false');assert.equal(r.heat_flux_W_m2,'');assert.equal(r.edge_temperature_K,'');assert.equal(r.recovery_temperature_K,'');assert.equal(r.solved_wall_temperature_K,'');assert.equal(r.input_wall_temperature_K,'');assert.equal(r.source_evidence_sha256,c.source.evidenceSHA256);if(frame.diagnostics.available)assert.equal(Number(r.diagnostic_input_time_s),frame.diagnostics.inputTimeS);else assert.equal(r.diagnostic_input_time_s,'');csvRows++;});
  }
 }
 // 模块二：完整箭体只显示未知，物理时刻留空，来源绑定原几何而非球鼻候选。
 select('view','rocket');const originalIndex=data.originalCaseIndex;
 assert.equal(query().get('case'),String(originalIndex));assert.ok(ids.time.textContent.includes('未提供'));
 for(const i of [0,17,203,data.originalGeometry.faces.length-1]){listeners.message(message(originalIndex,1,i));assert.equal(row('probe','流体压力 / Pa'),'未求解');assert.equal(row('methods','数据SHA256'),data.originalGeometrySource.evidenceSHA256);probes++;}
 assert.equal(ids.interfaceDiagnostics.hidden,true);
 ids.csv.events.click();const all=(await blob.text()).replace(/^\uFEFF/,'').trim().split(/\r?\n/),header=parse(all[0]);
 all.slice(1).forEach(line=>{const values=parse(line),r=Object.fromEntries(header.map((name,j)=>[name,values[j]]));assert.equal(r.saved_time_s,'');assert.equal(r.current_frame_time_s,'');assert.equal(r.wall_pressure_time_s,'');assert.equal(r.diagnostic_input_time_s,'');assert.equal(r.source_evidence_sha256,data.originalGeometrySource.evidenceSHA256);csvRows++;});
 // 模块三：源窗口、帧、版本和重心均须匹配；非法事件不能改变真实探针。
 select('caseChoice',data.defaultCaseIndex);select('view','outer');select('frame',1,'input');
 const valid=message(data.defaultCaseIndex,1,0),bad=[{...valid,origin:'https://invalid.test'},{...valid,source:{}},{...valid,data:{...valid.data,caseIndex:99}},{...valid,data:{...valid.data,frameIndex:0}}];
 for(const p of [{revision:0},{faceIndex:-1},{faceIndex:1e9},{barycentric:[1,1,1]},{barycentric:[NaN,0,1]},{barycentric:[-1,1,1]},{kind:'wrong'}])bad.push({...valid,data:{...valid.data,payload:{...valid.data.payload,...p}}});
 for(const event of bad){const before=JSON.stringify(ids.probe.children);listeners.message(event);assert.equal(JSON.stringify(ids.probe.children),before);}
 // 模块四：真实绘图适配器须接受类型5外壳、固定色标、准确逐面值和保相机交线。
 for(let k=0;k<data.cases.length;k++)for(const view of ['all','wall','outer']){
  const c=data.cases[k],f=c.frames[1],map=mapFor(c,view),events={};let payload;
  const parent={postMessage:()=>{}},w={ALE_STUDY:data,setup:bridge=>bridge.addEventListener('DataChanged',()=>{payload=bridge.Data;}),addEventListener:(kind,fn)=>events[kind]=fn};
  vm.runInNewContext(fs.readFileSync(path.join(folder,'geometry/ale-study-geometry.js'),'utf8'),{window:w,location:{origin:location.origin,search:'?case='+k+'&frame=1&revision=2&view='+view+'&field=wallPressurePa&layer=1&wire=0'},URLSearchParams,parent,document:{querySelector:()=>new Element('note')}});
  assert.equal(payload.faces.length,map.length);assert.equal(new Set(payload.faces.flat()).size,payload.vertices.length);
  assert.deepEqual(Array.from(payload.faceField.range),[0,45000]);
  payload.faces.forEach((ids,i)=>{ids.forEach((id,j)=>assert.deepEqual(Array.from(payload.vertices[id]),f.position[c.faces[map[i]][j]]));assert.equal(payload.faceField.values[i],f.wallPressurePa[map[i]]);});
  events.message({source:parent,origin:location.origin,data:{kind:'ale-section-overlay',caseIndex:k,frameIndex:1,revision:2,view,segments:[]}});assert.equal(payload.resetCamera,false);geometryChecks++;
 }
 fs.mkdirSync(target,{recursive:true});const sectionFiles=exported.map(({text,...r})=>{fs.writeFileSync(path.join(target,r.name),text);return{...r,bytes:Buffer.byteLength(text),sha256:crypto.createHash('sha256').update(text).digest('hex')};});
 const result={status:'passed-current-frame-wall-five-component-diagnostics-type5-and-CSV-times',workflowEvents:workflow,probedFaces:probes,csvRows,sectionExports:sectionFiles.length,downloads,rejectedMessages:bad.length,geometrySubsetsVerified:geometryChecks,sectionFiles,productionQualified:false,newProductionSolvedAreaM2:0};
 fs.writeFileSync(path.join(target,'node_audit.json'),JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify({...result,sectionFiles:sectionFiles.length}));
})().catch(e=>{console.error(e);process.exitCode=1;});
