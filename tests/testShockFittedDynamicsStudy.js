/* 模块用途：实际研究数据的流程、面元探针、真实周向CSV、截面叠加和消息隔离核对。 */
'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const folder=path.resolve(process.argv[2]),target=path.resolve(process.argv[3]),data=JSON.parse(fs.readFileSync(path.join(folder,'ale-data.json'),'utf8'));
const sectionCore=require(path.join(folder,'ale-section-core.js'));
assert.equal(data.kind,'shock-fitted-dynamics-study-v102');assert.equal(data.productionQualified,false);assert.equal(data.newProductionSolvedAreaM2,0);
const canvasContext=Object.fromEntries(['scale','clearRect','fillText','fillRect','beginPath','moveTo','lineTo','stroke'].map(name=>[name,(...args)=>assert.ok(args.filter(x=>typeof x==='number').every(Number.isFinite))]));
class E{
 constructor(id){this.id=id;this.events={};this.attrs={};this.children=[];this.value='';this.hidden=false;}
 append(...items){this.children.push(...items);}replaceChildren(...items){this.children=items;}
 addEventListener(k,fn){this.events[k]=fn;}setAttribute(k,v){this.attrs[k]=String(v);}
 getBoundingClientRect(){return{left:0,top:0,width:800,height:160};}getContext(){return canvasContext;}
}
const ids=Object.fromEntries(['caseChoice','title','scope','gasScope','details','frame','time','face','viewer','renderState','role','probe','methods','curve','view','field','csv','sectionControls','sliceX','slice','sectionCsv','sectionSummary','curveTitle','curveNote','csvPreviewDetails','csvPreview'].map(id=>[id,new E(id)]));
const overlayMessages=[],listeners={},buttons=Array.from({length:5},()=>new E('step'));
ids.viewer.contentWindow={postMessage:m=>overlayMessages.push(m)};
let blob=null,downloads=0;
const document={getElementById:id=>{assert.ok(ids[id],id);return ids[id];},querySelectorAll:()=>buttons,createElement:name=>{const e=new E(name);if(name==='a')e.click=()=>downloads++;return e;}};
const location={origin:'http://127.0.0.1:8782'},window={ALE_STUDY:data,RocketSection:sectionCore,devicePixelRatio:1,addEventListener:(k,fn)=>listeners[k]=fn};
vm.runInNewContext(fs.readFileSync(path.join(folder,'ale-study.js'),'utf8'),{window,document,location,Blob,URLSearchParams,URL:{createObjectURL:b=>{blob=b;return'blob:test';},revokeObjectURL:()=>{}}});
function select(id,value,event='change'){ids[id].value=String(value);ids[id].events[event]();}
function row(id,key){const i=ids[id].children.findIndex(x=>x.textContent===key);assert.ok(i>=0,key);return ids[id].children[i+1].textContent;}
function query(){return new URL(ids.viewer.attrs.src,location.origin).searchParams;}
function mapFor(c,view){return c.faces.map((_,i)=>i).filter(i=>view!=='wall'&&view!=='outer'||c.boundaryKind[i]===(view==='wall'?2:3));}
function message(k,f,i,payload={}){return{origin:location.origin,source:ids.viewer.contentWindow,data:{kind:'ale-view-event',caseIndex:k,frameIndex:f,payload:{kind:'geometry-face-probe',revision:Number(query().get('revision')),faceIndex:i,barycentric:[1/3,1/3,1/3],...payload}}};}
const fmt=v=>Number.isFinite(v)?Number(v.toPrecision(8)).toString():'未求解';
const parse=line=>[...line.matchAll(/"((?:[^"]|"")*)"(?:,|$)/g)].map(x=>x[1].replace(/""/g,'"'));
let workflow=0,probes=0,sections=0,csvRows=0,ignored=0,geometryChecks=0;
const exported=[];
(async()=>{
 assert.equal(ids.caseChoice.value,String(data.defaultCaseIndex));assert.equal(ids.view.value,'wall');assert.equal(ids.field.value,'wallPressurePa');
 assert.equal(ids.curveTitle.textContent,'真实截面周向分布');
 for(let k=0;k<data.cases.length;k++){
  const c=data.cases[k];select('caseChoice',k);
  for(let step=0;step<5;step++){buttons[step].events.click();assert.equal(query().get('layer'),step===4?'1':'0');workflow++;}
  const views=['all','wall',...(c.boundaryKind.includes(3)?['outer']:[])];
  for(const view of views){select('view',view);const map=mapFor(c,view);
   for(let f=0;f<c.frames.length;f++){
    select('frame',f,'input');
    // 所有局部边界面均点选；真实整箭灰色几何另以抽样检查，避免重复二次开销。
    for(let i=0;i<map.length;i++){
     listeners.message(message(k,f,i));assert.equal(ids.face.value,String(i));
     assert.equal(row('probe','流体压力 / Pa'),fmt(c.frames[f].pressurePa[map[i]]));
     assert.equal(row('probe','Riemann壁压 / Pa'),fmt(c.frames[f].wallPressurePa[map[i]]));
     assert.equal(row('probe','热流 / W/m²'),'未求解');assert.equal(row('probe','给定壁温 / K'),'未设置');
     assert.equal(row('methods','数据SHA256'),(c.source||data.source).evidenceSHA256);probes++;
    }
    for(const key of ['pressurePa','wallPressurePa','temperatureK','speedMS','heatFluxWm2','wallTemperatureK','edgeTemperatureK','recoveryTemperatureK','inputWallTemperatureK'])select('field',key);
    select('field',k>=2&&k<=4?'wallPressurePa':'pressurePa');
    const used=[...new Set(map.flatMap(i=>c.faces[i]))],xs=used.map(i=>c.frames[f].position[i][0]);
    const x=(Math.min(...xs)+Math.max(...xs))/2;select('sliceX',x);
    const last=overlayMessages.at(-1);assert.equal(last.kind,'ale-section-overlay');assert.equal(last.view,view);assert.equal(last.caseIndex,k);assert.equal(last.frameIndex,f);
    last.segments.forEach(segment=>segment.forEach(p=>assert.ok(Math.abs(p[0]-x)<1e-12)));
    ids.sectionCsv.events.click();const text=(await blob.text()).replace(/^\uFEFF/,''),lines=text.trim().split(/\r?\n/),header=parse(lines[0]);assert.equal(header.length,18);
    const frame=c.frames[f],selectedField=ids.field.value;
    lines.slice(1).forEach(line=>{const r=parse(line),face=Number(r[2])-1;assert.equal(r.length,18);assert.ok(face>=0&&face<map.length);assert.equal(Number(r[0]),x);const v=Array.isArray(frame[selectedField])?frame[selectedField][map[face]]:null;assert.equal(r[10],v===null?'':String(v));assert.equal(r[12],c.methods[selectedField]);assert.ok(r[14].includes((c.source||data.source).evidenceSHA256));csvRows++;});
    const name=c.id+'_'+view+'_frame'+f+'.csv';exported.push({name,caseIndex:k,frameIndex:f,view,field:selectedField,x,rows:lines.length-1,text});sections++;
    ids.csv.events.click();const all=(await blob.text()).replace(/^\uFEFF/,'').trim().split(/\r?\n/),columns=parse(all[0]);assert.equal(columns.length,30);assert.equal(all.length,map.length+1);
    all.slice(1).forEach((line,i)=>{const values=parse(line);assert.equal(values.length,columns.length);const r=Object.fromEntries(columns.map((name,j)=>[name,values[j]]));assert.equal(r.source_evidence_sha256,(c.source||data.source).evidenceSHA256);assert.equal(r.case_id,c.id);assert.equal(r.production_qualified,'false');assert.equal(r.heat_flux_W_m2,'');assert.equal(r.edge_temperature_K,'');assert.equal(r.recovery_temperature_K,'');assert.equal(r.input_wall_temperature_K,'');assert.equal(r.solved_wall_temperature_K,'');csvRows++;});
   }
  }
 }
 // 模块二：原整箭值保持未知；选择没有激波边界的波管时范围回到实际存在的边界。
 select('view','rocket');assert.equal(query().get('case'),'1');
 for(const i of [0,17,203,data.originalGeometry.faces.length-1]){listeners.message(message(1,1,i));assert.equal(row('probe','流体压力 / Pa'),'未求解');probes++;}
 select('caseChoice',5);select('view','outer');assert.equal(query().get('view'),'all');
 // 模块三：空截面叠加回执不得无限重发，非法消息不得改变探针。
 select('caseChoice',4);select('view','wall');select('field','wallPressurePa');select('sliceX',.05);
 assert.equal(overlayMessages.at(-1).segments.length,0);
 const before=overlayMessages.length;listeners.message(message(4,1,0,{kind:'geometry-rendered',webglError:0,visiblePixelSamples:10,sectionSegments:0}));assert.equal(overlayMessages.length,before);
 const valid=message(4,1,0),bad=[{...valid,origin:'https://invalid.test'},{...valid,source:{}},{...valid,data:{...valid.data,caseIndex:0}},{...valid,data:{...valid.data,frameIndex:0}}];
 for(const p of [{revision:0},{faceIndex:-1},{faceIndex:1e9},{barycentric:[1,1,1]},{barycentric:[NaN,0,1]},{barycentric:[-1,1,1]},{kind:'surface-flow-probe'}])bad.push({...valid,data:{...valid.data,payload:{...valid.data.payload,...p}}});
 for(const event of bad){const p=JSON.stringify(ids.probe.children);listeners.message(event);assert.equal(JSON.stringify(ids.probe.children),p);ignored++;}
 // 模块四：实际几何适配器压缩未引用顶点，接收截面时保留相机，拒绝非父窗口来源。
 for(let k=0;k<data.cases.length;k++)for(const view of ['all','wall',...(data.cases[k].boundaryKind.includes(3)?['outer']:[])]){
  const c=data.cases[k],f=c.frames[1],map=mapFor(c,view),geometryListeners={};let payload=null;
  const parent={postMessage:()=>{}},geoWindow={ALE_STUDY:data,setup:bridge=>bridge.addEventListener('DataChanged',()=>{payload=bridge.Data;}),addEventListener:(name,fn)=>geometryListeners[name]=fn};
  vm.runInNewContext(fs.readFileSync(path.join(folder,'geometry/ale-study-geometry.js'),'utf8'),{window:geoWindow,location:{origin:location.origin,search:'?case='+k+'&frame=1&revision=2&view='+view+'&field=pressurePa&layer=1&wire=0'},URLSearchParams,parent,document:{querySelector:()=>new E('note')}});
  assert.equal(new Set(payload.faces.flat()).size,payload.vertices.length);assert.equal(payload.faces.length,map.length);
  payload.faces.forEach((face,i)=>face.forEach((id,j)=>assert.deepEqual(Array.from(payload.vertices[id]),f.position[c.faces[map[i]][j]])));
  const empty={kind:'ale-section-overlay',caseIndex:k,frameIndex:1,revision:2,view,segments:[]};geometryListeners.message({source:parent,origin:location.origin,data:empty});assert.equal(payload.resetCamera,false);assert.equal(payload.sectionSegments.length,0);
  const saved=JSON.stringify(payload);geometryListeners.message({source:{},origin:location.origin,data:{...empty,segments:[[[0,0,0],[1,1,1]]]}});assert.equal(JSON.stringify(payload),saved);geometryChecks++;
 }
 fs.mkdirSync(target,{recursive:true});
 const records=exported.map(({text,...record})=>{const p=path.join(target,record.name);fs.writeFileSync(p,text);return{...record,bytes:Buffer.byteLength(text),sha256:crypto.createHash('sha256').update(text).digest('hex')};});
 const result={status:'passed-actual-dynamic-study-probes-sections-CSV-and-camera-preserving-overlays',workflowEvents:workflow,probedFaces:probes,sectionExports:sections,csvRows,downloads,rejectedMessages:ignored,geometrySubsetsVerified:geometryChecks,emptyOverlayAcknowledgementNoLoop:true,dataSHA256:crypto.createHash('sha256').update(fs.readFileSync(path.join(folder,'ale-data.json'))).digest('hex'),sectionFiles:records,productionQualified:false,newProductionSolvedAreaM2:0};
 fs.writeFileSync(path.join(target,'node_audit.json'),JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify({...result,sectionFiles:records.length}));
})().catch(e=>{console.error(e);process.exitCode=1;});
