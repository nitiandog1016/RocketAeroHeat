/* 模块用途：实际导出数据的五步骤、全部面元、未知身份、CSV和同窗消息拒绝验证。 */
'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const folder=path.resolve(process.argv[2]),file=path.join(folder,'ale-data.json'),raw=fs.readFileSync(file,'utf8'),data=JSON.parse(raw);
assert.equal(data.kind,'shock-fitted-kernel-study-v101');assert.equal(data.productionQualified,false);
assert.equal(fs.readFileSync(path.join(folder,'ale-data.js'),'utf8'),'window.ALE_STUDY='+raw+';\n');
const context=Object.fromEntries(['scale','clearRect','fillText','fillRect','beginPath','moveTo','lineTo','stroke'].map(k=>[k,(...args)=>assert.ok(args.filter(x=>typeof x==='number').every(Number.isFinite))]));
class E{constructor(id){this.id=id;this.children=[];this.events={};this.attrs={};this.value='';}append(...x){this.children.push(...x);}replaceChildren(...x){this.children=[...x];}addEventListener(k,f){this.events[k]=f;}setAttribute(k,v){this.attrs[k]=String(v);}getBoundingClientRect(){return{left:0,top:0,width:800,height:160};}getContext(){return context;}}
const ids=Object.fromEntries(['caseChoice','title','scope','gasScope','details','frame','time','face','viewer','renderState','role','probe','methods','curve','view','field','csv'].map(k=>[k,new E(k)]));ids.viewer.contentWindow={};
const buttons=Array.from({length:5},()=>new E('button')),listeners={};let blob=null,downloads=0;
const document={getElementById:k=>{assert.ok(ids[k],k);return ids[k];},querySelectorAll:()=>buttons,createElement:tag=>{const e=new E(tag);if(tag==='a')e.click=()=>downloads++;return e;}};
const location={origin:'http://127.0.0.1:8781'},window={ALE_STUDY:data,devicePixelRatio:1,addEventListener:(k,f)=>listeners[k]=f};
let geometryChecks=0;
for(let k=0;k<2;k++)for(const view of ['all','wall','outer','rocket']){
 const c=data.cases[k],frame=c.frames[1],map=c.faces.map((_,i)=>i).filter(i=>!['wall','outer'].includes(view)||c.boundaryKind[i]===(view==='wall'?2:3));
 const originalV=view==='rocket'?data.originalGeometry.vertices:frame.position,originalF=view==='rocket'?data.originalGeometry.faces:map.map(i=>c.faces[i]);let payload=null;
 const geoWindow={ALE_STUDY:data,setup:bridge=>bridge.addEventListener('DataChanged',()=>{payload=bridge.Data;})};
 const script=fs.readFileSync(path.join(folder,'geometry/ale-study-geometry.js'),'utf8');
 vm.runInNewContext(script,{window:geoWindow,location:{origin:location.origin,search:'?case='+k+'&frame=1&revision=1&view='+view+'&field=pressurePa&layer=1&wire=0'},URLSearchParams,parent:{postMessage:()=>{}},document:{querySelector:()=>new E('note')}});
 assert.equal(new Set(payload.faces.flat()).size,payload.vertices.length);assert.equal(payload.faces.length,originalF.length);
 payload.faces.forEach((f,i)=>f.forEach((id,j)=>assert.deepEqual(Array.from(payload.vertices[id]),originalV[originalF[i][j]])));
 assert.equal(payload.faceField.values.length,originalF.length);if(view==='rocket'||k===1)assert.ok(payload.faceField.values.every(v=>v===null));geometryChecks++;
}
vm.runInNewContext(fs.readFileSync(path.join(folder,'ale-study.js'),'utf8'),{window,document,location,Blob,URLSearchParams,URL:{createObjectURL:b=>{blob=b;return'blob:test';},revokeObjectURL:()=>{}}});
const fmt=v=>Number.isFinite(v)?Number(v.toPrecision(8)).toString():'未求解';
function row(id,key){const i=ids[id].children.findIndex(x=>x.textContent===key);assert.ok(i>=0,key);return ids[id].children[i+1].textContent;}
function select(id,value,event='change'){ids[id].value=String(value);ids[id].events[event]();}
function query(){return new URL(ids.viewer.attrs.src,location.origin).searchParams;}
function message(k,f,index,extra={}){return{origin:location.origin,source:ids.viewer.contentWindow,data:{kind:'ale-view-event',caseIndex:k,frameIndex:f,payload:{kind:'geometry-face-probe',revision:Number(query().get('revision')),faceIndex:index,barycentric:[1/3,1/3,1/3]}},...extra};}
function parse(line){return [...line.matchAll(/"((?:[^"]|"")*)"(?:,|$)/g)].map(x=>x[1].replace(/""/g,'"'));}
(async()=>{let workflow=0,faces=0,csvRows=0,ignored=0;
for(let k=0;k<2;k++){select('caseChoice',k);for(let step=0;step<5;step++){buttons[step].events.click();assert.equal(query().get('layer'),step===4?'1':'0');assert.equal(query().get('wire'),step===1||(step===4&&k===0)?'1':'0');workflow++;}
 for(const view of(k===0?['all','wall','outer']:['all','wall','outer','rocket'])){select('view',view);const c=data.cases[k],map=view==='rocket'?data.originalGeometry.faces.map((_,i)=>i):c.faces.map((_,i)=>i).filter(i=>!['wall','outer'].includes(view)||c.boundaryKind[i]===(view==='wall'?2:3));
  for(let f=0;f<2;f++){select('frame',f,'input');for(let i=0;i<map.length;i++){listeners.message(message(k,f,i));const p=view==='rocket'?null:c.frames[f].pressurePa[map[i]];assert.equal(row('probe','流体压力 / Pa'),fmt(p));assert.equal(row('probe','热流 / W/m²'),'未求解');assert.equal(ids.face.value,String(i));faces++;}
   for(const key of ['pressurePa','temperatureK','speedMS','wallPressurePa','heatFluxWm2','wallTemperatureK'])select('field',key);
   ids.csv.events.click();const lines=(await blob.text()).replace(/^\uFEFF/,'').trim().split(/\r?\n/),header=parse(lines[0]);assert.equal(lines.length,map.length+1);
   lines.slice(1).forEach((line,i)=>{const values=parse(line);assert.equal(values.length,header.length);const r=Object.fromEntries(header.map((name,j)=>[name,values[j]]));assert.equal(Number(r.view_face_id),i+1);assert.equal(Number(r.saved_time_s),c.frames[f].timeS);for(const key of ['heat_flux_W_m2','solved_wall_temperature_K','edge_temperature_K','recovery_temperature_K','input_wall_temperature_K'])assert.equal(r[key],'');assert.equal(r.production_qualified,'false');assert.equal(r.heat_flux_status,'unknown');const p=view==='rocket'?null:c.frames[f].pressurePa[map[i]];assert.equal(r.pressure_Pa,p===null?'':String(p));csvRows++;});
  }
 }
}
select('view','all');select('caseChoice',0);select('frame',1,'input');const valid=message(0,1,0),bad=[{...valid,origin:'https://invalid.test'},{...valid,source:{}},{...valid,data:{...valid.data,caseIndex:1}},{...valid,data:{...valid.data,frameIndex:0}}];
for(const change of [{revision:0},{faceIndex:-1},{faceIndex:1e9},{barycentric:[1,1,1]},{barycentric:[NaN,0,1]},{kind:'surface-flow-probe'}])bad.push({...valid,data:{...valid.data,payload:{...valid.data.payload,...change}}});
for(const event of bad){const before=JSON.stringify(ids.probe.children);listeners.message(event);assert.equal(JSON.stringify(ids.probe.children),before);ignored++;}
select('view','rocket');assert.equal(query().get('case'),'1');assert.ok(ids.role.textContent.includes('原完整箭体'));select('caseChoice',0);assert.equal(query().get('view'),'all');
const result={status:'passed-actual-ALE-workbench-workflow-all-face-probes-CSV-and-source-isolation',workflowEvents:workflow,probedFaces:faces,csvRows,downloads,rejectedMessages:ignored,geometrySubsetsVerified:geometryChecks,dataSHA256:crypto.createHash('sha256').update(raw).digest('hex'),productionQualified:false,newProductionSolvedAreaM2:0};if(process.argv[3])fs.writeFileSync(process.argv[3],JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify(result));
})().catch(e=>{console.error(e);process.exitCode=1;});
