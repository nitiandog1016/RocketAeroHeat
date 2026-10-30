/* 模块用途：实际导出上的五步骤、法向坐标、站点身份、消息隔离及逐值CSV验收。 */
'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const folder=path.resolve(process.argv[2]),dataFile=path.join(folder,'shock-layer-data.json'),raw=fs.readFileSync(dataFile,'utf8'),data=JSON.parse(raw);
assert.equal(data.kind,'shock-layer-study-v1');assert.equal(data.productionQualified,false);
assert.equal(fs.readFileSync(path.join(folder,'shock-layer-data.js'),'utf8'),'window.SHOCK_LAYER_STUDY='+raw+';\n');
const context=Object.fromEntries(['setTransform','clearRect','fillText','beginPath','moveTo','lineTo','stroke','arc','fill'].map(name=>[name,(...args)=>assert.ok(args.filter(x=>typeof x==='number').every(Number.isFinite))]));
class E{
 constructor(id){this.id=id;this.children=[];this.attrs={};this.events={};this.value='';this.clientWidth=800;this.clientHeight=140;}
 append(...items){this.children.push(...items);}replaceChildren(...items){this.children=[...items];}addEventListener(name,fn){this.events[name]=fn;}
 setAttribute(name,value){this.attrs[name]=String(value);}getAttribute(name){return this.attrs[name];}getContext(){return context;}
 getBoundingClientRect(){return{left:0,top:0,width:this.clientWidth,height:this.clientHeight};}
}
const ids=Object.fromEntries(['caseChoice','columnChoice','viewChoice','field','viewer','title','scope','details','methodScope','fitRows','renderState','station','stationNote','stateRole','probe','methods','curve','plotNote','saveCSV'].map(name=>[name,new E(name)]));
ids.field.value='pressurePa';ids.viewChoice.value='nose';ids.viewer.contentWindow={};ids.viewer.parentElement={scrollTop:200};
const buttons=Array.from({length:5},()=>new E('button')),listeners={};let blob=null,downloads=0;
const document={getElementById:name=>{assert.ok(ids[name],name);return ids[name];},querySelectorAll:()=>buttons,createElement:tag=>{const e=new E(tag);if(tag==='a')e.click=()=>downloads++;return e;}};
const location={origin:'http://127.0.0.1:8779'},window={SHOCK_LAYER_STUDY:data,devicePixelRatio:1,addEventListener:(name,fn)=>{listeners[name]=fn;}};
vm.runInNewContext(fs.readFileSync(path.join(folder,'shock-layer-study.js'),'utf8'),{window,document,location,Blob,URL:{createObjectURL:b=>{blob=b;return'blob:test';},revokeObjectURL:()=>{}}});
const fmt=x=>Number(x.toPrecision(7)).toString();
function row(node,name){const i=ids[node].children.findIndex(x=>x.textContent===name);assert.ok(i>=0,name);return ids[node].children[i+1].textContent;}
function message(caseIndex,column,station,extra={}){const revision=Number(new URL(ids.viewer.attrs.src,location.origin).searchParams.get('revision'));return{origin:location.origin,source:ids.viewer.contentWindow,data:{kind:'shock-layer-view-event',caseIndex,payload:{kind:'surface-flow-probe',revision,curveIndex:column,stationIndex:station}},...extra};}
function parse(line){return [...line.matchAll(/"((?:[^"]|"")*)"(?:,|$)/g)].map(m=>m[1].replace(/""/g,'"'));}
(async()=>{
 let workflow=0,stations=0,csvRows=0,ignored=0;
 for(let k=0;k<data.cases.length;k++){
  ids.caseChoice.value=String(k);ids.caseChoice.events.change();
  for(let s=0;s<5;s++){ids.viewer.parentElement.scrollTop=200;buttons[s].events.click();const q=new URL(ids.viewer.attrs.src,location.origin).searchParams;assert.equal(q.get('case'),String(k));assert.equal(q.get('layer'),s>=3?'1':'0');assert.equal(q.get('wire'),s===1?'1':'0');assert.equal(ids.viewer.parentElement.scrollTop,0);assert.equal(ids.viewChoice.value,q.get('view'));assert.equal(ids.viewChoice.disabled,s===0);workflow++;}
  const c=data.cases[k];
  for(let j=0;j<c.columns.length;j++){
   const col=c.columns[j];
   for(let i=0;i<col.eta.length;i++){
    listeners.message(message(k,j,i));assert.equal(ids.station.value,String(i));assert.equal(ids.columnChoice.value,String(j));
    assert.equal(row('probe','层内压力 / kPa'),fmt(col.pressurePa[i]*.001));assert.equal(row('probe','层内温度 / K'),fmt(col.temperatureK[i]));
    assert.equal(row('methods','温度方法'),col.methods.temperature);assert.equal(row('probe','热流 / W/m²'),'未求解');
    assert.ok(ids.stationNote.textContent.includes('η = '+fmt(col.eta[i])));if(col.eta[i]===0)assert.ok(ids.stateRole.textContent.includes('不代表真实壁面场'));stations++;
   }
   for(const field of ['pressurePa','temperatureK','speedMS','entropyJPerKgK','jacobianRatio']){ids.field.value=field;ids.field.events.change();}
   ids.saveCSV.events.click();const lines=(await blob.text()).replace(/^\uFEFF/,'').trim().split(/\r?\n/),header=parse(lines[0]);assert.equal(lines.length,col.eta.length+1);
   for(let i=0;i<col.eta.length;i++){
    const cells=parse(lines[i+1]);assert.equal(cells.length,header.length);const r=Object.fromEntries(header.map((name,n)=>[name,cells[n]]));
    assert.equal(Number(r.eta),col.eta[i]);assert.equal(Number(r.normal_distance_m),col.normalDistanceM[i]);assert.equal(Number(r.layer_pressure_Pa),col.pressurePa[i]);assert.equal(Number(r.layer_temperature_K),col.temperatureK[i]);
    for(const field of ['heat_flux_W_m2','wall_temperature_input_K','wall_temperature_solved_K','edge_temperature_K','recovery_temperature_K'])assert.equal(r[field],'');
    assert.equal(r.heat_flux_status,'unknown');assert.equal(r.wall_temperature_status,'unknown');assert.equal(r.pressure_method,col.methods.pressure);assert.equal(r.state_role,col.stateRole);csvRows++;
   }
  }
  const before=ids.station.value;
  for(const m of [message(k,0,0,{origin:'https://unrelated.invalid'}),message(k,0,0,{source:{}}),message(-1,0,0),message(k,999,0),message(k,0,999)]){listeners.message(m);assert.equal(ids.station.value,before);ignored++;}
  const stale=message(k,0,0);stale.data.payload.revision--;listeners.message(stale);assert.equal(ids.station.value,before);ignored++;
  const col=c.columns[Number(ids.columnChoice.value)],i=Number(ids.station.value),old=col.pressurePa[i];col.pressurePa[i]=null;ids.station.events.input();assert.equal(row('probe','层内压力 / kPa'),'未求解');col.pressurePa[i]=old;
 }
 const expectedStations=data.cases.reduce((sum,c)=>sum+c.columns.reduce((n,col)=>n+col.eta.length,0),0),expectedDownloads=data.cases.reduce((sum,c)=>sum+c.columns.length,0);
 assert.equal(stations,expectedStations);assert.equal(csvRows,expectedStations);assert.equal(downloads,expectedDownloads);
 const report={status:'passed-shock-layer-five-step-point-probe-message-isolation-methods-and-exact-CSV',workflowEvents:workflow,stationsChecked:stations,CSVrows:csvRows,CSVdownloads:downloads,ignoredUnrelatedOrStaleMessages:ignored,unknownHeatAndTemperaturesRemainUnknown:true,browserVisualOrWebGLVerified:false,newProductionSolvedAreaM2:0,dataSHA256:crypto.createHash('sha256').update(raw).digest('hex')};
 const target=path.resolve(process.argv[3]);assert.equal(fs.existsSync(target),false);fs.writeFileSync(target,JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report));
})().catch(e=>{console.error(e);process.exitCode=1;});
