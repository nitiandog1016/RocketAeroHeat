/* 模块用途：真实三维通量对照的原面数值、算法身份、CSV、截线和五步骤同宿主事件。 */
'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const folder=path.resolve(process.argv[2]),raw=fs.readFileSync(path.join(folder,'exact-study-data.json'),'utf8'),d=JSON.parse(raw);
const meaning=require('../web/viewer.js'),exp=require('../web/export.js'),section=require('../web/section.js');
assert.equal(d.cases.length,3);assert.equal(d.sourceFaces,19752);assert.equal(d.cellCount,217272);assert.equal(d.productionQualified,false);assert.equal(d.newProductionSolvedAreaM2,0);
assert.equal(d.cases[0].summary.residualEvaluations,0);assert.equal(d.cases[0].summary.sourceHLLCAcceptedIterations,6);
assert.equal(d.cases[0].summary.matrixVectorProducts,undefined);assert.equal(d.cases[0].summary.frozenSourceReport.totalAcceptedNonlinearIterations,6);
assert.equal(fs.readFileSync(path.join(folder,'exact-study-data.js'),'utf8'),'window.ROCKET_EXACT_STUDY='+raw+';\n');
const snapshots=[];let CSVrows=0,sectionSamples=0;
for(const c of d.cases){
 const data=JSON.parse(fs.readFileSync(path.join(folder,c.resultsJSONURL),'utf8'));snapshots.push(data);
 assert.deepEqual(data.vertices,snapshots[0].vertices);assert.deepEqual(data.faces,snapshots[0].faces);assert.deepEqual(data.faceArea,snapshots[0].faceArea);
 assert.equal(data.faces.length,19752);assert.equal(data.fields.length,9);assert.deepEqual(data.comparisonColorRanges,d.commonRanges);
 assert.equal(data.researchIdentity.interiorEulerFlux,c.interiorEulerFlux);assert.equal(data.researchIdentity.HLLCHistoryAppended,false);
 assert.equal(data.researchIdentity.productionQualified,false);assert.equal(data.resultCondition.physicalTimeAdvanced,0);assert.equal(data.resultCondition.BETimeStepsClaimed,0);
 const fields=Object.fromEntries(data.fields.map(f=>[f.key,f]));
 for(const key of ['edgeTemperature','recoveryTemperature','wallSolved']){assert.ok(fields[key].values.every(x=>x===null));assert.ok(fields[key].status.every(x=>x===0));assert.ok(fields[key].missingReason.every(x=>x.length>5));}
 assert.ok(fields.wallInput.status.every(x=>x===2));assert.ok(fields.wallInput.values.every(x=>Math.abs(x-300)<=4*Number.EPSILON*300));
 for(const key of ['pressure','heatFlux','wallAdjacentTemperature','temperatureRate']){
  const f=fields[key];assert.ok(f.values.every(x=>Number.isFinite(x)));assert.ok(f.status.every(x=>x===1));assert.equal(meaning.stateLabel(f,1),'未收敛研究数值');
  assert.ok(f.methodNames[0].includes(c.interiorEulerFlux==='hllc'?'HLLC':'Godunov'));
 }
 assert.equal(meaning.stateLabel(fields.freestreamNormalProjection,1),'几何诊断');
 let power=0;data.faces.forEach((f,i)=>{assert.ok(f.every(v=>Number.isInteger(v)&&v>=0&&v<data.vertices.length));power+=data.faceArea[i]*fields.heatFlux.values[i]*1000;});
 assert.ok(Math.abs(power-c.signedHeatPowerW)/Math.max(1,Math.abs(power))<3e-13);
 for(const f of data.fields){const csv=exp.csv(data,f,{caseNote:data.caseNote});assert.equal(csv.trim().split(/\r?\n/).length,data.faces.length+1);CSVrows+=data.faces.length;}
 for(const x of [-2.5,-.5,-.2,0]){
  const geometry=section.geometry(data.vertices,data.faces,x);assert.ok(geometry.segments.length>0);
  for(const key of ['pressure','heatFlux','wallAdjacentTemperature','temperatureRate','edgeTemperature']){
   const samples=section.sample(geometry,fields[key]);samples.forEach(s=>{assert.equal(s.value,fields[key].values[s.face]);assert.equal(s.status,fields[key].status[s.face]);});sectionSamples+=samples.length;
  }
 }
}
const fields=snapshots.map(s=>Object.fromEntries(s.fields.map(f=>[f.key,f])));
for(const key of ['pressure','heatFlux','wallAdjacentTemperature','wallInput'])assert.deepEqual(fields[0][key].values,fields[1][key].values);
assert.notDeepEqual(fields[0].temperatureRate.values,fields[1].temperatureRate.values);
assert.notDeepEqual(fields[1].heatFlux.values,fields[2].heatFlux.values);
assert.equal(d.cases[2].history.iteration[0],0);assert.equal(d.cases[2].summary.lineage.sourceHistoryAppended,false);
assert.ok(d.integrationReport.actualFacePressureStarPa>100);
class E{constructor(tag){this.tag=tag;this.children=[];this.attrs={};this.events={};this.value='';this.textContent='';}append(...e){this.children.push(...e);}replaceChildren(...e){this.children=[...e];}setAttribute(k,v){this.attrs[k]=String(v);}getAttribute(k){return this.attrs[k];}addEventListener(k,f){this.events[k]=f;}}
const ids=Object.fromEntries(['caseChoice','title','description','details','viewer','historyRows','fluxNote','stateNote'].map(x=>[x,new E(x)])),buttons=Array.from({length:5},()=>new E('button')),viewer=ids.viewer;
const document={getElementById:x=>{assert.ok(ids[x],x);return ids[x];},querySelectorAll:q=>{assert.equal(q,'[data-step]');return buttons;},createElement:t=>new E(t)};
vm.runInNewContext(fs.readFileSync(path.join(folder,'rocket-exact.js'),'utf8'),{window:{ROCKET_EXACT_STUDY:d},document});let events=0;
for(let j=0;j<3;j++){
 ids.caseChoice.value=String(j);ids.caseChoice.events.change();const h=d.cases[j].history;
 assert.equal(ids.historyRows.children.length,h&&Array.isArray(h.iteration)?h.iteration.length:0);
 for(let step=0;step<5;step++){buttons[step].events.click();events++;assert.equal(ids.viewer,viewer);assert.equal(viewer.attrs.src,step===4?d.cases[j].resultsURL:step===1?d.domainURL:d.geometryURL);assert.equal(buttons[step].attrs['aria-selected'],'true');}
 assert.ok(ids.fluxNote.textContent.includes('不是壁压'));assert.ok(ids.stateNote.textContent.includes('不计入正式合格面积'));
}
const report={status:'passed-three-real-exact-flux-cases-all-field-CSV-sections-and-same-host-events',processID:process.pid,
 cases:3,sourceFaces:19752,CSVrows,sectionSamples,sameHostEvents:events,browserVisualOrWebGLRuntimeVerified:false,newProductionSolvedAreaM2:0,
 dataSHA256:crypto.createHash('sha256').update(raw).digest('hex')};
fs.writeFileSync(path.join(folder,'node_exact_summary.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report));
