/* 真实20度同态功策略与迭代端点：原面、独立变量身份、CSV/截线和同窗事件。 */
'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const folder=path.resolve(process.argv[2]);const raw=fs.readFileSync(path.join(folder,'krylov-study-data.json'),'utf8'),d=JSON.parse(raw);
const meaning=require('../web/viewer.js'),exp=require('../web/export.js'),section=require('../web/section.js');
assert.equal(d.cases.length,5);assert.equal(d.sourceFaces,19752);assert.equal(d.cellCount,217272);assert.equal(d.productionQualified,false);assert.equal(d.newProductionSolvedAreaM2,0);
assert.equal(fs.readFileSync(path.join(folder,'krylov-study-data.js'),'utf8'),'window.ROCKET_KRYLOV_STUDY='+raw+';\n');
const snapshots=[];let CSVrows=0,sectionSegments=0;
for(let j=0;j<d.cases.length;j++){
 const c=d.cases[j],data=JSON.parse(fs.readFileSync(path.join(folder,c.resultsJSONURL),'utf8'));snapshots.push(data);
 if(j>0){assert.deepEqual(data.vertices,snapshots[0].vertices);assert.deepEqual(data.faces,snapshots[0].faces);assert.deepEqual(data.faceArea,snapshots[0].faceArea);}
 assert.equal(data.faces.length,19752);assert.equal(data.fields.length,9);assert.equal(data.researchIdentity.productionQualified,false);assert.equal(data.resultCondition.physicalTimeAdvanced,0);assert.equal(data.resultCondition.BETimeStepsClaimed,0);
 assert.deepEqual(data.comparisonColorRanges,d.commonRanges);assert.equal(meaning.initialField(data),'heatFlux');assert.equal(data.researchIdentity.arithmeticWork,c.arithmeticWork);
 const fields=Object.fromEntries(data.fields.map(f=>[f.key,f]));
 for(const key of ['edgeTemperature','recoveryTemperature','wallSolved']){assert.ok(fields[key].values.every(x=>x===null));assert.ok(fields[key].status.every(x=>x===0));assert.ok(fields[key].missingReason.every(x=>x.length>5));}
 assert.ok(fields.wallInput.status.every(x=>x===2));assert.ok(fields.wallInput.values.every(x=>Math.abs(x-300)<=4*Number.EPSILON*300));
 for(const key of ['pressure','heatFlux','wallAdjacentTemperature','temperatureRate']){
  assert.ok(fields[key].values.every(x=>Number.isFinite(x)));assert.ok(fields[key].status.every(x=>x===1));assert.equal(meaning.stateLabel(fields[key],1),'未收敛研究数值');
 }
 assert.equal(meaning.stateLabel(fields.freestreamNormalProjection,1),'几何诊断');
 let power=0;data.faces.forEach((f,i)=>{assert.ok(f.every(v=>Number.isInteger(v)&&v>=0&&v<data.vertices.length));power+=data.faceArea[i]*fields.heatFlux.values[i]*1000;});
 assert.ok(Math.abs(power-c.signedHeatPowerW)/Math.max(1,Math.abs(power))<3e-13);
 for(const field of data.fields){const csv=exp.csv(data,field,{caseNote:data.caseNote});assert.equal(csv.trim().split(/\r?\n/).length,data.faces.length+1);CSVrows+=data.faces.length;if(field.status[0]===1){assert.ok(!csv.includes('"工程计算"'));}}
 for(const x of [-2.5,-.5,-.2,0]){
  const geometry=section.geometry(data.vertices,data.faces,x);assert.ok(geometry.segments.length>0);
  for(const key of ['pressure','heatFlux','wallAdjacentTemperature','temperatureRate','edgeTemperature']){
   const samples=section.sample(geometry,fields[key]);samples.forEach(s=>{assert.equal(s.value,fields[key].values[s.face]);assert.equal(s.status,fields[key].status[s.face]);});sectionSegments+=samples.length;
   if(key==='heatFlux')assert.ok(section.csv(geometry,fields[key],'未收敛研究',data.caseNote).includes('未收敛研究数值'));
  }
 }
}
const fields=snapshots.map(s=>Object.fromEntries(s.fields.map(f=>[f.key,f])));
for(const key of ['pressure','heatFlux','wallAdjacentTemperature','wallInput'])assert.deepEqual(fields[0][key].values,fields[1][key].values);
assert.ok(fields[0].temperatureRate.values[16339]<0);assert.ok(fields[1].temperatureRate.values[16339]>0);
assert.ok(d.workReport.maximumOldWorkSpeed>3e6);assert.ok(d.workReport.maximumNewWorkSpeed<1600);assert.ok(d.workReport.negativePairCount>0);
class E{constructor(tag){this.tag=tag;this.children=[];this.attrs={};this.events={};this.value='';this.textContent='';}append(...e){this.children.push(...e);}replaceChildren(...e){this.children=[...e];}setAttribute(k,v){this.attrs[k]=String(v);}getAttribute(k){return this.attrs[k];}addEventListener(k,f){this.events[k]=f;}}
const ids=Object.fromEntries(['caseChoice','title','description','details','viewer','historyRows','workNote','stateNote'].map(x=>[x,new E(x)])),buttons=Array.from({length:5},()=>new E('button')),viewer=ids.viewer;
const document={getElementById:x=>{assert.ok(ids[x],x);return ids[x];},querySelectorAll:q=>{assert.equal(q,'[data-step]');return buttons;},createElement:t=>new E(t)};
vm.runInNewContext(fs.readFileSync(path.join(folder,'rocket-krylov.js'),'utf8'),{window:{ROCKET_KRYLOV_STUDY:d},document});let events=0;
for(let j=0;j<5;j++){
 ids.caseChoice.value=String(j);ids.caseChoice.events.change();const h=d.cases[j].history;assert.equal(ids.historyRows.children.length,h&&Array.isArray(h.iteration)?h.iteration.length:0);
 for(let step=0;step<5;step++){buttons[step].events.click();events++;assert.equal(ids.viewer,viewer);assert.equal(viewer.attrs.src,step===4?d.cases[j].resultsURL:step===1?d.domainURL:d.geometryURL);assert.equal(buttons[step].attrs['aria-selected'],'true');}
 assert.ok(ids.workNote.textContent.includes('不能宣称一般正性'));assert.ok(ids.stateNote.textContent.includes('不是壁温历程'));
}
const report={status:'passed-five-actual-facet-cases-work-identity-all-field-CSV-sections-and-same-host-events',processID:process.pid,
 cases:5,sourceFaces:19752,CSVrows,sectionSegments,sameHostEvents:events,browserVisualOrWebGLRuntimeVerified:false,newProductionSolvedAreaM2:0,
 dataSHA256:crypto.createHash('sha256').update(raw).digest('hex')};
fs.writeFileSync(path.join(folder,'node_krylov_summary.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report));
