/* 真实三维研究原面、逐变量身份、CSV/截线及同窗事件；不替代浏览器视觉验收。 */
'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const folder=path.resolve(process.argv[2]),root=path.resolve(__dirname,'..');
const raw=fs.readFileSync(path.join(folder,'shell-study-data.json'),'utf8'),d=JSON.parse(raw);
const meaning=require('../web/viewer.js'),exp=require('../web/export.js'),section=require('../web/section.js');
assert.equal(d.kind,'rocket-real-20deg-three-dimensional-startup-v1');assert.equal(d.cases.length,3);
assert.equal(d.cellCount,217272);assert.equal(d.sourceFaces,19752);assert.equal(d.history.length,4);
assert.equal(d.productionQualified,false);assert.equal(d.newProductionSolvedAreaM2,0);
assert.equal(d.advanceReport.steadyConverged,false);assert.equal(d.advanceReport.physicalSolution,false);
assert.equal(fs.readFileSync(path.join(folder,'shell-study-data.js'),'utf8'),'window.ROCKET_SHELL_STUDY='+raw+';\n');
let vertices,faces,areas,CSVrows=0,sectionSegments=0;const snapshots=[];
for(let index=0;index<3;index++){
 const c=d.cases[index],data=JSON.parse(fs.readFileSync(path.join(folder,c.resultsJSONURL),'utf8'));snapshots.push(data);
 if(index===0){vertices=data.vertices;faces=data.faces;areas=data.faceArea;}
 else{assert.deepEqual(data.vertices,vertices);assert.deepEqual(data.faces,faces);assert.deepEqual(data.faceArea,areas);}
 assert.equal(faces.length,d.sourceFaces);assert.equal(data.researchIdentity.newProductionSolvedAreaM2,0);assert.equal(data.researchIdentity.productionQualified,false);
 assert.equal(data.resultCondition.startupTimeIsFlightTime,false);assert.equal(data.resultCondition.steadyConverged,false);
 assert.equal(data.resultCondition.angleOfAttackDeg,20);assert.equal(meaning.initialField(data),'heatFlux');
 assert.deepEqual(data.comparisonColorRanges,d.commonRanges);assert.ok(data.caseNote.includes('未收敛'));
 const fields=Object.fromEntries(data.fields.map(f=>[f.key,f]));assert.equal(data.fields.length,8);
 for(const key of ['edgeTemperature','recoveryTemperature','wallSolved']){
  assert.ok(fields[key].values.every(x=>x===null));assert.ok(fields[key].status.every(x=>x===0));assert.ok(fields[key].missingReason.every(x=>x.length>5));
 }
 // 原墙片面积积分再除面积产生至多几个ulp；保留原映射值，不为精确300改写结果。
 assert.ok(fields.wallInput.values.every(x=>Math.abs(x-300)<=4*Number.EPSILON*300));assert.ok(fields.wallInput.status.every(x=>x===2));
 assert.ok(fields.pressure.status.every(x=>x===(index<2?2:1)));assert.ok(fields.wallAdjacentTemperature.status.every(x=>x===(index<2?2:1)));
 assert.ok(fields.heatFlux.status.every(x=>x===1));assert.equal(meaning.stateLabel(fields.heatFlux,1),'未收敛研究数值');
 assert.equal(meaning.stateLabel(fields.freestreamNormalProjection,1),'几何诊断');
 let heatPower=0,wrong=0;data.faces.forEach((f,i)=>{
  assert.ok(f.every(v=>Number.isInteger(v)&&v>=0&&v<vertices.length));heatPower+=areas[i]*fields.heatFlux.values[i]*1000;
  if(fields.heatFlux.values[i]>0)wrong++;
  if(index<2)assert.ok(Math.abs(fields.wallAdjacentTemperature.values[i]-226.50908361133046)<1e-8);
 });
 assert.ok(Math.abs(heatPower-c.signedHeatPowerW)/Math.max(1,Math.abs(heatPower))<2e-13);
 if(index===0)assert.equal(wrong,184);if(index===1)assert.equal(wrong,0);
 for(const field of data.fields){
  const csv=exp.csv(data,field,{caseNote:data.caseNote});assert.equal(csv.trim().split(/\r?\n/).length,faces.length+1);CSVrows+=faces.length;
  if(field.key==='heatFlux'){assert.ok(csv.includes('未收敛研究数值'));assert.ok(!csv.includes('"工程计算"'));}
  if(['edgeTemperature','recoveryTemperature','wallSolved'].includes(field.key))assert.ok(csv.includes('未求解'));
 }
 for(const x of [-2.5,-.5,-.2,0]){
  const geometry=section.geometry(vertices,faces,x);assert.ok(geometry.segments.length>0);
  for(const key of ['pressure','heatFlux','wallAdjacentTemperature','wallInput','edgeTemperature']){
   const sampled=section.sample(geometry,fields[key]);assert.equal(sampled.length,geometry.segments.length);
   sampled.forEach(s=>{assert.equal(s.value,fields[key].values[s.face]);assert.equal(s.status,fields[key].status[s.face]);});
   if(key==='heatFlux'){const csv=section.csv(geometry,fields[key],'未收敛数值启动',data.caseNote);assert.ok(csv.includes('未收敛研究数值'));assert.ok(!csv.includes('"工程计算"'));}
   sectionSegments+=sampled.length;
  }
 }
 const plan=exp.figurePlan(data,fields.heatFlux,{caseNote:data.caseNote,sceneWidth:270,sceneHeight:105,range:d.commonRanges.heatFlux,rangeHasData:true},(s,size)=>Array.from(s).length*size);
 const text=plan.texts.map(x=>x.text).join('');assert.ok(text.includes('未收敛研究数值'));assert.ok(text.includes('不是飞行时间'));
}
assert.deepEqual(snapshots[0].fields.find(f=>f.key==='pressure').values,snapshots[1].fields.find(f=>f.key==='pressure').values);
class E{constructor(tag){this.tag=tag;this.children=[];this.attrs={};this.events={};this.value='';this.textContent='';}append(...e){this.children.push(...e);}replaceChildren(...e){this.children=[...e];}setAttribute(k,v){this.attrs[k]=String(v);}getAttribute(k){return this.attrs[k];}addEventListener(k,f){this.events[k]=f;}}
const ids=Object.fromEntries(['caseChoice','title','description','details','viewer','historyRows','costNote','heatNote'].map(x=>[x,new E(x)]));
const buttons=Array.from({length:5},()=>new E('button')),viewer=ids.viewer;
const document={getElementById:x=>{assert.ok(ids[x],x);return ids[x];},querySelectorAll:q=>{assert.equal(q,'[data-step]');return buttons;},createElement:t=>new E(t)};
vm.runInNewContext(fs.readFileSync(path.join(folder,'rocket-shell.js'),'utf8'),{window:{ROCKET_SHELL_STUDY:d},document});
assert.equal(ids.historyRows.children.length,4);let events=0;
for(let i=0;i<3;i++){
 ids.caseChoice.value=String(i);ids.caseChoice.events.change();
 for(let step=0;step<5;step++){
  buttons[step].events.click();events++;assert.equal(ids.viewer,viewer);
  assert.equal(viewer.attrs.src,step===4?d.cases[i].resultsURL:step===1?d.domainURL:d.geometryURL);
  assert.equal(buttons[step].attrs['aria-selected'],'true');
 }
}
assert.equal(meaning.stateLabel({computedStateLabel:'  '},1),'工程计算');assert.equal(meaning.stateLabel({},2),'给定条件');
const report={status:'passed-three-real-original-facet-numerical-identity-CSV-sections-and-same-host-events',processID:process.pid,
 cases:3,sourceFaces:d.sourceFaces,CSVrows,sectionSegments,sameHostEvents:events,browserVisualOrWebGLRuntimeVerified:false,
 newProductionSolvedAreaM2:0,dataSHA256:crypto.createHash('sha256').update(raw).digest('hex'),
 sourceSHA256:Object.fromEntries(['web/viewer.js','web/export.js','web/index.html','web/rocket-shell.js','tests/testRocketShellWorkbench.js'].map(p=>[p,crypto.createHash('sha256').update(fs.readFileSync(path.join(root,p))).digest('hex')]))};
fs.writeFileSync(path.join(folder,'node_shell_summary.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report));
