/* 原图实验对照、逐场量身份及同窗五步骤事件；不冒充浏览器/WebGL视觉验收。 */
'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const folder=path.resolve(process.argv[2]),root=path.resolve(__dirname,'..');
const raw=fs.readFileSync(path.join(folder,'boattail-study-data.json'),'utf8'),d=JSON.parse(raw);
assert.equal(fs.readFileSync(path.join(folder,'boattail-study-data.js'),'utf8'),'window.NASA_BOATTAIL_STUDY='+raw+';\n');
assert.equal(d.cases.length,3);assert.equal(d.comparison.pointCount,15);assert.equal(d.comparison.experimentalComparisonPerformed,true);
assert.equal(d.comparison.completeNSRANSExperimentalValidation,false);assert.equal(d.newProductionSolvedAreaM2,0);
const source=JSON.parse(fs.readFileSync(path.join(root,'tests/data/nasa_tnd2761/pressure_digitization.json'),'utf8'));
const meaning=require(path.join(root,'web/viewer.js')),exp=require(path.join(root,'web/export.js')),section=require(path.join(root,'web/section.js'));
let intersected=0,CSVrows=0;
for(const c of d.cases){
 const data=JSON.parse(fs.readFileSync(path.join(folder,c.resultsJSONURL),'utf8'));assert.equal(data.faces.length,c.faceCount);
 assert.equal(meaning.initialField(data),'pressureRatio');assert.deepEqual(data.comparisonColorRanges.pressureRatio,d.commonRatioRange);
 assert.equal(data.researchIdentity.productionQualified,false);assert.equal(data.researchIdentity.experimentalSymbolsAssignedToFacets,false);
 const ratio=data.fields.find(f=>f.key==='pressureRatio');let area=0,count=0;
 data.faces.forEach((face,i)=>{
  assert.ok(face.every(v=>Number.isInteger(v)&&v>=0&&v<data.vertices.length));
  if(data.region[i]==='boattail'){assert.equal(ratio.values[i],c.comparison.predictionRatio);assert.equal(ratio.status[i],1);area+=data.faceArea[i];count++;}
  else {assert.equal(ratio.values[i],null);assert.equal(ratio.status[i],0);}
 });
 assert.ok(count>0&&Math.abs(area-c.candidateRatioEvaluatedFacetAreaM2)<1e-13);
 for(const f of data.fields.filter(f=>f.key!=='pressureRatio')){
  assert.ok(f.values.every(x=>x===null)&&f.status.every(x=>x===0),f.key+' must remain unknown');
  assert.ok(f.missingReason.every(x=>typeof x==='string'&&x.length));
 }
 c.comparison.points.forEach(p=>{
  const s=source.points.find(x=>x.id===p.id);assert.ok(s);
  for(const key of ['pressureRatio','pixelXY','sourceFigureXprimeOverD','ratioDigitizationHalfWidth','matchedOrificeID'])assert.deepEqual(p[key],s[key]);
 });
 const surfaceCSV=exp.csv(data,ratio,{caseNote:data.caseNote});assert.equal(surfaceCSV.trim().split(/\r?\n/).length,data.faces.length+1);CSVrows+=data.faces.length;
 assert.ok(surfaceCSV.includes('条件')&&surfaceCSV.includes('未求解'));
 // 使用同一实际三角网格的平面截线，核对候选尾锥曲线与上游未知曲线没有互相填补。
 const tail=section.geometry(data.vertices,data.faces,c.geometry.boattailLengthM/2);
 assert.ok(tail.segments.length>0);assert.ok(section.sample(tail,ratio).every(s=>s.status===1&&s.value===c.comparison.predictionRatio));
 const upstream=section.geometry(data.vertices,data.faces,c.geometry.boattailLengthM+.02);
 assert.ok(upstream.segments.length>0&&section.sample(upstream,ratio).every(s=>s.status===0&&s.value===null));
 intersected+=tail.segments.length;
}
class E{constructor(tag){this.tag=tag;this.children=[];this.attrs={};this.events={};this.value='';this.textContent='';this.checked=false;}append(...e){this.children.push(...e);}replaceChildren(...e){this.children=[...e];}setAttribute(k,v){this.attrs[k]=String(v);}getAttribute(k){return this.attrs[k];}addEventListener(k,f){this.events[k]=f;}}
const ids=Object.fromEntries(['caseChoice','title','description','details','viewer','allCases','pressureCurve','pointDetails','comparisonRows'].map(x=>[x,new E(x)]));ids.allCases.checked=true;
const buttons=Array.from({length:5},()=>new E('button')),viewer=ids.viewer;
const document={getElementById:x=>{assert.ok(ids[x],x);return ids[x];},querySelectorAll:q=>{assert.equal(q,'[data-step]');return buttons;},createElement:t=>new E(t),createElementNS:(ns,t)=>{assert.equal(ns,'http://www.w3.org/2000/svg');return new E(t);}};
vm.runInNewContext(fs.readFileSync(path.join(folder,'nasa-boattail.js'),'utf8'),{window:{NASA_BOATTAIL_STUDY:d},document});
let events=0,clicks=0;
for(let i=0;i<3;i++){
 ids.caseChoice.value=String(i);ids.caseChoice.events.change();
 for(let step=0;step<5;step++){buttons[step].events.click();events++;assert.equal(ids.viewer,viewer);assert.equal(viewer.attrs.src,step===4?d.cases[i].resultsURL:d.cases[i].geometryURL);}
 for(const all of [true,false]){
  ids.allCases.checked=all;ids.allCases.events.change();const circles=ids.pressureCurve.children.filter(x=>x.tag==='circle');assert.equal(circles.length,all?15:5);
  const candidates=ids.pressureCurve.children.filter(x=>x.attrs['data-role']==='candidate');assert.equal(candidates.length,all?3:1);
  circles.forEach(e=>{for(const key of ['cx','cy'])assert.ok(Number.isFinite(Number(e.attrs[key])));e.events.click();assert.ok(ids.pointDetails.textContent.includes('绝对压力/Tw未知'));clicks++;});
 }
}
const report={status:'passed-three-geometry-pressure-ratio-source-identity-CSV-section-and-same-host-events',processID:process.pid,cases:3,
 experimentalSymbols:15,sameHostEvents:events,observationClicks:clicks,surfaceCSVrows:CSVrows,actualSectionSegments:intersected,
 explicitDefaultField:'pressureRatio',browserVisualOrWebGLRuntimeVerified:false,newProductionSolvedAreaM2:0,
 dataSHA256:crypto.createHash('sha256').update(raw).digest('hex'),sourceSHA256:Object.fromEntries(['web/viewer.js','web/nasa-boattail.js','tests/testNasaBoattailWorkbench.js'].map(p=>[p,crypto.createHash('sha256').update(fs.readFileSync(path.join(root,p))).digest('hex')]))};
fs.writeFileSync(path.join(folder,'node_boattail_summary.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report));
