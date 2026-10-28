/* 四档真实数据、统一色标、同窗事件和分布曲线核对；没有启动浏览器或执行WebGL。 */
'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),crypto=require('node:crypto'),assert=require('node:assert/strict'),cp=require('node:child_process');
const root=path.resolve(__dirname,'..'),folder=path.resolve(process.argv[2]||'outputs/v76_validation/final_workbench');
const sha=p=>crypto.createHash('sha256').update(fs.readFileSync(p)).digest('hex');
const raw=fs.readFileSync(path.join(folder,'sensitivity-data.json'),'utf8'),d=JSON.parse(raw),processes=[];
assert.equal(fs.readFileSync(path.join(folder,'sensitivity-data.js'),'utf8'),'window.ROCKET_SENSITIVITY='+raw+';\n');
assert.equal(d.cases.length,4);assert.equal(d.gridIndependentClaimed,false);assert.equal(d.newProductionSolvedAreaM2,0);
const fields=new Map();
for(const [i,c] of d.cases.entries()){
 const f=path.join(folder,path.dirname(c.resultsURL)),study=path.dirname(f);
 for(const [script,target] of [['testRocketAxisymmetricView.js',f],['testRocketAxisymmetricWorkbench.js',study]]){
  const result=cp.spawnSync(process.execPath,[path.join(__dirname,script),target],{encoding:'utf8',windowsHide:true});
  assert.equal(result.status,0,result.stdout+result.stderr);processes.push({script,processID:result.pid,exitCode:result.status});
 }
 const data=JSON.parse(fs.readFileSync(path.join(f,'results.json'),'utf8'));
 assert.deepEqual(data.comparisonColorRanges,d.comparisonColorRanges);
 assert.equal(data.checkpointValidation.acceptedEndpointsReplayed,3);
 assert.equal(c.metrics.BEAcceptedSteps,3);assert.ok(c.metrics.originalResidual<=2e-6);
 const near=(a,b)=>assert.ok(Math.abs(a-b)<=2e-12*Math.max(1,Math.abs(a),Math.abs(b)),`${a} != ${b}`);
 near(c.metrics.totalSignedHeatW,data.surfaceIntegralAudit.revolvedWallHeatW);
 for(let j=0;j<data.wallBinding.sourceWallIndex.length;j++){
  const index=data.wallBinding.sourceWallIndex[j]-1;
  near(c.metrics.fields.pressurePa[index],data.wallBinding.sourcePressurePa[j]);
  near(c.metrics.fields.signedHeatWm2[index],data.wallBinding.sourceHeatFluxIntoWallWm2[j]);
  near(c.metrics.fields.wallAdjacentGasProxyK[index],data.wallBinding.sourceWallAdjacentGasProxyK[j]);
  near(c.metrics.fields.wallNormalDistanceM[index],data.wallThermalBinding.normalDistanceM[j]);
 }
 assert.equal(data.initialization.seed.initialGuessIsBoundaryCondition,false);
 if(i>0){assert.equal(data.researchIdentity.caseFactory,'createRocketAxisymmetricSensitivityCase');assert.equal(data.initialization.seed.kind,'v76-sensitivity-spatial-root');assert.equal(data.initialization.seed.priorNewtonIterations,c.metrics.NewtonIterations);}
 const meaning=require(path.join(f,'viewer.js'));
 for(const field of data.fields){
  let v=fields.get(field.key)||[];v.push(...field.values.filter((x,k)=>field.status[k]>0&&Number.isFinite(x)));fields.set(field.key,v);
  const common=meaning.comparisonRange(data,field.key);
  if(field.key==='wallSolved'||field.key==='edgeTemperature'||field.key==='recoveryTemperature')assert.equal(common,null);
  if(common)assert.deepEqual(common,d.comparisonColorRanges[field.key]);
 }
 assert.ok(fs.existsSync(path.join(folder,c.geometryURL)));
}
for(const [key,v] of fields){if(!v.length){assert.equal(d.comparisonColorRanges[key],undefined);continue;}
 const [lo,hi]=d.comparisonColorRanges[key];assert.ok(lo<hi&&lo<=Math.min(...v)&&hi>=Math.max(...v));
}
/* 模块用途：最小DOM仅验证真实页面代码的选择事件、同一iframe及曲线几何；不冒称浏览器。 */
class E{constructor(tag){this.tag=tag;this.children=[];this.attrs={};this.events={};this.value='';this.textContent='';this.hidden=false;}append(...e){this.children.push(...e);}replaceChildren(...e){this.children=[...e];}setAttribute(k,v){this.attrs[k]=String(v);}getAttribute(k){return this.attrs[k];}addEventListener(k,f){this.events[k]=f;}}
const ids=Object.fromEntries(['caseChoice','title','description','details','viewer','iterations','newtonRows','beRows','comparisonRows','comparisonCurve','curveField'].map(k=>[k,new E(k)]));
ids.curveField.value='signedHeatWm2';const buttons=Array.from({length:5},()=>new E('button'));const viewer=ids.viewer;
const document={getElementById:id=>{assert.ok(ids[id],id);return ids[id];},querySelectorAll:q=>{assert.equal(q,'[data-step]');return buttons;},createElement:t=>new E(t),createElementNS:(ns,t)=>{assert.equal(ns,'http://www.w3.org/2000/svg');return new E(t);}};
const source=fs.readFileSync(path.join(root,'web/rocket-sensitivity.js'),'utf8');
vm.runInNewContext(source,{window:{ROCKET_SENSITIVITY:d},document});
assert.equal(ids.comparisonRows.children.length,4);
for(let i=0;i<4;i++){
 ids.caseChoice.value=String(i);ids.caseChoice.events.change();
 for(let stage=0;stage<5;stage++){buttons[stage].events.click();assert.equal(ids.viewer,viewer);assert.equal(viewer.attrs.src,stage===4?d.cases[i].resultsURL:d.cases[i].geometryURL);assert.equal(ids.iterations.hidden,stage!==3);}
 for(const key of ['signedHeatWm2','pressurePa','wallAdjacentGasProxyK','yPlusWallDensityProxy','signedWallTractionAlongWallPa']){
  ids.curveField.value=key;ids.curveField.events.change();const lines=ids.comparisonCurve.children.filter(e=>e.tag==='polyline');assert.equal(lines.length,4);
  lines.forEach((line,k)=>{assert.equal(line.attrs.points.split(' ').length,d.cases[k].metrics.fields[key].length);assert.ok(!/NaN|Infinity|undefined/.test(line.attrs.points));});
 }
}
const report={status:'passed-four-real-case-data-common-color-ranges-existing-probes-curves-CSV-and-same-host-event-contract',
 processID:process.pid,subprocesses:processes,dataSHA256:sha(path.join(folder,'sensitivity-data.json')),
 sourceSHA256:Object.fromEntries(['tests/testRocketSensitivityWorkbench.js','tests/testRocketAxisymmetricView.js','tests/testRocketAxisymmetricWorkbench.js','web/rocket-sensitivity.js','web/viewer.js'].map(p=>[p,sha(path.join(root,p))])),
 cases:4,BEAcceptedEndpoints:12,sameIframeEvents:20,curveVariableSelections:20,commonColorRanges:d.comparisonColorRanges,
 browserVisualOrWebGLRuntimeVerified:false,physicalSolution:false,newProductionSolvedAreaM2:0};
const out=path.join(folder,'node_sensitivity_summary.json');assert.ok(!fs.existsSync(out));fs.writeFileSync(out,JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report));
