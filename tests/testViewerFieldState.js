/* 逐变量科学身份合同：旧JSON兼容、未知上游/本面区别、壁温输入与未执行求解、冻结工况与科研导出。 */
'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),meaning=require('../web/viewer.js'),exp=require('../web/export.js');
const f={key:'heatFlux',label:'对流热流密度',values:[0,null,null],status:[1,0,0],methodIndex:[0,1,1],methodNames:['计算得到精确零值','未计算'],missingReason:['','本面元外缘热力状态未知/域外占位','上游外缘热力状态未知；即使本面有效也不放行']};
const model={missingReason:['已计算','热流总原因','热流总原因'],resultCondition:{caption:'冻结快照：高度30 km；p∞1197 Pa；T∞226.5 K；统一壁温边界输入300 K（不是壁温解）'}};
assert.equal(meaning.state(f,0),1);assert.equal(meaning.missingReason(model,f,0),'');
assert.match(meaning.missingReason({missingReason:['已计算']},{...f,values:[NaN],status:[1],missingReason:undefined},0),/没有有限数值/);
assert.match(meaning.missingReason({missingReason:['给定条件']},{...f,values:[null],status:[2],missingReason:undefined},0),/按未求解显示/);
assert.match(meaning.missingReason(model,f,1),/本面元/);assert.match(meaning.missingReason(model,f,2),/上游/);
assert.equal(meaning.coverageTitle(f),'对流热流密度 · 面积覆盖');assert.equal(meaning.showHeatGaps(f),true);
assert.deepEqual(meaning.coverageLabels(model,f),['计算面积','给定面积','未解面积']);
const researchModel={...model,researchIdentity:{productionQualified:false}};
assert.deepEqual(meaning.coverageLabels(researchModel,f),['研究数值面积','给定面积','未解面积']);
assert.equal(meaning.coverageLabels(researchModel,{...f,computedStateLabel:'几何诊断'})[0],'诊断数值面积');
assert.match(meaning.coverageNote(researchModel,f),/不代表正式合格求解面积/);
assert.equal(meaning.coverageNote(model,f),meaning.stateNote(f));
assert.equal(meaning.showHeatGaps({...f,history:[]}),true,'MATLAB空history数组不是时间历程');
const edge={...f,key:'edgeTemperature',label:'边界层外缘温度',values:[400,null,100],status:[1,0,1],methodNames:['NASA7静温','NASA7尝试'],missingReason:['','外缘域外占位，不是已接受静温','']};
assert.equal(meaning.showHeatGaps(edge),false);assert.equal(meaning.missingReason(model,edge,2),'');
assert.equal(meaning.method(edge,1),'未接受的尝试模型：NASA7尝试');assert.equal(meaning.method(edge,2),'NASA7尝试');
const input={...f,key:'wallInput',label:'热流计算所用壁温',values:[300,null,300],status:[2,0,2],missingReason:undefined};
assert.equal(meaning.state(input,0),2);assert.match(meaning.stateNote(input),/不是本次求解壁温/);assert.match(meaning.missingReason(model,input,1),/输入壁温未知/);
assert.equal(meaning.showHeatGaps(input),false);
const solved={...f,key:'wallSolved',label:'求解壁温',values:[null,null,null],status:[0,0,0],availability:'notProvidedInSnapshot',missingReason:undefined};
assert.match(meaning.stateNote(solved),/未执行/);assert.equal(meaning.showHeatGaps(solved),false);assert.match(meaning.missingReason(model,solved,1),/独立壁温解/);
assert.equal(meaning.detailTitle(solved,1),'未执行壁温求解的说明');assert.equal(meaning.detailTitle(input,1),'输入缺失与使用范围');assert.equal(meaning.detailTitle(edge,2),'算法与适用范围');
assert.equal(meaning.caseCaption(model,f),model.resultCondition.caption);
const legacy={...edge,missingReason:undefined};assert.match(meaning.missingReason(model,legacy,1),/本字段具体未解原因/);assert.ok(!meaning.missingReason(model,legacy,1).includes('热流总原因'));
assert.match(meaning.caseCaption({caseNote:'旧M5工况'},legacy),/缺少完整冻结工况元数据/);
const timed={...f,history:{timeNote:['历史站：高度20 km','历史站：高度40 km'],faceNote:['失效历程','失效历程','失效历程']}};
assert.equal(meaning.showHeatGaps(timed),false);assert.match(meaning.caseCaption(model,timed,1),/高度40 km/);assert.ok(!meaning.caseCaption(model,timed,1).includes('1197'));
assert.equal(meaning.missingReason(model,timed,1),'失效历程');
// 小画布PNG和CSV均通过已有工况列继承冻结信息，不修改原场量或状态。
const mesh={...model,vertices:[[0,0,0],[1,0,0],[0,1,0]],faces:[[0,1,2]],faceArea:[.5],region:['侧壁'],coordinateNote:'SI',caseNote:model.resultCondition.caption};
const one={...input,unit:'K',values:[300],status:[2],methodIndex:[0],methodNames:['壁温输入'],scope:'不是壁温解'};
const csv=exp.csv(mesh,one,{caseNote:meaning.caseCaption(mesh,one)});assert.match(csv,/高度30 km/);assert.match(csv,/p∞1197 Pa/);assert.match(csv,/不是壁温解/);
const plan=exp.figurePlan(mesh,one,{caseNote:meaning.caseCaption(mesh,one),sceneWidth:270,sceneHeight:105,range:[299,301],rangeHasData:true,colors:[[0,0,0],[255,255,255]]},(s,size)=>Array.from(s).length*size);
const text=plan.texts.map(x=>x.text).join('');assert.match(text,/高度30 km/);assert.match(text,/统一壁温边界输入300 K/);
// 可选真实新JSON和旧JSON对照：每场量数值/状态逐值不变，新增说明不得篡改求解覆盖。
if(process.argv[2]){
  const actual=JSON.parse(fs.readFileSync(process.argv[2],'utf8'));const fields=Object.fromEntries(actual.fields.map(f=>[f.key,f]));
  assert.ok(actual.resultCondition&&actual.resultCondition.altitudeM===30000);assert.equal(actual.resultCondition.Mach,5);
  assert.equal(actual.resultCondition.wallInputType,'uniform');assert.equal(actual.resultCondition.wallUniformValueK,300);
  assert.match(actual.caseNote,/p∞/);assert.match(actual.caseNote,/T∞/);assert.match(actual.caseNote,/高度 30 km/);
  assert.equal(fields.wallInput.prescribedAreaFraction,1);assert.equal(fields.wallSolved.missingAreaFraction,1);assert.equal(meaning.showHeatGaps(fields.wallSolved),false);
  assert.equal(meaning.showHeatGaps(fields.heatFlux),true);
  assert.match(meaning.stateNote(fields.wallSolved),/未执行/);
  assert.match(meaning.missingReason(actual,fields.heatFlux,14260),/本面元外缘/);assert.match(meaning.missingReason(actual,fields.heatFlux,21),/上游外缘/);
  assert.equal(meaning.state(fields.edgeTemperature,21),1);assert.equal(meaning.state(fields.heatFlux,21),0);
  assert.match(meaning.method(fields.edgeTemperature,14260),/未接受/);
  if(process.argv[3]){
    const before=JSON.parse(fs.readFileSync(process.argv[3],'utf8'));assert.deepEqual(actual.vertices,before.vertices);assert.deepEqual(actual.faces,before.faces);assert.deepEqual(actual.faceArea,before.faceArea);
    for(const field of actual.fields){const old=before.fields.find(f=>f.key===field.key);assert.ok(old);assert.deepEqual(field.values,old.values,field.key+' values');assert.deepEqual(field.status,old.status,field.key+' status');assert.equal(field.computedAreaFraction,old.computedAreaFraction);}
  }
  console.log('Actual saved snapshot metadata verified: '+actual.faces.length+' faces; unchanged field values/status; direct/upstream unknown distinguished; snapshot conditions recorded. No physical recomputation.');
}
// 专用研究默认先显示实际可用变量；已选变量继续保留，旧快照和非法偏好保持原行为。
const preference={fields:[{key:'pressureRatio'},{key:'heatFlux'},{key:'wallInput'}],defaultField:'pressureRatio'};
assert.equal(meaning.initialField(preference),'pressureRatio');
assert.equal(meaning.initialField(preference,'wallInput'),'wallInput');
assert.equal(meaning.initialField({...preference,defaultField:'invalid'}),'heatFlux');
assert.equal(meaning.initialField({fields:[{key:'pressure'}]}),'pressure');
assert.equal(meaning.initialField({fields:preference.fields}),'heatFlux');
console.log('Viewer field-state contract passed: per-variable coverage, zero/missing/input identity, rejected attempted model, legacy JSON fallback, history separation, explicit default field and frozen PNG/CSV conditions.');
