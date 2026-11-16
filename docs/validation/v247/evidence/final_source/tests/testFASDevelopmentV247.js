/* 模块职责：真实三阶段原面数据、固定色标、算法来源和未知固体温度的前端合同。 */
'use strict';
const fs=require('fs'),path=require('path'),assert=require('assert');
const root=path.resolve(__dirname,'..'),viewer=require(path.join(root,'web/viewer.js'));
const base=path.join(root,'outputs/v247_validation'),m=JSON.parse(fs.readFileSync(path.join(base,'workbench01/rocket_fas_control.json'),'utf8'));
assert.strictEqual(m.faces.length,19752);assert.strictEqual(m.researchIdentity.physicalQualification,false);
assert.strictEqual(m.researchIdentity.initialization,'new-compatible-uniform');const stages=m.researchIdentity.successfulUpdates;assert.strictEqual(viewer.initialField(m),'heatFlux_'+stages[stages.length-1]);
assert(Math.abs(m.faceArea.reduce((a,b)=>a+b,0)-1.9646292565715553)<1e-12);
for(const step of stages){
  const data=fs.readFileSync(path.join(base,'audit03','original_faces_'+String(step).padStart(6,'0')+'.csv'),'utf8').trim().split(/\r?\n/).slice(1).map(s=>s.split(',').map(Number));
  for(const [key,col] of [['pressure',2],['heatFlux',3],['wallFluid',4],['wallDeviation',null]]){
    const f=m.fields.find(x=>x.key===key+'_'+step);assert(!f.history);assert(f.status.every(x=>x===1));
    assert.deepStrictEqual(viewer.comparisonRange(m,f.key),viewer.comparisonRange(m,key+'_0'));
    assert.strictEqual(viewer.stateLabel(f,1),'未收敛研究数值');assert(viewer.stateNote(f).includes('不是物理时间'));
    assert(viewer.method(f,0).includes('FAS')&&viewer.method(f,0).includes('SA-neg'));
    data.forEach((r,i)=>assert.strictEqual(f.values[i],col===null?r[4]-r[5]:r[col]));
  }
}
const given=m.fields.find(f=>f.key==='wallInput'),unknown=m.fields.find(f=>f.key==='wallSolved');
assert(given.values.every(x=>x===300)&&given.status.every(x=>x===2));assert(unknown.values.every(x=>x===null)&&unknown.status.every(x=>x===0));
assert(viewer.missingReason(m,unknown,0).includes('未'));
console.log(JSON.stringify({passed:true,faces:19752,fields:m.fields.length,values:19752*4*stages.length,browserInteractionTested:false}));
