/* 模块职责：检查四阶段共同色标、方法/迭代身份与全部原面数据；不冒称浏览器验收。 */
'use strict';
const fs=require('fs'),path=require('path'),assert=require('assert');
const root=path.resolve(__dirname,'..'),viewer=require(path.join(root,'web/viewer.js'));
const m=JSON.parse(fs.readFileSync(path.join(root,'outputs/v223_validation/workbench01/rocket_phase_comparison.json'),'utf8'));
assert.strictEqual(m.faces.length,19752);assert.strictEqual(m.fields.length,14);
assert.strictEqual(m.researchIdentity.productionQualified,false);
assert.strictEqual(viewer.initialField(m),'wallDeviation_200');
assert(Math.abs(m.faceArea.reduce((a,b)=>a+b,0)-1.9646292565715553)<1e-12);
for(const kind of ['pressure','heatFlux','wallDeviation']){
  const group=m.fields.filter(f=>f.comparisonVariable===kind);
  assert.deepStrictEqual(group.map(f=>f.successfulUpdates),[51,101,151,200]);
  const range=viewer.comparisonRange(m,group[0].key);
  for(const f of group){
    assert(!f.history);assert.deepStrictEqual(viewer.comparisonRange(m,f.key),range);
    assert.strictEqual(f.values.length,19752);assert(f.values.every(x=>Number.isFinite(x)&&x>=range[0]&&x<=range[1]));
    assert(f.status.every(x=>x===1));assert(f.label.includes(String(f.successfulUpdates)));
    assert.strictEqual(viewer.stateLabel(f,1),'未收敛研究数值');assert(viewer.method(f,0).includes('SA-neg'));
    assert(viewer.stateNote(f).includes('不是物理时间'));
    const file=path.join(root,'outputs/v222_validation/audit01','original_faces_'+String(f.nativeIteration).padStart(6,'0')+'.csv');
    const rows=fs.readFileSync(file,'utf8').trim().split(/\r?\n/).slice(1);
    rows.forEach((r,i)=>{const a=r.split(',').map(Number);assert.strictEqual(f.values[i],kind==='pressure'?a[2]:kind==='heatFlux'?a[3]:a[4]-a[5]);});
  }
  if(kind==='wallDeviation')assert.strictEqual(range[0],-range[1]);
}
const given=m.fields.find(f=>f.key==='wallInput'),unknown=m.fields.find(f=>f.key==='wallSolved');
assert(given.values.every(x=>x===300)&&given.status.every(x=>x===2));
assert(unknown.values.every(x=>x===null)&&unknown.status.every(x=>x===0));
assert(viewer.missingReason(m,unknown,0).includes('未'));
console.log(JSON.stringify({passed:true,faces:19752,comparisonFields:12,physicalQualification:false,browserInteractionTested:false}));
