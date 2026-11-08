/* 模块职责：核对真实导出字段的科研身份、面积与未知/给定壁温，不启动浏览器。 */
'use strict';
const fs=require('fs'),path=require('path'),assert=require('assert');
const root=path.resolve(__dirname,'..');
const viewer=require(path.join(root,'web/viewer.js'));
const file=path.join(root,'outputs/v222_validation/workbench01/rocket_phase_research.json');
const m=JSON.parse(fs.readFileSync(file,'utf8'));
assert.strictEqual(m.kind,'aerothermal-model');assert.strictEqual(m.faces.length,19752);
assert(Math.abs(m.faceArea.reduce((s,a)=>s+a,0)-1.9646292565715553)<1e-12);
assert.strictEqual(m.researchIdentity.productionQualified,false);
assert.strictEqual(viewer.initialField(m,null),'heatFlux');
for(const f of m.fields){
  assert.strictEqual(f.values.length,m.faces.length);assert.strictEqual(f.status.length,m.faces.length);
  assert.strictEqual(f.methodIndex.length,m.faces.length);
  if(f.key==='wallSolved'){
    assert(f.values.every(x=>x===null));assert(f.status.every(x=>x===0));
    assert(viewer.missingReason(m,f,0).includes('未'));
  }else if(f.key==='wallInput'){
    assert(f.values.every(x=>x===300));assert(f.status.every(x=>x===2));
    assert.strictEqual(viewer.stateLabel(f,2),'给定条件');
  }else{
    assert(f.values.every(Number.isFinite));assert(f.status.every(x=>x===1));
    assert.strictEqual(viewer.stateLabel(f,1),'未收敛研究数值');
    assert.strictEqual(viewer.coverageLabels(m,f)[0],'研究数值面积');
    assert(viewer.method(f,0).includes('SA-neg'));
  }
}
assert(m.fields.find(f=>f.key==='heatFlux').values.some(x=>x<0));
console.log(JSON.stringify({passed:true,faces:m.faces.length,fields:m.fields.length,physicalQualification:false,browserInteractionTested:false}));
