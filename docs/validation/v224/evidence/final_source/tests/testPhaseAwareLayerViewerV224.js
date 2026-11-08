/* 模块职责：确认主导出路径的物理身份、原面积和未知厚度/壁温，不启动浏览器。 */
'use strict';
const fs=require('fs'),path=require('path'),assert=require('assert');
const root=path.resolve(__dirname,'..'),view=require(path.join(root,'web/viewer.js'));
const m=JSON.parse(fs.readFileSync(path.join(root,'outputs/v224_validation/matlab01/workbench.json'),'utf8'));
assert.strictEqual(m.kind,'aerothermal-model');assert.strictEqual(m.faces.length,81);
assert.strictEqual(m.fields.length,8);assert.strictEqual(m.researchIdentity.productionQualified,false);
assert(m.resultCondition.caption.includes('prescribed-control'));
assert.strictEqual(view.initialField(m),'heatFlux');
assert(Math.abs(m.faceArea.reduce((s,v)=>s+v,0)-m.totalArea)<1e-14);
for(const f of m.fields){
  assert.strictEqual(f.values.length,81);assert.strictEqual(f.status.length,81);
  assert.strictEqual(f.methodIndex.length,81);assert(!f.history);
  if(['wallSolved','energyThickness'].includes(f.key)){
    assert(f.values.every(x=>x===null)&&f.status.every(x=>x===0));
    assert(view.missingReason(m,f,0));
  }else if(f.key==='wallInput'){
    assert(f.values.every(x=>x===300)&&f.status.every(x=>x===2));
  }else{
    assert.strictEqual(f.values[0],null);assert.strictEqual(f.status[0],0);
    assert(f.values.slice(1).every(Number.isFinite));assert(f.status.slice(1).every(x=>x===1));
    assert.strictEqual(view.stateLabel(f,1),'未验证工程候选');
    assert.strictEqual(view.coverageLabels(m,f)[0],'研究数值面积');
  }
}
assert(!m.fields.some(f=>f.key==='pressure')); // 旧压力对象不能与新相态链静默拼接。
const q=m.fields.find(f=>f.key==='heatFlux');
const rows=fs.readFileSync(path.join(root,'outputs/v224_validation/matlab01/uniform.csv'),'utf8').trim().split(/\r?\n/).slice(1);
rows.slice(1).forEach((r,i)=>assert(Math.abs(q.values[i+1]-Number(r.split(',')[13]))<1e-9));
console.log(JSON.stringify({passed:true,fields:8,faces:81,physicalQualification:false,browserInteractionTested:false}));
