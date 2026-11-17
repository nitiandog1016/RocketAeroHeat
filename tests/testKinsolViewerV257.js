/* 模块职责：逐面核对返回态/差值/算法及未知温度身份，拒绝把失败控制显示为新解。 */
'use strict';
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert'),root=path.resolve(__dirname,'..'),out=process.argv[2];
const api=require(path.join(root,'web/kinsol-coupled-v257.js')),data=JSON.parse(fs.readFileSync(path.join(out,'coupled.json'),'utf8'));
const context={window:{}};vm.runInNewContext(fs.readFileSync(path.join(root,'outputs/v252_validation/workbench02/case-data.js'),'utf8'),context);
const m=context.window.ROCKET_CASE,F=JSON.stringify(m.faces),V=JSON.stringify(m.vertices),oldFields=JSON.stringify(m.fields),n=m.fields.length;
api.applyKinsolCoupled(m,data);assert.equal(JSON.stringify(m.faces),F);assert.equal(JSON.stringify(m.vertices),V);assert.equal(JSON.stringify(m.fields.slice(0,n)),oldFields);
assert.equal(m.fields.length,n+data.fields.length);let checked=0;
for(const input of data.fields){const f=m.fields.find(v=>v.key===input.key);assert.equal(f.methodNames[0],input.method);assert(f.methodNames[0].length>20);
  assert(f.computedStateLabel.includes('未收敛'));for(let i=0;i<19752;i++){assert.equal(f.values[i],input.values[i]);assert.equal(f.status[i],1);checked++;}}
const q=m.fields.find(f=>f.key==='kinsolReturnedQ');assert(q.values.some(x=>x<0));assert(q.values.some(x=>x>0));
assert.equal(data.identity.actualObservedStateChanges,0);assert.equal(data.identity.nativeNonlinearAttempts,1);assert.equal(data.identity.linearGateRefused,true);
for(const key of ['kinsolDeltaQ','kinsolDeltaP'])assert(m.fields.find(f=>f.key===key).values.every(x=>x===0));
assert.equal(data.identity.solvedSolidWallTemperatureK,null);assert.equal(data.identity.profileReadyForCFD,false);assert.equal(data.identity.newQualifiedAreaM2,0);
assert.equal(data.identity.nativeFlag,-13);assert.equal(data.identity.physicalQualification,false);assert(m.caseNote.includes('接受更新0'));
console.log(JSON.stringify({passed:true,checkedOriginalFaceValues:checked,oldGeometryAndFieldsUnchanged:true,zeroDifferenceIsPhysicalComparison:true,
  signedHeatPreserved:true,rejectedLinearStepNotPublished:true,unknownSolidTemperaturePreserved:true,browserInteractionTested:false}));
