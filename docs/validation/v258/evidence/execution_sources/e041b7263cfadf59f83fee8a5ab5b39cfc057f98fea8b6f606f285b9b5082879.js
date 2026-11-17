/* 模块职责：只核对新增方向层和资格传播；不重复旧工作台交互验收。 */
'use strict';
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert'),root=path.resolve(__dirname,'..'),out=process.argv[2];
const api=require(path.join(root,'web/direction-consistency-v258.js')),data=JSON.parse(fs.readFileSync(path.join(out,'direction.json'),'utf8'));
const context={window:{}};vm.runInNewContext(fs.readFileSync(path.join(root,'outputs/v252_validation/workbench02/case-data.js'),'utf8'),context);
const model=context.window.ROCKET_CASE;require(path.join(root,'web/kinsol-coupled-v257.js')).applyKinsolCoupled(model,JSON.parse(fs.readFileSync(path.join(root,'outputs/v257_validation/workbench01/coupled.json'),'utf8')));
const geometry=JSON.stringify([model.vertices,model.faces]),old=JSON.stringify(model.fields),n=model.fields.length;api.applyDirectionConsistency(model,data);
assert.equal(JSON.stringify([model.vertices,model.faces]),geometry);assert.equal(JSON.stringify(model.fields.slice(0,n)),old);assert.equal(model.fields.length,n+3);
let count=0;for(const field of data.fields){const f=model.fields.find(v=>v.key===field.key);assert.equal(f.methodNames[0],field.method);assert(f.stateNote.includes('非精确导数'));
  assert(f.computedStateLabel.includes('未通过'));for(let i=0;i<19752;i++){assert.equal(f.values[i],field.values[i]);assert.equal(f.status[i],1);count++;}}
const delta=data.fields.find(v=>v.key==='directionEnergyDifference');assert(delta.values.some(x=>x<0)&&delta.values.some(x=>x>0));
assert.equal(data.identity.fixedDirectionHalfStepConsistent,false);assert.equal(data.identity.entireJacobianAccuracyProved,false);
assert.equal(data.identity.newQualifiedAreaM2,0);assert.equal(data.identity.nonlinearUpdates,0);assert.equal(data.identity.KINSOLExecuted,false);
assert.equal(data.identity.profileReadyForCFD,false);assert.equal(data.identity.terminalIsRestoredInitial,true);
console.log(JSON.stringify({passed:true,newDiagnosticValuesChecked:count,originalGeometryAndExistingFieldsUnchanged:true,
  unreliableDifferenceNotPublishedAsExactDerivative:true,signedDifferencePreserved:true,browserInteractionTested:false}));
