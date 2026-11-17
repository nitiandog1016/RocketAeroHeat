/* 模块职责：全原面诊断值与未知身份核对，不冒称浏览器交互或物理解。 */
'use strict';
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert'),root=path.resolve(__dirname,'..'),out=process.argv[2];
const api=require(path.join(root,'web/full-residual-v256.js')),data=JSON.parse(fs.readFileSync(path.join(out,'residual.json'),'utf8'));
const context={window:{}};vm.runInNewContext(fs.readFileSync(path.join(root,'outputs/v252_validation/workbench02/case-data.js'),'utf8'),context);
const model=context.window.ROCKET_CASE,F=JSON.stringify(model.faces),V=JSON.stringify(model.vertices),count=model.fields.length;api.applyFullResidual(model,data);
assert.equal(JSON.stringify(model.faces),F);assert.equal(JSON.stringify(model.vertices),V);assert.equal(model.fields.length,count+5);
let checked=0;for(const field of data.fields){const f=model.fields.find(x=>x.key===field.key);assert.equal(f.methodNames[0],field.method);assert.equal(f.values.length,19752);
 for(let i=0;i<19752;i++){assert.equal(f.values[i],field.values[i]);assert.equal(f.status[i],1);checked++;}assert(f.computedStateLabel.includes('未收敛'));}
const energy=data.fields.find(f=>f.key==='nativePhysicalEnergyResidual').values;assert(energy.some(x=>x<0)&&energy.some(x=>x>0));
assert.equal(data.identity.physicalQualification,false);assert.equal(data.identity.nonlinearIterations,0);assert.equal(data.identity.timeUpdates,0);assert.equal(data.identity.solvedSolidWallTemperatureK,null);
assert.equal(data.identity.newQualifiedAreaM2,0);assert.equal(data.identity.profileReadyForCFD,false);assert(model.caseNote.includes('3904'));
console.log(JSON.stringify({passed:true,checkedDiagnosticValues:checked,originalGeometryUnchanged:true,signedValuesPreserved:true,zeroUpdatesExplicit:true,browserInteractionTested:false}));
