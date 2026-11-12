/* 模块职责：核对三阶段原面科研字段、共同色标、来源与未知热量，不冒充浏览器验收。 */
'use strict';
const fs=require('fs'),path=require('path'),assert=require('assert');
const root=path.resolve(__dirname,'..'),viewer=require(path.join(root,'web/viewer.js'));
const base=path.join(root,'outputs/v232_validation'),m=JSON.parse(fs.readFileSync(path.join(base,'workbench01/rocket_euler_research.json'),'utf8'));
assert.strictEqual(m.fields.length,13);assert.strictEqual(m.faces.length,19752);assert.strictEqual(viewer.initialField(m),'pressure_800');
assert.strictEqual(m.researchIdentity.physicalQualification,false);assert.strictEqual(m.researchIdentity.normalCompletion,true);
assert.deepStrictEqual(m.researchIdentity.selectedUpdates,[50,400,800]);assert.strictEqual(m.researchIdentity.lastPreprocessingRejectedPoints,8);
const columns={pressure:2,fluidTemperature:3,totalEnthalpy:4,entropy:5,velocityX:6,velocityY:7,velocityZ:8},frames={};
for(const s of [50,400,800])frames[s]=fs.readFileSync(path.join(base,'audit01','original_faces_'+String(s).padStart(6,'0')+'.csv'),'utf8').trim().split(/\r?\n/).slice(1).map(x=>x.split(',').map(Number));
for(const f of m.fields.slice(0,11)){
  assert(f.status.every(s=>s===1)&&f.values.length===19752);assert(!f.history);assert(viewer.method(f,0).includes('Euler'));
  f.values.forEach((v,i)=>assert.strictEqual(v,frames[f.successfulUpdates][i][columns[f.comparisonVariable]]));
}
for(const k of ['pressure','fluidTemperature'])for(const s of [400,800])assert.deepStrictEqual(viewer.comparisonRange(m,k+'_50'),viewer.comparisonRange(m,k+'_'+s));
for(const f of m.fields.slice(11)){assert(f.values.every(v=>v===null)&&f.status.every(s=>s===0));assert(viewer.missingReason(m,f,0).includes('未'));}
assert(Math.abs(m.faceArea.reduce((a,b)=>a+b,0)-1.9646292565715553)<1e-12);
console.log(JSON.stringify({passed:true,fields:13,computedValuesVerified:11*19752,heatFluxUnknown:true,browserInteractionTested:false}));
