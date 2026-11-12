/* 模块职责：实际原面研究数据的前端解释/共同色标/未知热量与全量CSV回读。 */
'use strict';
const fs=require('fs'),path=require('path'),assert=require('assert');
const root=path.resolve(__dirname,'..'),viewer=require(path.join(root,'web/viewer.js'));
const base=path.join(root,'outputs/v230_validation');
const m=JSON.parse(fs.readFileSync(path.join(base,'workbench01/rocket_euler_research.json'),'utf8'));
assert.strictEqual(m.fields.length,11);assert.strictEqual(m.faces.length,19752);
assert.strictEqual(viewer.initialField(m),'pressure_1');assert.strictEqual(m.researchIdentity.physicalQualification,false);
assert.deepStrictEqual(m.researchIdentity.savedUpdates,[0,1]);assert.strictEqual(m.researchIdentity.completedUpdates,30);
assert(Math.abs(m.faceArea.reduce((a,b)=>a+b,0)-1.9646292565715553)<1e-12);
const cols={pressure:2,fluidTemperature:3,totalEnthalpy:4,entropy:5,velocityX:6,velocityY:7,velocityZ:8};
const frames={};for(const s of [0,1])frames[s]=fs.readFileSync(path.join(base,'audit02','original_faces_'+String(s).padStart(6,'0')+'.csv'),'utf8').trim().split(/\r?\n/).slice(1).map(r=>r.split(',').map(Number));
for(const f of m.fields.slice(0,9)){
  assert.strictEqual(f.values.length,19752);assert(f.status.every(s=>s===(f.successfulUpdates===0?2:1)));
  assert(viewer.method(f,0).includes('Euler'));assert(!f.history);
  f.values.forEach((v,i)=>assert.strictEqual(v,frames[f.successfulUpdates][i][cols[f.comparisonVariable]]));
}
for(const key of ['pressure','fluidTemperature'])assert.deepStrictEqual(viewer.comparisonRange(m,key+'_0'),viewer.comparisonRange(m,key+'_1'));
for(const f of m.fields.slice(9)){assert(f.values.every(v=>v===null)&&f.status.every(s=>s===0));assert(viewer.missingReason(m,f,0).includes('未'));}
console.log(JSON.stringify({passed:true,fields:11,faces:19752,computedValuesVerified:9*19752,heatFluxUnknown:true,browserInteractionTested:false}));
