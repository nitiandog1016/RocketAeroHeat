/* 模块职责：调用真实前端纯函数，核对区域有值/未知、全分母、方法及阶段；不控制浏览器。 */
'use strict';
const fs=require('fs'),path=require('path'),assert=require('assert');
const root=path.resolve(__dirname,'..'),v=require(path.join(root,'web/viewer.js'));
const base=path.join(root,'outputs/v242_validation');
const m=JSON.parse(fs.readFileSync(path.join(base,'workbench01/forebody_research.json'),'utf8'));
assert.strictEqual(m.faces.length,19752);assert.strictEqual(m.researchIdentity.physicalQualification,false);
assert.strictEqual(v.initialField(m),'pressure_400');
const ids=new Set(m.researchIdentity.originalFaceIDs.map(i=>i-1));assert.strictEqual(ids.size,4434);
assert(Math.abs(m.faceArea.reduce((a,b)=>a+b,0)-1.9646292565715553)<1e-12);
let compared=0;
for(const step of m.researchIdentity.successfulUpdates){
 const rows=fs.readFileSync(path.join(base,'flow_audit01','original_faces_'+String(step).padStart(6,'0')+'.csv'),'utf8').trim().split(/\r?\n/).slice(1).map(r=>r.split(',').map(Number));
 for(const [key,col] of [['pressure',2],['fluid',3],['H',4],['S',5],['Mach',9]]){
  const f=m.fields.find(x=>x.key===key+'_'+step);assert(!f.history);
  assert.deepStrictEqual(v.comparisonRange(m,f.key),v.comparisonRange(m,key+'_0'));
  assert(v.stateNote(f).includes('不是物理时间'));
  rows.forEach(r=>{const i=r[0]-1;assert.strictEqual(f.values[i],r[col]);assert.strictEqual(v.state(f,i),1);assert(v.method(f,i).includes('Euler'));compared++;});
  f.values.forEach((x,i)=>{if(!ids.has(i)){assert.strictEqual(x,null);assert.strictEqual(v.state(f,i),0);assert(v.missingReason(m,f,i).includes('区域之外'));}});
 }
}
for(const key of ['heatFlux','wallSolved']){const f=m.fields.find(x=>x.key===key);assert(f.values.every(x=>x===null));assert(f.status.every(x=>x===0));}
const given=m.fields.find(x=>x.key==='wallInput');assert(given.values.every(x=>x===300)&&given.status.every(x=>x===2));
console.log(JSON.stringify({passed:true,faces:19752,regionFaces:4434,compared,fields:m.fields.length,browserInteractionTested:false}));
