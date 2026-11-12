/* 模块职责：核对实际15字段、热覆盖、给定/未知温度与工作台纯函数。 */
'use strict';
const fs=require('fs'),path=require('path'),assert=require('assert');
const root=path.resolve(__dirname,'..'),viewer=require(path.join(root,'web/viewer.js'));
const m=JSON.parse(fs.readFileSync(path.join(process.argv[2],'rocket_workbench.json'),'utf8'));
assert.equal(m.faces.length,19752);assert.equal(m.fields.length,15);assert.equal(viewer.initialField(m),'heatFlux');
const q=m.fields.find(f=>f.key==='heatFlux');let count=0,power=0,area=0;
for(let i=0;i<m.faces.length;i++){
  if(q.status[i]===1){count++;power+=q.values[i]*m.faceArea[i];area+=m.faceArea[i];}
  for(const f of m.fields){
    assert.equal(f.values.length,19752);assert.equal(f.status.length,19752);assert(f.methodNames[0]);
    if(f.key==='wallInput'){assert.equal(f.status[i],2);assert.equal(f.values[i],300);}
    else if(f.key==='wallSolved'){assert.equal(f.status[i],0);assert.equal(f.values[i],null);}
    else if(q.status[i]===1){assert.equal(f.status[i],1);assert(Number.isFinite(f.values[i]));}
    else{assert.equal(f.status[i],0);assert.equal(f.values[i],null);assert(viewer.missingReason(m,f,i));}
  }
}
assert.equal(count,1306);assert(Math.abs(power-m.researchIdentity.heatSummary.signedPowerW)<1e-8);
assert.equal(viewer.stateLabel(q,1),'未验证鼻部热工程候选');assert.equal(m.researchIdentity.physicalQualification,false);
assert(m.fields.find(f=>f.key==='recoveryTemperature').methodNames[0].includes('sqrt(Pr'));
console.log(JSON.stringify({passed:true,fields:15,originalFaces:19752,heatCandidateFaces:count,area,power,browserInteractionTested:false}));
