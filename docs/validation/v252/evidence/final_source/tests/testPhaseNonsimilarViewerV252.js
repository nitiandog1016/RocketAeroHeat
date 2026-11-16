/* 模块职责：真实三分支的面积/温度身份、逐值回读及剖面CSV验证，无浏览器假验收。 */
'use strict';
const fs=require('fs'),path=require('path'),assert=require('assert');
const root=path.resolve(__dirname,'..'),folder=process.argv[2],viewer=require(path.join(root,'web/viewer.js'));
const profiles=require(path.join(root,'web/phase-profiles-v252.js'));
const read=n=>JSON.parse(fs.readFileSync(path.join(folder,n),'utf8')),m=read('rocket_workbench.json'),data=read('profiles.json');
assert.equal(m.faces.length,19752);assert.equal(viewer.initialField(m),'heatFlux_base');assert.equal(m.researchIdentity.physicalQualification,false);
const q=m.fields.find(f=>f.key==='heatFlux_base');let count=0,area=0,power=0;
for(let i=0;i<19752;i++){
  if(q.status[i]===1){count++;area+=m.faceArea[i];power+=q.values[i]*m.faceArea[i];}
  for(const f of m.fields){assert.equal(f.values.length,19752);assert(f.methodNames[0]);
    if(f.status[i]===0){assert.equal(f.values[i],null);assert(viewer.missingReason(m,f,i));}
    if(f.key==='wallSolved'){assert.equal(f.status[i],0);assert.equal(f.values[i],null);}
    if(f.key==='wallInput'){assert.equal(f.status[i],2);assert.equal(f.values[i],300);}
  }
}
assert.equal(count,1306);assert(Math.abs(power-m.researchIdentity.heatSummary.signedPowerW)<1e-8);
assert.equal(viewer.stateLabel(q,1),'未验证非相似工程候选');
for(const [branch,b] of Object.entries(data)){
  assert.equal(b.stations.length,branch==='fine'?61:31);
  for(const s of b.stations){assert.equal(s.y.length,101);assert(s.y.every(Number.isFinite));assert(s.y.every((v,i)=>i===0||v>s.y[i-1]));
    const csv=profiles.profileCSV(s,branch).trim().split('\n');assert.equal(csv.length,102);
    for(let i=0;i<101;i++){const row=csv[i+1].split(',');assert.equal(Number(row[2]),s.y[i]);assert.equal(Number(row[4]),s.temperature[i]);assert.equal(row[8],s.sourceChainSHA256);}
  }
}
console.log(JSON.stringify({passed:true,fields:m.fields.length,originalFaces:19752,candidateFaces:count,area,power,branches:3,browserInteractionTested:false}));
