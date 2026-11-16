/* 模块职责：实际138入口的缺失身份、73923剖面样本及CSV逐值核对。 */
'use strict';
const fs=require('fs'),path=require('path'),assert=require('assert'),root=path.resolve(__dirname,'..');
const api=require(path.join(root,'web/entropy-inlet-v253.js')),m=JSON.parse(fs.readFileSync(path.join(process.argv[2],'inlets.json'),'utf8'));
assert.equal(m.profileReadyForCFD,false);assert.equal(m.originalFaces,19752);let available=0,missing=0,samples=0;
for(const [key,b] of Object.entries(m.branches)){
  assert.equal(b.stations.length,46);
  for(const s of b.stations){assert.equal(s.branch,key);assert.equal(s.profileReadyForCFD,false);assert(s.sourceChainSHA256.length===64);
    if(!s.available){missing++;assert.equal(s.profile,null);assert.equal(s.matchedEdge,null);assert.throws(()=>api.inletCSV(s),/UNAVAILABLE/);continue;}
    available++;assert.equal(s.profile.y.length,601);const csv=api.inletCSV(s).trim().split('\n');assert.equal(csv.length,602);
    for(let j=0;j<601;j++){const row=csv[j+1].split(',');assert.equal(Number(row[5]),s.profile.y[j]);assert.equal(Number(row[12]),s.profile.temperature[j]);assert.equal(Number(row[13]),s.profile.density[j]);assert.equal(Number(row[14]),s.profile.rhoE[j]);assert.equal(row[17],'false');samples++;}
  }
  const seg=api.finiteSegments(b.stations.map(s=>s.phiDeg),b.stations.map(s=>s.available?s.edgeSpeedRelativeCorrection:null));assert.equal(seg.reduce((n,s)=>n+s.length,0),41);assert.equal(seg.length,2);
}
assert.equal(available,123);assert.equal(missing,15);assert.equal(samples,73923);
console.log(JSON.stringify({passed:true,available,missing,samples,browserInteractionTested:false,unknownNotZero:true}));
