/* 模块职责：原生目标逐面映射、逐变量公式、未知值及138入口CSV精确回读。 */
'use strict';
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert'),root=path.resolve(__dirname,'..'),out=process.argv[2];
const read=p=>JSON.parse(fs.readFileSync(p,'utf8')),api=require(path.join(root,'web/regularized-entropy-v255.js')),data=read(path.join(out,'coupled.json'));
const audit=read(path.join(out,'../audit01/report.json')),base=read(path.join(out,'..',audit.matlabFolder,'base.json')),source=read(path.join(out,'../prepare01/source.json'));
const context={window:{}};vm.runInNewContext(fs.readFileSync(path.join(root,'outputs/v252_validation/workbench02/case-data.js'),'utf8'),context);
const model=context.window.ROCKET_CASE,faces=JSON.stringify(model.faces),vertices=JSON.stringify(model.vertices);api.applyCoupledSurface(model,data);
assert.equal(JSON.stringify(model.faces),faces);assert.equal(JSON.stringify(model.vertices),vertices);assert.equal(model.faces.length,19752);
const map=new Map(data.faceIDs.map((id,i)=>[id,i]));let count=0;
for(const spec of data.surfaceFields){const f=model.fields.find(f=>f.key===spec.key);assert.equal(f.methodNames[0],spec.method);assert.equal(f.values.length,19752);
 for(let i=0;i<19752;i++){if(map.has(i+1)){assert.equal(f.values[i],spec.values[map.get(i+1)]);assert.equal(f.status[i],1);}else{assert.equal(f.values[i],null);assert.equal(f.status[i],0);assert(f.missingReason[i]);}count++;}}
assert(model.fields.find(f=>f.key==='entropyIntegralPressure').methodNames[0].includes('修正牛顿'));
for(const [key,field] of [['entropyIntegralHeat','heatFluxWm2'],['entropyIntegralMassThickness','massThicknessM'],['entropyIntegralTe','edgeTemperatureK'],['entropyIntegralUe','edgeSpeedMps']]){
 const f=model.fields.find(f=>f.key===key);for(const id of data.faceIDs){const i=source.surfaceFaceIDs.indexOf(id),target=source.surfaceTargetIndicesZeroBased[i],j=base.targetIndices[target]-1;assert.equal(f.values[id-1],base[field][j]);}}
assert.equal(data.faceIDs.length,1306);assert.equal(data.withheldFaceIDs.length,0);assert.equal(data.identity.recoveredPriorFaceIDs.length,31);
for(const id of data.identity.recoveredPriorFaceIDs)assert(Number.isFinite(model.fields.find(f=>f.key==='entropyIntegralHeat').values[id-1]));
let rows=0;for(const [key,b] of Object.entries(data.branches)){const csv=api.integralCSV(data,key).trim().split('\n'),headers=csv[0].split(',');assert.equal(csv.length,47);
 b.stations.forEach((s,i)=>{const r=csv[i+1].split(',');for(const k of ['upstreamFaceID','downstreamFaceID','pressurePa','heatFluxWm2','massThicknessM'])assert.equal(Number(r[headers.indexOf(k)]),s[k]);assert.equal(r[headers.indexOf('solvedSolidWallTemperatureK')],'');assert.equal(r[headers.indexOf('profileReadyForCFD')],'false');assert.equal(s.actualProfile,null);rows++;});}
assert.equal(data.identity.newQualifiedAreaM2,0);assert.equal(data.identity.actualProfileSolved,false);assert.equal(data.identity.positiveBelowSeedSolved,false);
assert.equal(data.legacyAlong.heatFluxWm2.length,data.legacyAlong.angleDeg.length);assert(fs.readFileSync(path.join(out,'index.html'),'utf8').includes('iframe'));
console.log(JSON.stringify({passed:true,exactSurfaceSlots:count,exactCSVRows:rows,directNativeTargets:true,originalMeshUnchanged:true,explicitPerFieldMethod:true,browserInteractionTested:false}));
