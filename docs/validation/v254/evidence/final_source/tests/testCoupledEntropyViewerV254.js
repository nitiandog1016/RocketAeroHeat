/* 模块职责：核对稀疏原面映射、未知身份和三分支138处CSV；不冒称浏览器交互测试。 */
'use strict';
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert'),root=path.resolve(__dirname,'..'),out=process.argv[2];
const api=require(path.join(root,'web/coupled-entropy-v254.js')),data=JSON.parse(fs.readFileSync(path.join(out,'coupled.json'),'utf8'));
const context={window:{}};vm.runInNewContext(fs.readFileSync(path.join(root,'outputs/v252_validation/workbench02/case-data.js'),'utf8'),context);
const model=context.window.ROCKET_CASE,faceJSON=JSON.stringify(model.faces),vertexJSON=JSON.stringify(model.vertices),oldCount=model.fields.length;
api.applyCoupledSurface(model,data);assert.equal(JSON.stringify(model.faces),faceJSON);assert.equal(JSON.stringify(model.vertices),vertexJSON);
assert.equal(model.fields.length,oldCount+8);const map=new Map(data.faceIDs.map((id,j)=>[id,j]));let checked=0;
for(const spec of data.surfaceFields){const f=model.fields.find(f=>f.key===spec.key);assert.equal(f.values.length,19752);
 for(let i=0;i<19752;i++){if(map.has(i+1)){assert.equal(f.values[i],spec.values[map.get(i+1)]);assert.equal(f.status[i],1);}else{assert.equal(f.values[i],null);assert.equal(f.status[i],0);assert(f.missingReason[i]);}checked++;}}
let inlets=0;for(const [key,b]of Object.entries(data.branches)){assert.equal(b.stations.length,46);const csv=api.integralCSV(data,key).trim().split('\n'),headers=csv[0].split(',');assert.equal(csv.length,47);
 b.stations.forEach((s,i)=>{const row=csv[i+1].split(',');for(const name of ['upstreamFaceID','downstreamFaceID','pressurePa','heatFluxWm2','edgeTemperatureK','edgeSpeedMps','massThicknessM'])assert.equal(Number(row[headers.indexOf(name)]),s[name]);
 assert.equal(row[headers.indexOf('solvedSolidWallTemperatureK')],'');assert.equal(row[headers.indexOf('profileReadyForCFD')],'false');assert.equal(row[headers.indexOf('physicalQualification')],'false');assert.equal(s.actualProfile,null);inlets++;});}
assert.equal(api.integralSegments([0,1,2,3],[0,null,-1,2]).length,2);assert.throws(()=>api.integralCSV(data,'missing'));
assert.equal(data.identity.newQualifiedAreaM2,0);assert.equal(data.identity.actualProfileSolved,false);assert.equal(data.identity.solvedSolidWallTemperatureK,null);
assert(fs.readFileSync(path.join(out,'index.html'),'utf8').includes('iframe'));assert.equal(data.faceIDs.length+data.withheldFaceIDs.length,1306);
assert.equal(data.identity.fullCapNumericalGatePassed,false);assert(data.withheldFaceIDs.length>=4);
for(const id of data.withheldFaceIDs){assert.equal(model.fields.find(f=>f.key==='entropyIntegralHeat').values[id-1],null);assert(model.fields.find(f=>f.key==='entropyIntegralHeat').missingReason[id-1].includes(data.withheldReasons[id]));}
console.log(JSON.stringify({passed:true,surfaceValuesChecked:checked,integralCSVRowsChecked:inlets,originalGeometryUnchanged:true,unknownNotZero:true,browserInteractionTested:false}));
