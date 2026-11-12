/* 模块职责：真实原面九字段与已有工作台纯函数核对，不调用浏览器。 */
'use strict';
const fs=require('fs'),path=require('path'),assert=require('assert');
const root=path.resolve(__dirname,'..'),viewer=require(path.join(root,'web/viewer.js'));
const folder=process.argv[2],m=JSON.parse(fs.readFileSync(path.join(folder,'rocket_workbench.json'),'utf8'));
assert.equal(m.faces.length,19752);assert.equal(m.fields.length,9);assert.equal(viewer.initialField(m),'pressure');
assert.equal(m.researchIdentity.productionQualified,false);
const p=m.fields.find(f=>f.key==='pressure'),valid=p.status.map(s=>s===1);
assert(valid.some(Boolean));let candidateArea=0;
for(let i=0;i<19752;i++){
  if(valid[i])candidateArea+=m.faceArea[i];
  for(const f of m.fields){
    assert.equal(f.values.length,19752);assert.equal(f.methodIndex.length,19752);
    if(f.key==='wallInput'){assert.equal(f.status[i],2);assert.equal(f.values[i],300);}
    else if(['heatFlux','wallSolved'].includes(f.key)){
      assert.equal(f.status[i],0);assert.equal(f.values[i],null);assert(viewer.missingReason(m,f,i));
    }else if(valid[i]){assert.equal(f.status[i],1);assert(Number.isFinite(f.values[i]));}
    else{assert.equal(f.status[i],0);assert.equal(f.values[i],null);}
  }
}
assert.equal(viewer.stateLabel(p,1),'未验证鼻部工程候选');
assert(Math.abs(m.faceArea.reduce((a,b)=>a+b,0)-m.totalArea)<1e-12);
console.log(JSON.stringify({passed:true,fields:9,faces:19752,candidateAreaM2:candidateArea,
  unknownHeatFaces:19752,browserInteractionTested:false}));
