/* 完整箭体真实端点的离线数据、闭合几何、源面保守投影和周向曲线独立核对。 */
'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const folder=path.resolve(process.argv[2]||'outputs/v68_validation/rocket_view');
const jsonText=fs.readFileSync(path.join(folder,'results.json'),'utf8');
const script=fs.readFileSync(path.join(folder,'case-data.js'),'utf8');
assert.equal(script,'window.ROCKET_CASE='+jsonText+';\n','离线内嵌数据必须与JSON逐字节相同');
const d=JSON.parse(jsonText),fields=Object.fromEntries(d.fields.map(f=>[f.key,f]));
const section=require(path.join(folder,'section.js')),exportsAPI=require(path.join(folder,'export.js'));
const meaning=require(path.join(folder,'viewer.js'));
const id=d.researchIdentity,a=d.surfaceIntegralAudit,b=d.wallBinding,n=d.sourceSurfaceProjection;
if(d.numericalSolver?.solver==='advanceAxisymmetricSteadyLocal'){
  const c=d.resultCondition.inputContract;
  assert.ok(['createRocketAxisymmetricLocalCase','createRocketAxisymmetricRootCase','createRocketAxisymmetricSensitivityCase','createRocketAxisymmetricWallResolutionCase','createRocketAxisymmetricContactCase','createRocketAxisymmetricContactRefinedCase'].includes(id.caseFactory));
  assert.equal(d.numericalSolver.jacobianDifferencePolicy,c.jacobianDifferencePolicy);
  assert.equal(d.numericalSolver.differenceStep,c.jacobianRelativeStep);
  assert.equal(d.initialization.seed.initialGuessIsBoundaryCondition,false);
  assert.deepEqual(d.initialization.seed,c.internalInitialGuess);
  assert.equal(c.newChainIncludesPriorAcceptedSteps,false);
  if(id.caseFactory==='createRocketAxisymmetricRootCase'){
    assert.equal(d.initialization.seed.kind,'v75-spatial-root');
    assert.equal(d.numericalSolver.initializationSolver,'solveAxisymmetricSteadyNewton');
    assert.equal(d.initialization.seed.priorAcceptedSteps,0);
    assert.equal(d.initialization.seed.priorNewtonIterations,3);
  }
}
if(d.interiorConvectiveFlux){
  const q=d.interiorConvectiveFlux,policy=d.resultCondition.boundaryInput.interiorConvectiveFlux;
  assert.equal(q.policy,policy);assert.equal(id.interiorConvectiveFluxPolicy,policy);
  assert.equal(q.stateClipped,false);assert.equal(q.adaptiveFluxFallback,false);
  if(policy.startsWith('hlle-')){
    assert.match(d.resultCondition.sourceMethod,/first-order HLLE/);
    assert.doesNotMatch(d.resultCondition.sourceMethod,/first-order HLLC/);
    for(const key of ['pressure','heatFlux','wallAdjacentGas']){
      const methods=fields[key].methodNames.join('\n');
      assert.match(methods,/first-order HLLE/);
      assert.doesNotMatch(methods,/first-order HLLC/);
      assert.ok(methods.includes('explicit interior Euler flux='+policy));
    }
  }
}
if(d.pressureSourceQuadrature){
  const q=d.pressureSourceQuadrature,policy=d.resultCondition.boundaryInput.pressureSourceQuadrature;
  assert.equal(q.policy,policy);assert.equal(id.pressureSourceQuadraturePolicy,policy);
  assert.equal(q.cellStateClipped,false);assert.equal(q.negativePressureClipped,false);assert.equal(q.adaptivePolicyFallback,false);
  if(policy==='cell-constant')assert.deepEqual(q.cellPressureProxyPa,q.sourcePressureProxyPa);
}
assert.equal(id.physicalSolution,false);assert.equal(id.eligibleForSurfaceCoverage,false);
assert.equal(id.experimentalValidation,false);assert.equal(id.newProductionSolvedAreaM2,0);
assert.equal(id.pseudoTimeIsPhysical,false);assert.equal(id.physicalTimeAdvancedS,0);
assert.ok(id.acceptedSteps>0);assert.equal(id.nativeOriginalTriangleProjectionCompleted,true);
assert.equal(id.displayGeometryIsOriginalNativeTriangles,false);
assert.equal(d.checkpointValidation.acceptedEndpointsReplayed,id.acceptedSteps);
assert.equal(d.checkpointValidation.immutableChainProvided,true);
if(!id.steadyConverged)assert.match(d.caseNote,/未达稳态/);
assert.match(d.caseNote,/新增生产求解面积0/);
const near=(x,y,rel=2e-11,abs=1e-10)=>assert.ok(Math.abs(x-y)<=abs+rel*Math.max(Math.abs(x),Math.abs(y)),`${x} != ${y}`);
const sum=x=>x.reduce((a,b)=>a+b,0);

// 从真实三角连接独立重建全部共享边；两个极点须唯一，不能加虚盖补漏洞。
const counts=new Map(),directions=new Map();let volume=0;
for(const face of d.faces){
  assert.equal(new Set(face).size,3);
  for(let k=0;k<3;k++){
    const i=face[k],j=face[(k+1)%3];assert.ok(i>=0&&i<d.vertices.length);
    const key=Math.min(i,j)+':'+Math.max(i,j);
    counts.set(key,(counts.get(key)||0)+1);directions.set(key,(directions.get(key)||0)+(i<j?1:-1));
  }
  const [v0,v1,v2]=face.map(i=>d.vertices[i]);
  const cross=[v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0]];
  volume+=sum(v0.map((x,j)=>x*cross[j]))/6;
}
assert.ok([...counts.values()].every(v=>v===2));assert.ok([...directions.values()].every(v=>v===0));
assert.equal(d.vertices.filter(v=>v[1]===0&&v[2]===0).length,2);assert.ok(volume>0);
assert.equal(a.fictitiousEndcapsAdded,false);near(sum(d.faceArea),a.displayFacetedAreaM2);

for(const key of ['heatFlux','pressure','wallInput'])assert.equal(fields[key].values.length,d.faces.length);
// 模块用途：初期全负或后期全正都可是真实signed结果；逐面与源值比较符号，禁止abs热流。
assert.ok(fields.heatFlux.values.every((q,i)=>(q<0)===(b.sourceHeatFluxIntoWallWm2[i]<0)&&(q>0)===(b.sourceHeatFluxIntoWallWm2[i]>0)));
assert.ok(fields.heatFlux.status.every(s=>s===1));assert.ok(fields.pressure.status.every(s=>s===1));
assert.ok(fields.wallInput.status.every(s=>s===2));
for(let i=0;i<d.faces.length;i++){
  near(fields.heatFlux.values[i]*1000,b.sourceHeatFluxIntoWallWm2[i]);
  near(fields.pressure.values[i]*1000,b.sourcePressurePa[i]);
  assert.equal(fields.wallInput.values[i],b.sourceTemperatureInputK[i]);
  assert.equal(meaning.state(fields.wallInput,i),2);
}
for(const key of ['wallSolved','edgeTemperature','recoveryTemperature']){
  const f=fields[key];assert.ok(f.values.every(v=>v===null));assert.ok(f.status.every(s=>s===0));
  assert.ok(f.methodNames.length&&f.scope.length);assert.ok(meaning.missingReason(d,f,0).length>0);
}
// 模块用途：邻壁气体代理另列；三种未知温度仍保持null，不能借代理填充。
if(fields.wallAdjacentGas){
  assert.ok(fields.wallAdjacentGas.status.every(s=>s===1));
  assert.match(fields.wallAdjacentGas.methodNames[0],/不是壁温、外缘温度或恢复温度/);
  fields.wallAdjacentGas.values.forEach((t,i)=>near(t,b.sourceWallAdjacentGasProxyK[i]));
  assert.ok(n.wallAdjacentGasProxyStatus.every(s=>s===1));
  assert.equal(n.gasTemperatureProjectionProvesEnergyConservation,false);
}
near(sum(d.faceArea.map((area,i)=>area*fields.heatFlux.values[i]*1000)),a.displayFacetedHeatIntegralW);

// 独立将JSON中的三元组还原为逐原始面平均值与面积账，不使用MATLAB稀疏运算。
assert.equal(n.indexBase,0);assert.equal(n.originalSourceFacesPreserved,true);assert.equal(n.candidateMappingOnly,false);
assert.equal(n.pressureForceOrMomentConservationEstablished,false);
const byWall=new Map();b.sourceFaceID.forEach((f,i)=>byWall.set(f,{p:b.sourcePressurePa[i],q:b.sourceHeatFluxIntoWallWm2[i],tw:b.sourceTemperatureInputK[i],gas:b.sourceWallAdjacentGasProxyK?.[i]}));
const p=Array(n.faces.length).fill(0),q=p.slice(),tw=p.slice(),gas=p.slice(),rowSum=p.slice(),rowArea=p.slice();
for(let k=0;k<n.weight.length;k++){
  const i=n.sourceFaceIndex[k],j=n.wallSegmentIndex[k],w=n.weight[k],wall=byWall.get(n.wallBoundaryFaceID[j]);
  assert.ok(w>0&&n.overlapAreaM2[k]>0&&wall);rowSum[i]+=w;rowArea[i]+=n.overlapAreaM2[k];
  p[i]+=w*wall.p;q[i]+=w*wall.q;tw[i]+=w*wall.tw;
  if(fields.wallAdjacentGas)gas[i]+=w*wall.gas;
}
for(let i=0;i<n.faces.length;i++){
  near(rowSum[i],1);near(rowArea[i],n.sourceCurvedAreaM2[i]);near(p[i],n.pressurePa[i]);
  near(q[i],n.heatFluxIntoWallWm2[i]);near(tw[i],n.wallInputK[i]);
  if(fields.wallAdjacentGas)near(gas[i],n.wallAdjacentGasProxyK[i]);
  assert.equal(n.pressureStatus[i],1);assert.equal(n.heatFluxStatus[i],1);assert.equal(n.wallInputStatus[i],2);
  assert.equal(n.wallSolvedK[i],null);assert.equal(n.edgeTemperatureK[i],null);assert.equal(n.recoveryTemperatureK[i],null);
}
near(sum(n.sourceCurvedAreaM2.map((area,i)=>area*q[i])),a.revolvedWallHeatW);
near(sum(n.sourceFacetAreaM2.map((area,i)=>area*q[i])),n.originalFacetHeatIntegralW);

// 任意X处只读取实际三角面交线及本面常值；未解温度不能在曲线或CSV中变零。
const xmin=Math.min(...d.vertices.map(v=>v[0])),xmax=Math.max(...d.vertices.map(v=>v[0]));
const sections=[];
for(const fraction of [.08,.27,.55,.86,.995]){
  const x=xmin+(xmax-xmin)*fraction,g=section.geometry(d.vertices,d.faces,x);
  assert.equal(g.components.length,1);assert.equal(g.components[0].closed,true);assert.ok(g.segments.length>0);
  const samples=section.sample(g,fields.heatFlux);
  for(const s of samples){assert.equal(s.value,fields.heatFlux.values[s.face]);assert.equal(s.status,1);}
  const unknown=section.sample(g,fields.wallSolved);assert.ok(unknown.every(s=>s.value===null&&s.status===0));
  if(fields.wallAdjacentGas){
    const temperature=section.sample(g,fields.wallAdjacentGas);
    temperature.forEach(s=>assert.equal(s.value,fields.wallAdjacentGas.values[s.face]));
    assert.match(section.csv(g,fields.wallAdjacentGas,'伪时间研究',d.caseNote),/不是壁温、外缘温度或恢复温度/);
  }
  const csv=section.csv(g,fields.heatFlux,'伪时间研究',d.caseNote);
  assert.match(csv,/SA/);assert.match(csv,/kW\/m/);assert.ok(!csv.includes('undefined'));
  const missingCSV=section.csv(g,fields.wallSolved,'伪时间研究',d.caseNote);assert.match(missingCSV,/未求解/);
  sections.push({x,segments:g.segments.length,perimeter:g.components[0].perimeter});
}
const csv=exportsAPI.csv(d,fields.heatFlux,{caseNote:d.caseNote});
assert.match(csv,/新增生产求解面积0/);assert.match(csv,/Fourier/);
fs.writeFileSync(path.join(folder,'verified-heat-flux.csv'),csv);
const report={status:'passed-node-data-geometry-projection-section-and-csv-contract',displayFaces:d.faces.length,
  nativeFaces:n.faces.length,acceptedEndpoints:id.acceptedSteps,sections,signedHeatPreserved:true,
  jsonSHA256:crypto.createHash('sha256').update(jsonText).digest('hex'),
  browserVisualOrWebGLRuntimeVerified:false,physicalSolution:false,experimentalValidation:false,newProductionSolvedAreaM2:0};
fs.writeFileSync(path.join(folder,'node_view_summary.json'),JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify(report));
