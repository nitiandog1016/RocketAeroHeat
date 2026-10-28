/* 五步骤同窗研究宿主的真实数据绑定与事件核对；不是浏览器视觉/WebGL验收。 */
'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const folder=path.resolve(process.argv[2]||'outputs/v68_validation/workbench');
const text=fs.readFileSync(path.join(folder,'study-data.json'),'utf8'),d=JSON.parse(text);
assert.equal(fs.readFileSync(path.join(folder,'study-data.js'),'utf8'),'window.ROCKET_STUDY='+text+';\n');
const r=JSON.parse(fs.readFileSync(path.join(folder,d.resultsJSONURL),'utf8'));
if(d.wallAdjacentGasProxyIncluded!==undefined)assert.equal(d.wallAdjacentGasProxyIncluded,r.fields.some(f=>f.key==='wallAdjacentGas'));
assert.equal(d.readOnlyFrozenStudy,true);assert.equal(d.physicalSolution,false);assert.equal(d.newProductionSolvedAreaM2,0);
assert.equal(d.acceptedSteps,r.researchIdentity.acceptedSteps);assert.equal(d.residual,r.resultCondition.residualNorm);
assert.equal(d.steadyConverged,r.researchIdentity.steadyConverged);
assert.equal(d.sourceFaces,r.sourceSurfaceProjection.faces.length);assert.equal(d.displayFaces,r.faces.length);
assert.equal(d.nativeProjectionIncluded,true);assert.equal(d.steps.length,d.acceptedSteps+1);
assert.ok(fs.existsSync(path.join(folder,d.geometryURL))&&fs.existsSync(path.join(folder,d.resultsURL)));
/* 模块用途：实际导出策略须绑定结果合同与终态审计，不能由新工厂默认值代填。 */
const strategies=d.numericalStrategies;
if(strategies){
  const c=r.resultCondition.inputContract,bc=r.resultCondition.boundaryInput;
  if(strategies.jacobianDifferencePolicy!==undefined){
    assert.equal(strategies.jacobianDifferencePolicy,r.numericalSolver.jacobianDifferencePolicy);
    assert.equal(strategies.jacobianRelativeStep,r.numericalSolver.differenceStep);
    assert.equal(strategies.numericalSolver,r.numericalSolver.solver);
    assert.deepEqual(d.initialization,r.initialization);
    if(strategies.numericalSolver==='advanceAxisymmetricSteadyLocal'){
      assert.equal(c.jacobianDifferencePolicy,strategies.jacobianDifferencePolicy);
      assert.equal(c.jacobianRelativeStep,strategies.jacobianRelativeStep);
      assert.equal(d.initialization.seed.initialGuessIsBoundaryCondition,false);
      assert.equal(c.newChainIncludesPriorAcceptedSteps,false);
      if(d.initialization.seed.kind==='v73-uniform23'){
        assert.equal(d.initialization.seed.priorAcceptedSteps,23);
        assert.equal(d.initialization.allFiveComponentsUniform,false);
        assert.match(d.initialization.seed.sourceEvidenceSHA256,/^[a-f0-9]{64}$/);
      }
    }
  }
  const wall=bc.jMin.velocityGradient||'legacy',work=bc.viscousWorkVelocity||'legacy-wls';
  assert.equal(strategies.wallVelocityGradientPolicy,wall);assert.equal(strategies.viscousWorkVelocityPolicy,work);
  assert.equal(strategies.caseFactory,r.researchIdentity.caseFactory||'createRocketAxisymmetricCase');
  if(c.wallVelocityGradientPolicy)assert.equal(wall,c.wallVelocityGradientPolicy);
  if(strategies.interiorConvectiveFluxPolicy!==undefined){
    const interior=bc.interiorConvectiveFlux||'hllc';
    assert.equal(strategies.interiorConvectiveFluxPolicy,interior);
    assert.equal(strategies.interiorConvectiveFluxAuditPresent,Boolean(r.interiorConvectiveFlux));
    if(r.interiorConvectiveFlux){
      assert.equal(r.interiorConvectiveFlux.policy,interior);
      assert.equal(r.interiorConvectiveFlux.stateClipped,false);
      assert.equal(r.interiorConvectiveFlux.adaptiveFluxFallback,false);
    }
  }
  if(strategies.pressureSourceQuadraturePolicy!==undefined){
    const pressure=bc.pressureSourceQuadrature||'linear-planar';
    assert.equal(strategies.pressureSourceQuadraturePolicy,pressure);
    assert.equal(strategies.pressureSourceQuadratureAuditPresent,Boolean(r.pressureSourceQuadrature));
    if(r.pressureSourceQuadrature){
      assert.equal(r.pressureSourceQuadrature.policy,pressure);
      assert.equal(r.pressureSourceQuadrature.negativePressureClipped,false);
      assert.equal(r.pressureSourceQuadrature.adaptivePolicyFallback,false);
    }
  }
  if(c.viscousWorkVelocityPolicy)assert.equal(work,c.viscousWorkVelocityPolicy);
  assert.equal(strategies.wallVelocityGradientAuditPresent,Boolean(r.wallVelocityGradient));
  assert.equal(strategies.viscousWorkVelocityAuditPresent,Boolean(r.viscousWorkVelocity));
  assert.equal(strategies.wallVelocityGradientAuditKind,r.wallVelocityGradient?.kind||'');
  assert.equal(strategies.viscousWorkVelocityAuditKind,r.viscousWorkVelocity?.kind||'');
  assert.equal(strategies.physicalAccuracyValidated,false);assert.equal(strategies.generalSkewGridPositivityProved,false);
  assert.equal(strategies.secondOrderArithmeticReconstructionClaimed,false);
  assert.equal(strategies.gridPolicy,c.gridPolicy||'star-ray-reference');
  if(c.gridPolicy){
    assert.equal(r.researchIdentity.gridPolicy,c.gridPolicy);
    assert.deepEqual(strategies.effectiveCollarOptions,c.effectiveCollarOptions);
  }
}
const source=fs.readFileSync(path.join(folder,'rocket-study.js'),'utf8');assert.ok(!source.includes('window.open'));
function el(){return{children:[],events:{},attributes:{},hidden:false,textContent:'',
  append(...nodes){this.children.push(...nodes);},replaceChildren(...nodes){this.children=nodes;},
  addEventListener(k,fn){this.events[k]=fn;},setAttribute(k,v){this.attributes[k]=v;},getAttribute(k){return this.attributes[k];}};}
function runUI(data){
  const ids=Object.fromEntries(['stage','stepTitle','description','details','iterationPanel','viewer','footer','iterationRows','rootIterationGroup','rootIterationRows'].map(id=>[id,el()]));
  const buttons=Array.from({length:5},el),frame=ids.viewer,panels=[];
  const document={getElementById:id=>ids[id],querySelectorAll:()=>buttons,createElement:()=>el()};
  vm.runInNewContext(source,{window:{ROCKET_STUDY:data},document});
  assert.equal(ids.viewer.getAttribute('src'),data.resultsURL);assert.equal(ids.iterationRows.children.length,data.acceptedSteps+1);
  assert.equal(ids.rootIterationGroup.hidden,!data.rootIterationHistory);
  assert.equal(ids.rootIterationRows.children.length,data.rootIterationHistory?.length||0);
  for(let i=0;i<5;i++){
    buttons[i].events.click();assert.equal(ids.viewer,frame);
    assert.equal(ids.viewer.getAttribute('src'),i===4?data.resultsURL:data.geometryURL);
    assert.equal(ids.iterationPanel.hidden,i!==3);
    buttons.forEach((b,j)=>assert.equal(b.getAttribute('aria-selected'),String(j===i)));
    assert.ok(ids.details.children.length>=6);
    const rows={};for(let j=0;j<ids.details.children.length;j+=2)rows[ids.details.children[j].textContent]=ids.details.children[j+1].textContent;
    panels.push(rows);
  }
  assert.match(ids.description.textContent,/面元/);assert.match(ids.footer.textContent,/新增生产面积0/);
  return panels;
}
const panels=runUI(d);assert.ok(panels[2]['壁速度梯度']&&panels[2]['内面黏性功速度']);
if(d.wallAdjacentGasProxyIncluded)assert.match(panels[4]['温度说明'],/邻壁气体温度另列为单元代理/);
assert.ok(panels[1]['网格构造']&&panels[1]['近壁节点设置']);
if(strategies?.gridPolicy==='normal-collar')assert.match(panels[1]['近壁节点设置'],/节点偏移不等于体质心墙距/);
assert.ok(panels[3]['壁应力梯度策略']&&panels[3]['黏性功速度策略']);
assert.ok(panels[3]['压力源积分']);
assert.ok(panels[3]['内部Euler通量']);
if(strategies?.interiorConvectiveFluxPolicy?.startsWith('hlle-')){
  assert.match(panels[3]['求解方法'],/RANS \/ HLLE/);
  assert.doesNotMatch(panels[3]['求解方法'],/RANS \/ HLLC/);
}
if(strategies?.pressureSourceQuadraturePolicy==='cell-constant')assert.match(panels[3]['压力源积分'],/单元压力分片常数.*不是裁零/);
else assert.match(panels[3]['压力源积分'],/平面质心/);
assert.match(panels[3]['结果身份'],/不代表物理精度已验证/);
if(strategies?.numericalSolver==='advanceAxisymmetricSteadyLocal'){
  assert.match(panels[3]['Jacobian扰动'],/局部密度.*声学动量.*SA舍入底限/);
  if(d.initialization.seed.kind==='v73-uniform23'){
    assert.match(panels[3]['初猜来源'],/内部冻结Q23.*此前23步单列.*本链从0计步.*不是给定入口/);
    assert.ok(panels[3]['初猜来源'].includes(d.initialization.seed.sourceEvidenceSHA256));
  }
  if(d.initialization.seed.kind==='v75-spatial-root'){
    assert.match(panels[3]['初猜来源'],/内部定常Newton解.*本链从0计BE步.*不是给定入口/);
    assert.equal(d.rootIterationHistory.length,d.initialization.seed.priorNewtonIterations+1);
    assert.equal(d.rootIterationHistory.at(-1).residual,d.initialization.seed.spatialRootResidualNorm);
    assert.ok(d.rootIterationHistory.at(-1).residual<=d.steadyTolerance);
    assert.match(panels[3]['内部通量算术'],/double/);
    assert.equal(d.steadyConverged,true);
  }
}
/* 模块用途：旧v68缺字段必须仍标原方法；两项新策略可独立选择，标签不能互相冒用。 */
const legacy=JSON.parse(JSON.stringify(d));delete legacy.numericalStrategies;
const oldPanels=runUI(legacy);assert.match(oldPanels[2]['壁速度梯度'],/原WLS与位移修正/);
assert.match(oldPanels[2]['内面黏性功速度'],/原内部面WLS速度迹/);
assert.match(oldPanels[3]['压力源积分'],/平面质心/);
assert.match(oldPanels[3]['内部Euler通量'],/^HLLC/);
for(const [wall,work] of [['dirichlet-tangent','legacy-wls'],['legacy','arithmetic'],['dirichlet-tangent','arithmetic']]){
  const fixture={...legacy,numericalStrategies:{wallVelocityGradientPolicy:wall,viscousWorkVelocityPolicy:work}};
  const p=runUI(fixture);
  if(wall==='dirichlet-tangent'){
    assert.match(p[2]['壁速度梯度'],/两速度分量沿壁导数为0/);
    assert.match(p[3]['壁应力梯度策略'],/仅改选中墙面.*内场WLS、SA与Fourier保持/);
  }else assert.match(p[2]['壁速度梯度'],/原WLS与位移修正/);
  if(work==='arithmetic'){
    assert.match(p[2]['内面黏性功速度'],/两侧单元速度代理.*算术平均.*边界速度迹保持/);
    assert.match(p[3]['黏性功速度策略'],/仅改内部面黏性功.*牵引与动量通量保持/);
    assert.match(p[3]['黏性功速度策略'],/未证明.*正性或二阶精度/);
  }else assert.match(p[2]['内面黏性功速度'],/原内部面WLS速度迹/);
}
const report={status:'passed-frozen-five-step-single-frame-workbench-data-and-events',
  steps:5,acceptedEndpoints:d.acceptedSteps,sameIframeReused:true,inputViewDoesNotSolve:true,
  independentNumericalStrategyLabels:true,legacyV68MissingStrategyFieldsPreserved:true,
  strategyBindingToActualResultsChecked:Boolean(strategies),
  dataSHA256:crypto.createHash('sha256').update(text).digest('hex'),
  browserVisualOrWebGLRuntimeVerified:false,physicalSolution:false,newProductionSolvedAreaM2:0};
fs.writeFileSync(path.join(folder,'node_host_summary.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report));
