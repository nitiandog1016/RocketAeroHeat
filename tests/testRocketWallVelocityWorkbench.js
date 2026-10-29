/* 真实原面静止壁梯度：数值身份、给定/未知、CSV/截线与同宿主五步骤事件。 */
'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const folder=path.resolve(process.argv[2]),raw=fs.readFileSync(path.join(folder,'wall-velocity-study-data.json'),'utf8'),d=JSON.parse(raw);
const meaning=require('../web/viewer.js'),exp=require('../web/export.js'),section=require('../web/section.js');
assert.equal(d.cases.length,3);assert.equal(d.sourceFaces,19752);assert.equal(d.cellCount,217272);assert.equal(d.productionQualified,false);
assert.deepEqual(d.convergenceTargets,{maximumResidualPerSecond:2e-6,steadyBoundaryFluxBalance:5e-10,spatialConservationRelative:1e-9});
assert.ok(fs.readFileSync(path.join(folder,'index.html'),'utf8').includes('<th>定常边界不平衡</th>'));
const linearColumnCount=fs.readFileSync(path.join(folder,'index.html'),'utf8').includes('<th>粗层代数残差</th>')?9:8;
assert.equal(d.cases[0].summary.residualEvaluations,0);assert.equal(d.cases[0].summary.sourceAcceptedIterations,3);
assert.equal(fs.readFileSync(path.join(folder,'wall-velocity-study-data.js'),'utf8'),'window.ROCKET_WALL_VELOCITY_STUDY='+raw+';\n');
const snapshots=[];let CSVrows=0,sectionSamples=0;
for(const c of d.cases){
 const data=JSON.parse(fs.readFileSync(path.join(folder,c.resultsJSONURL),'utf8'));snapshots.push(data);
 assert.deepEqual(data.vertices,snapshots[0].vertices);assert.deepEqual(data.faces,snapshots[0].faces);assert.deepEqual(data.faceArea,snapshots[0].faceArea);
 assert.equal(data.faces.length,19752);assert.equal(data.fields.length,10);assert.deepEqual(data.comparisonColorRanges,d.commonRanges);
 assert.equal(data.researchIdentity.wallVelocityGradient,c.wallVelocityGradient);assert.equal(data.researchIdentity.sourceHistoryAppended,false);
 assert.equal(data.researchIdentity.productionQualified,false);assert.equal(data.researchIdentity.stationaryWallExternalPowerW,0);
 const fields=Object.fromEntries(data.fields.map(f=>[f.key,f]));
 assert.equal(meaning.coverageLabels(data,fields.heatFlux)[0],'研究数值面积');
 assert.match(meaning.coverageNote(data,fields.heatFlux),/不代表正式合格求解面积/);
 for(const key of ['edgeTemperature','recoveryTemperature','wallSolved']){assert.ok(fields[key].values.every(x=>x===null));assert.ok(fields[key].status.every(x=>x===0));assert.ok(fields[key].missingReason.every(x=>x.length>5));}
 assert.ok(fields.wallInput.status.every(x=>x===2));assert.ok(fields.wallInput.values.every(x=>Math.abs(x-300)<1e-10));
 for(const key of ['pressure','heatFlux','wallAdjacentTemperature','temperatureRate','wallMechanicalExchange']){
  const f=fields[key];assert.ok(f.values.every(Number.isFinite));assert.ok(f.status.every(x=>x===1));assert.equal(meaning.stateLabel(f,1),'未收敛研究数值');
 }
 assert.ok(fields.wallMechanicalExchange.methodNames[0].includes('静止壁外功仍为0'));
 const negative=fields.wallMechanicalExchange.values.filter(x=>x<0).length;assert.equal(negative,c.mechanicalCoolingWallCount);
 if(c.wallVelocityGradient==='dirichlet-tangent')assert.equal(negative,0);
 let power=0;data.faceArea.forEach((a,i)=>power+=a*fields.heatFlux.values[i]*1000);assert.ok(Math.abs(power-c.signedHeatPowerW)/Math.max(1,Math.abs(power))<3e-13);
 for(const f of data.fields){assert.equal(exp.csv(data,f,{caseNote:data.caseNote}).trim().split(/\r?\n/).length,data.faces.length+1);CSVrows+=data.faces.length;}
 for(const x of [-2.5,-.5,-.2,0]){
  const geometry=section.geometry(data.vertices,data.faces,x);assert.ok(geometry.segments.length>0);
  for(const key of ['pressure','heatFlux','wallAdjacentTemperature','temperatureRate','wallMechanicalExchange','edgeTemperature']){
   const samples=section.sample(geometry,fields[key]);samples.forEach(s=>{assert.equal(s.value,fields[key].values[s.face]);assert.equal(s.status,fields[key].status[s.face]);});sectionSamples+=samples.length;
  }
 }
}
const fields=snapshots.map(s=>Object.fromEntries(s.fields.map(f=>[f.key,f])));
if(d.cases[2].summary.openBoundaryDiscretization==='exact-Riemann-reservoir'){
 assert.equal(snapshots[2].researchIdentity.openBoundaryDiscretization,'exact-Riemann-reservoir');
 assert.equal(snapshots[2].researchIdentity.priorProblemHistoryAppended,false);
 assert.ok(d.cases[2].summary.priorAcceptedNonlinearIterations>=0);
 assert.equal(d.cases[2].summary.totalNewEquationAcceptedIterations,d.cases[2].summary.priorAcceptedNonlinearIterations+d.cases[2].summary.acceptedNonlinearIterations);
 assert.equal(d.cases[2].history.iteration[0],0);
 for(const key of ['pressure','heatFlux','wallAdjacentTemperature','temperatureRate'])assert.ok(fields[2][key].methodNames.every(x=>x.includes('开放储库为精确Riemann')));
}
for(const key of ['pressure','heatFlux','wallAdjacentTemperature','wallInput'])assert.deepEqual(fields[0][key].values,fields[1][key].values);
assert.notDeepEqual(fields[0].temperatureRate.values,fields[1].temperatureRate.values);assert.notDeepEqual(fields[0].wallMechanicalExchange.values,fields[1].wallMechanicalExchange.values);
assert.ok(d.cases[0].mechanicalCoolingWallCount>0);assert.equal(d.cases[1].mechanicalCoolingWallCount,0);assert.equal(d.cases[2].history.iteration[0],0);
class E{constructor(tag){this.tag=tag;this.children=[];this.attrs={};this.events={};this.value='';this.textContent='';}append(...e){this.children.push(...e);}replaceChildren(...e){this.children=[...e];}setAttribute(k,v){this.attrs[k]=String(v);}getAttribute(k){return this.attrs[k];}addEventListener(k,f){this.events[k]=f;}}
const ids=Object.fromEntries(['caseChoice','title','description','details','viewer','historyRows','linearRows','trialRows','diagnosticNote','layerRows','layerNote','wallNote','stateNote','linearResearchRows','linearResearchNote'].map(x=>[x,new E(x)])),buttons=Array.from({length:5},()=>new E('button')),viewer=ids.viewer;
const document={getElementById:x=>{assert.ok(ids[x],x);return ids[x];},querySelectorAll:q=>{assert.equal(q,'[data-step]');return buttons;},createElement:t=>new E(t)};
const program=process.argv[3]?path.resolve(process.argv[3]):path.join(folder,'rocket-wall-velocity.js');
vm.runInNewContext(fs.readFileSync(program,'utf8'),{window:{ROCKET_WALL_VELOCITY_STUDY:d},document});let events=0,linearDiagnosticRows=0,candidateDiagnosticRows=0;
for(let j=0;j<3;j++){
 ids.caseChoice.value=String(j);ids.caseChoice.events.change();const h=d.cases[j].history;assert.equal(ids.historyRows.children.length,h&&Array.isArray(h.iteration)?h.iteration.length:0);
 const trend=d.cases[j].solverDiagnostics&&d.cases[j].solverDiagnostics.segmentResidualTrend;
 const coordinates=d.cases[j].solverDiagnostics&&d.cases[j].solverDiagnostics.trialStateCoordinates;
 const preparation=d.cases[j].solverDiagnostics&&d.cases[j].solverDiagnostics.directionPreparation;
 if(preparation){
  assert.ok(!preparation.hardUpdateLimitClaimed);
  assert.ok(preparation.usedTotalDirectionExactlyVerified||(preparation.proposedTotalDirectionBoundToCurrentState&&preparation.displayedStateUnchanged));
  assert.equal(preparation.internalCTargetMet,preparation.internalCRelativeResidual<=preparation.internalCTarget);
  if(preparation.usedTotalDirectionExactlyVerified)assert.ok(preparation.originalNSDirectionDefect<=preparation.originalNSDirectionTarget);
  buttons[3].events.click();assert.ok(ids.details.children.some(e=>e.textContent===preparation.methodLabel));
  assert.ok(ids.details.children.some(e=>e.textContent.includes('不保证硬幅度上限')));
  if(preparation.kind==='temperature-comparison-QR-LSQR-v1'){
   assert.equal(preparation.physicalTemperatureBoundaryAdded,false);assert.equal(preparation.nonlinearTemperatureBoundGuaranteed,false);
   const study=d.temperatureProjectionResearch;assert.ok(study&&study.physicalSolution===false&&study.newProductionSolvedAreaM2===0);
   const flatten=e=>[e,...(e.children||[]).flatMap(flatten)],nodes=ids.details.children.flatMap(flatten);
   assert.ok(nodes.some(e=>e.tag==='details'));assert.ok(nodes.some(e=>e.textContent===study.constraintLabel));
   const links=nodes.filter(e=>e.tag==='a');const link=links.find(e=>e.attrs.href===study.reportURL);assert.ok(link);
   let prevented=false;link.events.click({preventDefault(){prevented=true;}});assert.ok(prevented);assert.equal(ids.viewer.attrs.src,study.reportURL);
   buttons[4].events.click();assert.equal(ids.viewer.attrs.src,d.cases[j].resultsURL);
  }
  if(preparation.kind==='current-band-active-temperature-v1'){
   assert.equal(preparation.oldActiveSetInherited,false);assert.equal(preparation.physicalTemperatureBoundaryAdded,false);
   assert.equal(preparation.nonlinearTemperatureBoundGuaranteed,false);assert.equal(preparation.constrainedOptimumClaimed,false);
   if(preparation.lookaheadSelection){
    const q=preparation.lookaheadSelection,asArray=v=>Array.isArray(v)?v:v==null?[]:[v];
    const originals=asArray(q.originalBandFaceIDs),required=asArray(q.requiredOriginalStateFaceIDs);
    if(q.selectionPolicy==='original-and-predicted-required-faces'){
     assert.deepEqual(asArray(preparation.bandFaceIDs),[...new Set([...originals,...required])].sort((a,b)=>a-b));
     assert.equal(q.constraintRankUsedToSelectFaces,false);assert.equal(q.allPredictedViolatingCellsRepresented,true);
    }
    if(q.selectionPolicy==='original-and-predicted-cell-and-interface-faces'){
     const preview=q.interfacePreview,interfaces=asArray(q.predictedViolatingInterfaceFaceIDs);
     assert.deepEqual(asArray(preparation.bandFaceIDs),[...new Set([...originals,...required,...interfaces])].sort((a,b)=>a-b));
     const assessed=asArray(preview.assessedFaceIDs),blocked=asArray(preview.blockedFaceIDs);
     assert.equal(new Set([...assessed,...blocked]).size,preview.completeInteriorFaceCount);
     assert.equal(assessed.length+blocked.length,preview.completeInteriorFaceCount);
     assert.equal(preview.allInteriorFacesAssessed,blocked.length===0);
     assert.equal(preview.temperatureDomainRelaxed,false);
     assert.ok(ids.details.children.some(e=>e.textContent.includes('因单元超域阻挡 '+blocked.length+' 面')));
     assert.ok(ids.details.children.some(e=>e.textContent.includes('已检查 '+assessed.length+' 面')));
     const diagnostic=fs.readFileSync(path.join(folder,d.temperatureActiveResearch.reportURL),'utf8');
     assert.ok(diagnostic.includes('预测越温界面')&&diagnostic.includes('阻挡面'));
     // 尚未执行的预览不能显示为“已检查零面”，保留缺失与零的区别。
     q.interfacePreview=[];buttons[3].events.click();
     assert.ok(ids.details.children.some(e=>e.textContent.includes('未执行：目标候选')));
     q.interfacePreview=preview;buttons[3].events.click();
    }
    if(q.selectionPolicy==='original-and-iteratively-predicted-cell-and-interface-faces'){
     const r=preparation.repreview,h=Array.isArray(r.history)?r.history:[r.history];
     assert.equal(h.length,r.previewCount);assert.ok(r.previewCount<=4&&preparation.activeSolveCount<=3);
     assert.deepEqual(asArray(preparation.bandFaceIDs),[...new Set([...originals,...required,...asArray(q.predictedViolatingInterfaceFaceIDs)])].sort((a,b)=>a-b));
     assert.ok(ids.details.children.some(e=>e.textContent.includes('仅指定步长的全部单元和内部Riemann温域')));
     const page=fs.readFileSync(path.join(folder,d.temperatureActiveResearch.reportURL),'utf8');
     assert.ok(page.includes('投影后逐轮温域复查')&&page.includes('累计必要面'));
    }
    assert.ok(ids.details.children.some(e=>e.textContent.includes('原10 K带内 '+originals.length+' 面')));
    assert.ok(ids.details.children.some(e=>e.textContent.includes('实际选择 '+asArray(preparation.bandFaceIDs).length+' 面')));
    if(preparation.freeCheckpointReplay){
     assert.equal(preparation.newFreeNativeIterations,0);assert.equal(preparation.reusedFreeNativeIterations,128);
     assert.ok(ids.details.children.some(e=>e.textContent.includes('复用自由方向 128 步')));
    }
   }
   if(preparation.projectionRepresentation){
    assert.equal(preparation.projectionRepresentation.length,preparation.activeSolveCount);
    preparation.projectionRepresentation.forEach(p=>{
     assert.ok(p.originalRows>=p.distinctConditions&&p.distinctConditions>0);
     assert.equal(p.allOriginalRowsExactlyRepresented,true);assert.equal(p.rankTruncationUsed,false);
     assert.ok(p.basisArrayBytes>0&&p.equivalentFullBasisArrayBytes>=p.basisArrayBytes);
     if(p.compactStorage)assert.ok(p.storedCoordinates>0);
    });
    if(preparation.projectionRepresentation.length){
     const last=preparation.projectionRepresentation.at(-1),faces=preparation.activeFaceIDs;
     assert.equal(last.originalRows,Array.isArray(faces)?faces.length:faces==null?0:1);
    }
    const page=fs.readFileSync(path.join(folder,d.temperatureActiveResearch.reportURL),'utf8');
    assert.ok(page.includes('方向约束与表示核对')&&page.includes('等价完整基 KiB'));
   }
   const study=d.temperatureActiveResearch;assert.ok(study&&study.physicalSolution===false&&study.newProductionSolvedAreaM2===0);
   const flatten=e=>[e,...(e.children||[]).flatMap(flatten)],nodes=ids.details.children.flatMap(flatten);
   assert.ok(nodes.some(e=>e.tag==='details'));assert.ok(nodes.some(e=>e.textContent===study.constraintLabel));
   const link=nodes.find(e=>e.tag==='a'&&e.attrs.href===study.reportURL);assert.ok(link);
   let prevented=false;link.events.click({preventDefault(){prevented=true;}});assert.ok(prevented);assert.equal(ids.viewer.attrs.src,study.reportURL);
   assert.ok(fs.existsSync(path.join(folder,study.reportURL)));buttons[4].events.click();assert.equal(ids.viewer.attrs.src,d.cases[j].resultsURL);
  }
 }
 if(coordinates){
  assert.ok(['conservative','primitive-tangent'].includes(coordinates));buttons[3].events.click();
  assert.ok(ids.details.children.some(e=>e.textContent.includes(coordinates==='primitive-tangent'?'原始变量切向更新后恢复守恒量':'守恒变量直线更新')));
 }
 if(trend){
  assert.equal(trend.rmsDecreased,trend.finalRmsMerit<trend.initialRmsMerit);
  assert.equal(trend.maximumIncreased,trend.finalMaximumResidual>trend.initialMaximumResidual);
  assert.equal(trend.finalMaximumResidual,d.cases[j].summary.maximumResidual);
  assert.equal(trend.finalRmsMerit,d.cases[j].summary.rmsMerit);
  if(trend.rmsDecreased&&trend.maximumIncreased){
   assert.ok(ids.stateNote.textContent.includes('RMS下降，但最大局部残差升高'));
   // 计算和查看结果两步都须显示，不能只藏在展开诊断里。
   for(const step of [3,4]){buttons[step].events.click();assert.ok(ids.details.children.some(e=>e.textContent.includes('一次接受不能认定稳定收敛')));}
  }
 }
 if(d.linearResearch&&j===2){
  assert.equal(d.linearResearch.sourceCase,d.cases[j].solverDiagnostics.sourceFolder);
  assert.equal(d.linearResearch.physicalStateChanged,false);assert.equal(d.linearResearch.productionQualified,false);
  assert.equal(ids.linearResearchRows.children.length,d.linearResearch.studies.length);assert.ok(ids.linearResearchNote.textContent.includes('未更新当前云图'));
  ids.linearResearchRows.children.forEach((tr,i)=>{const a=d.linearResearch.studies[i];assert.equal(tr.children.length,7);assert.equal(tr.children[0].textContent,a.label);assert.equal(tr.children[4].textContent,a.displayState||'拒绝');assert.equal(tr.children[5].textContent,'0');assert.equal(tr.children[6].textContent,a.note||'内部方向检查');
   ['innerRelativeCResidual','basisOrthogonalityError','modelRelationRelativeDifference'].forEach((key,k)=>{const lsqr=[d.linearResearch.fullFLeastSquaresDiagnostics,d.linearResearch.fullFImageScaleDiagnostics,d.linearResearch.parentLeastSquaresDiagnostics].some(q=>q&&a.caseID===q.caseID);if(a[key]===null)assert.equal(tr.children[k+1].textContent,lsqr&&k>0?'不适用（LSQR）':'未知');else assert.ok(Math.abs(Number(tr.children[k+1].textContent)-a[key])/Math.max(1e-300,Math.abs(a[key]))<1e-5);});
   const parent=d.linearResearch.flexibleReductionDiagnostics;const internalC=parent&&a.caseID===parent.caseID?Boolean(parent.internalCTargetMet):false;
   assert.equal(a.candidateUsable,internalC);assert.equal(a.orthogonalityThreshold,1e-8);assert.equal(a.modelRelationThreshold,1e-8);
  });
 }else assert.equal(ids.linearResearchRows.children.length,0);
 if(h&&Array.isArray(h.iteration))ids.historyRows.children.forEach((tr,i)=>{
  assert.equal(tr.children.length,8);assert.ok(Math.abs(Number(tr.children[6].textContent)-h.steadyBoundaryFluxBalance[i])/Math.max(1e-30,Math.abs(h.steadyBoundaryFluxBalance[i]))<1e-5);
  assert.ok(Number(tr.children[5].textContent)<d.convergenceTargets.spatialConservationRelative);
 });
 for(let step=0;step<5;step++){buttons[step].events.click();events++;assert.equal(ids.viewer,viewer);assert.equal(viewer.attrs.src,step===4?d.cases[j].resultsURL:step===1?d.domainURL:d.geometryURL);assert.equal(buttons[step].attrs['aria-selected'],'true');
  if(step===3){assert.ok(ids.details.children.some(e=>e.textContent==='定常边界不平衡'));assert.ok(ids.details.children.some(e=>e.textContent.includes('目标≤5e−10')));
   if(d.cases[j].summary.stopReason==='bounded-line-search-and-pseudostep-retries-exhausted')assert.ok(ids.details.children.some(e=>e.textContent==='有限线搜索与伪步重试已用尽（未收敛）'));
   if(d.cases[j].summary.stopReason==='steady-newton-linear-target-not-met')assert.ok(ids.details.children.some(e=>e.textContent==='直接稳态线性目标未达到（未提交更新）'));
  }
 }
 assert.ok(ids.wallNote.textContent.includes('固定初猜'));assert.ok(ids.stateNote.textContent.includes('不是静止壁外功'));
 const diag=d.cases[j].solverDiagnostics,rows=x=>Array.isArray(x)?x:x?[x]:[];
 if(diag&&diag.localSegmentOnly){
  if(Array.isArray(diag.linearRemainingLayerFractions)&&diag.linearRemainingLayerFractions.length){
   assert.equal(diag.linearRemainingLayerFractions.length,d.layerCount);assert.equal(diag.linearRemainingEquationFractions.length,5);
   assert.ok(Math.abs(diag.linearRemainingLayerFractions.reduce((a,b)=>a+b,0)-1)<1e-12);
   assert.equal(ids.layerRows.children.length,d.layerCount);assert.ok(ids.layerNote.textContent.includes('不是压力、热流误差或物理解'));
   ids.layerRows.children.forEach((tr,i)=>{assert.equal(tr.children[0].textContent,String(i+1));assert.ok(Math.abs(Number(tr.children[1].textContent)/100-diag.linearRemainingLayerFractions[i])<1e-6);});
  }
  const linear=rows(diag.linearAttempts),trials=rows(diag.nonlinearTrials);linearDiagnosticRows+=linear.length;candidateDiagnosticRows+=trials.length;
  assert.equal(ids.linearRows.children.length,linear.length);assert.equal(ids.trialRows.children.length,trials.length);
  assert.ok(ids.diagnosticNote.textContent.includes(diag.sourceFolder));assert.ok(ids.diagnosticNote.textContent.includes('完整NS线性缺陷'));
  const total=d.cases[j].summary.totalNewEquationResidualEvaluations??d.cases[j].summary.residualEvaluations;
  assert.equal(diag.residualEvaluationsAfterLastAcceptedHistory,total-h.residualEvaluations[h.residualEvaluations.length-1]);
  assert.ok(ids.diagnosticNote.textContent.includes('另有 '+diag.residualEvaluationsAfterLastAcceptedHistory+' 次NS调用'));
  if(diag.normalEulerLinearization){
   assert.ok(ids.diagnosticNote.textContent.includes(diag.normalEulerLinearization==='exact-branch'?'精确当前波系分支导数':'冻结Rusanov导数'));
   assert.ok(ids.diagnosticNote.textContent.includes('仅为预条件，正式NS不改'));
  }
  if(diag.linearEquation==='steady-newton'){assert.equal(diag.preconditionerTimeStepOnly,true);assert.ok(ids.diagnosticNote.textContent.includes('(-J)δ=R'));}
  else assert.ok(ids.diagnosticNote.textContent.includes('(I/τ-J)δ=R'));
  if(['lateral-colored-line','two-level-minres'].includes(diag.rightPreconditioner)){
   assert.ok(ids.diagnosticNote.textContent.includes('固定 '+diag.graphSweeps+' 次对称着色法向线扫描'));
   assert.ok(ids.diagnosticNote.textContent.includes('不是精确图解'));
  }
  if(diag.rightPreconditioner==='two-level-minres'){
   assert.ok(diag.coarseUnknowns>0&&diag.coarseImageBytes>0);
   assert.ok(ids.diagnosticNote.textContent.includes('粗层代数残差不能替代完整NS缺陷'));
   assert.ok(ids.diagnosticNote.textContent.includes('均不计NS'));
   for(const key of ['coarseConstructionModelApplications','coarsePreconditionerModelApplications','coarseReuseVerificationModelApplications'])assert.equal(diag[key],d.cases[j].summary[key]);
   if(diag.coarseBasis==='conservative-targeted-spatial'){
    assert.ok(diag.parentCoarseUnknowns<diag.coarseUnknowns&&diag.coarseTargetCellCount>0);
    assert.ok(ids.diagnosticNote.textContent.includes('不是物理网格加密'));
   }
   if(diag.coarseBasis==='mixed-conservative-primitive'){
    assert.equal(diag.coarseUnknowns,diag.coarseBaseConservativeColumns+diag.coarseAdditionalPrimitiveColumns);
    assert.ok(diag.coarseAdditionalPrimitiveColumns>0&&ids.diagnosticNote.textContent.includes('逐列提升身份固定'));
   }
   if(diag.coarseBasis==='local-conservative-enriched-mixed'){
    assert.equal(diag.coarseUnknowns,diag.coarseRetainedMixedColumns+diag.coarseLocalConservativeColumns);
    assert.equal(diag.coarseRetainedMixedColumns,diag.coarseBaseConservativeColumns+diag.coarseAdditionalPrimitiveColumns);
    assert.ok(diag.coarseLocalDistinctCells>0&&diag.coarseLocalConservativeColumns>=diag.coarseLocalDistinctCells);
    assert.ok(ids.diagnosticNote.textContent.includes('原CV与方程不删除')&&ids.diagnosticNote.textContent.includes('逐列提升身份固定'));
   }
  }
  if(diag.jacobianAction==='selected-analytic'){
   assert.ok(diag.selectedAnalyticApplications>0);
   assert.ok(ids.diagnosticNote.textContent.includes('组合方向已由原NS差分重新检查'));
   assert.ok(ids.diagnosticNote.textContent.includes('解析模型作用不计作NS调用'));
  }
  if(diag.wlsNormalCorrection){
   const correction=diag.wlsNormalCorrection;
   assert.ok(correction.stressHeatAndBothWorkTermsIncluded&&!correction.crossOriginalFacetWlsBlocksIncluded&&!correction.physicalEquationChanged);
   assert.ok(correction.additionalPersistentBlockUpperBoundBytes<=correction.maximumAdditionalBytes);
   assert.ok(ids.diagnosticNote.textContent.includes('WLS五带耦合')&&ids.diagnosticNote.textContent.includes('其他原面之间仍用既有图近似'));
  }
  if(diag.steadyInitialGuessModelApplications>0){
   assert.ok(diag.steadyInitialGuessReport.sourceStateAndGeometryVerified);
   assert.ok(ids.diagnosticNote.textContent.includes('分母仍是原始右端')&&ids.diagnosticNote.textContent.includes('不复用旧导数'));
   const guesses=linear.filter(a=>a.initialGuessUsed);assert.ok(guesses.length>0);
   guesses.forEach(a=>assert.ok(a.effectiveCorrectionTolerance>0&&a.effectiveCorrectionTolerance<1));
  }
  if(diag.krylovMethod==='SVD-selected-augmented-restarted-right-GMRES'){
   assert.ok(ids.diagnosticNote.textContent.includes('SVD多方向增强')&&ids.diagnosticNote.textContent.includes('不是进程峰值'));
   const recycled=linear.filter(a=>a.recycled&&a.recycled.method);
   assert.ok(recycled.length>0);assert.equal(recycled.reduce((n,a)=>n+a.recycled.arnoldiSteps,0),diag.recycledKrylovArnoldiSteps);
   recycled.forEach(a=>{assert.ok(a.recycled.basisStorageUpperBoundBytes<=a.recycled.maximumBasisBytes);assert.ok(a.recycled.modelApplications<=a.recycled.maximumModelApplications);});
  }
  linear.forEach((a,i)=>{
   const tr=ids.linearRows.children[i];assert.equal(tr.children.length,linearColumnCount);
   if(linearColumnCount===9){
    if(a.coarseAlgebraicRelativeResidual==null)assert.equal(tr.children[8].textContent,'未知');
    else assert.ok(Math.abs(Number(tr.children[8].textContent)-a.coarseAlgebraicRelativeResidual)<1e-8);
   }
   if(a.normalizedMeritDerivative===null){assert.equal(tr.children[6].textContent,'未知');assert.equal(tr.children[7].textContent,'未知');}
   else{assert.ok(Math.abs(Number(tr.children[6].textContent)-a.normalizedMeritDerivative)/Math.max(1,Math.abs(a.normalizedMeritDerivative))<1e-5);assert.ok(a.normalizedMeritDerivative<=a.inexactNewtonUpperBound+1e-12);}
  });
  linear.forEach((a,i)=>{assert.equal(ids.linearRows.children[i].children[3].textContent,a.fullSystemTargetMet?'通过':'拒绝');assert.equal(a.fullSystemTargetMet,a.fullSystemRelativeDefect<=a.target);});
  trials.forEach((a,i)=>{assert.equal(ids.trialRows.children[i].children[4].textContent,a.accepted?'接受':'拒绝');assert.equal(ids.trialRows.children[i].children[5].textContent,a.errorIdentifier||'无护栏错误');});
 }else{assert.equal(ids.linearRows.children.length,0);assert.equal(ids.trialRows.children.length,0);assert.equal(ids.layerRows.children.length,0);}
}
// 模块用途：用实际小域输出检查“温域通过/方向拒绝/预算不足”呈现，不将制造系统映射为真实流场。
let repreviewPresentationCases=0;
if(process.argv[4]){
 const fixture=JSON.parse(fs.readFileSync(process.argv[4],'utf8')),base=d.cases[2].solverDiagnostics.directionPreparation;
 const source=path.resolve(__dirname,'..',fixture.sourceFile);
 assert.equal(crypto.createHash('sha256').update(fs.readFileSync(source)).digest('hex'),fixture.sourceSHA256);
 ids.caseChoice.value='2';ids.caseChoice.events.change();
 for(const p of fixture.preparations){
  d.cases[2].solverDiagnostics.directionPreparation={...base,...p};buttons[3].events.click();
  const text=ids.details.children.map(e=>e.textContent).join('\n'),r=p.repreview;
  assert.ok(text.includes('预览 '+r.previewCount+' 次')&&text.includes('受限重算 '+p.activeSolveCount+' 次'));
  assert.ok(text.includes('不包含开放边界、完整NS、稳态或物理验证'));
  assert.equal(text.includes('指定步长全部单元与内部界面温域通过'),r.targetCandidateThermallyFeasible);
  assert.equal(text.includes('通过原128eps点积界'),p.bandDirectionFeasible);
  if(r.targetCandidateThermallyFeasible&&!p.bandDirectionFeasible)assert.ok(text.includes('内部温度合法也不能放行方向'));
  repreviewPresentationCases++;
 }
 d.cases[2].solverDiagnostics.directionPreparation=base;buttons[3].events.click();
}
const sourceSHA256=Object.fromEntries([__filename,program,path.join(folder,'index.html'),require.resolve('../web/viewer.js'),require.resolve('../web/export.js'),require.resolve('../web/section.js')].map(p=>[path.relative(path.resolve(__dirname,'..'),p),crypto.createHash('sha256').update(fs.readFileSync(p)).digest('hex')]));
if(d.finiteStepResearch){
 const f=d.finiteStepResearch;assert.equal(f.physicalSolution,false);assert.equal(f.newProductionSolvedAreaM2,0);
 buttons[3].events.click();assert.ok(ids.details.children.some(e=>e.textContent===f.patchComparisonLabel));
 const links=ids.details.children.flatMap(e=>e.children||[]).filter(e=>e.tag==='a');
 assert.ok(links.some(e=>e.attrs.href===f.reportURL));assert.ok(fs.existsSync(path.join(folder,f.reportURL)));
 assert.ok(f.fullComparison.branches.every(b=>b.history.every(r=>!r.steadyConverged)));
}
const report={status:'passed-three-real-wall-velocity-cases-CSV-sections-and-same-host-events',processID:process.pid,cases:3,sourceFaces:19752,
 CSVrows,sectionSamples,sameHostEvents:events,linearDiagnosticRows,candidateDiagnosticRows,repreviewPresentationCases,sourceSHA256,inPlaceDiagnosticNavigationVerified:Boolean(d.temperatureProjectionResearch||d.temperatureActiveResearch),browserVisualOrWebGLRuntimeVerified:false,newProductionSolvedAreaM2:0,dataSHA256:crypto.createHash('sha256').update(raw).digest('hex')};
fs.writeFileSync(path.join(folder,'node_wall_velocity_summary.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report));
