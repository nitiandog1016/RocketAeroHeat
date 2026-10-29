/* 模块用途：真实原面静止墙梯度三组对照，保留给定量、算法和独立迭代历史。 */
'use strict';
(function(){
 const d=window.ROCKET_WALL_VELOCITY_STUDY,$=id=>document.getElementById(id),buttons=[...document.querySelectorAll('[data-step]')];
 const num=x=>Number.isFinite(x)?Number(x.toPrecision(6)).toString():'未知';let selected=d.cases.length-1,step=4;
 d.cases.forEach((c,i)=>{const o=document.createElement('option');o.value=String(i);o.textContent=c.label;$('caseChoice').append(o);});$('caseChoice').value=String(selected);
 function panel(){
  const c=d.cases[selected],s=c.summary,scheme=c.wallVelocityGradient==='legacy'?'旧WLS加斜位移约束':'静止壁零切向导数与真实正dn约束';
  const stopReasons={'maximum-nonlinear-iterations':'达到本段修正次数上限','bounded-line-search-and-pseudostep-retries-exhausted':'有限线搜索与伪步重试已用尽（未收敛）','steady-newton-linear-target-not-met':'直接稳态线性目标未达到（未提交更新）','original-spatial-residual-and-steady-boundary-balance-met':'已满足原稳态数值门槛'};
  const stopReason=stopReasons[s.stopReason]||(s.stopReason||'基线或冻结状态，本页未迭代');
  const panels=[
   ['输入模型 · 真实原箭体','保持原19752面及20°完整三维来流。',[['原三角面',d.sourceFaces],['坐标','SI / X鼻尖正向 / 零攻角来流-X']]],
   ['生成网格 · 原三维层域','真实面法向层域；法向线只作右预条件，完整NS单独验收。',[['控制体 / 层数',d.cellCount+' / '+d.layerCount],['网格资格','几何已验；近壁/网格/外域物理敏感性待验证']]],
   ['添加边界条件 · 静止等温墙','相同三维特征储库，Uw=[0,0,0]与Tw300K始终为输入。',[['墙速度梯度',scheme],['给定切向导数','常数零墙速的精确导数；没有猜未知导数'],['气体','常γ='+d.gas.gamma+' / 冻结μ / Pr='+d.gas.Pr]]],
   ['计算 · '+c.label,'改变墙应力离散后从新基线计数，不拼接旧方程历史。',[['墙速度梯度',scheme],['本段修正 / NS调用',s.acceptedNonlinearIterations+' / '+s.residualEvaluations],['累计新方程修正 / NS调用',(s.totalNewEquationAcceptedIterations??s.acceptedNonlinearIterations)+' / '+(s.totalNewEquationResidualEvaluations??s.residualEvaluations)],['停止原因',stopReason],['最大残差',num(s.maximumResidual)+' s⁻¹（目标≤2e−6）'],['RMS残差',num(s.rmsMerit)+' s⁻¹'],['空间通量账误差',num(s.spatialConservationRelative)+'（目标≤1e−9）'],['定常边界不平衡',num(s.steadyBoundaryFluxBalance)+'（目标≤5e−10）'],['稳态','空间账通过不代表残差或定常边界收敛']]],
   ['查看结果 · '+c.label,'压力、热流、温度和代理耗散诊断共用当前状态、原面及色标。',[['邻壁气体温度',c.adjacentTemperatureRangeK.map(num).join(' ～ ')+' K'],['有符号入壁热功率',num(c.signedHeatPowerW/1000)+' kW（研究值）'],['负机械耗散墙面',c.mechanicalCoolingWallCount],['静止壁外功','0 W；与代理黏性内能转换分开'],['新增合格面积','0 m²']]]
  ];
  if(s.openBoundaryDiscretization==='exact-Riemann-reservoir'){
   panels[2][1]='显式三维精确Riemann储库，在xi=0采样真实波系；黏性按实际质量方向闭合。';
   panels[2][2].push(['开放边界方法','精确Riemann储库（数值研究）']);
   panels[3][1]='本边界问题从0单独计数，旧特征问题的13次接受和258次NS只作来源。';
   panels[3][2].push(['历史身份','新储库问题独立历史；未拼接旧特征问题']);
  }
  // 模块用途：同一段初末态的整体和局部趋势并列，避免小步接受掩盖局部恶化。
  const trend=c.solverDiagnostics&&c.solverDiagnostics.segmentResidualTrend;
  const coordinates=c.solverDiagnostics&&c.solverDiagnostics.trialStateCoordinates;
  if(coordinates)panels[3][2].push(['非线性候选路径',coordinates==='primitive-tangent'?'原始变量切向更新后恢复守恒量；每步仍由完整NS验收':'守恒变量直线更新；每步由完整NS验收']);
  const preparation=c.solverDiagnostics&&c.solverDiagnostics.directionPreparation;
  if(preparation){
   panels[3][2].push(['本次方向构造',preparation.methodLabel],
    ['近似父矩阵C剩余',num(preparation.internalCRelativeResidual)+(preparation.internalCTargetMet?'；内部目标通过':'；未达内部目标 '+num(preparation.internalCTarget))],
    ['完整模型 / 原NS缺陷',num(preparation.completeModelRelativeDefect)+' / '+num(preparation.originalNSDirectionDefect)+'（原NS目标≤'+num(preparation.originalNSDirectionTarget)+'）'],
    ['幅度控制身份','峰值软惩罚；不保证硬幅度上限，实际候选仍由完整NS验收']);
   if(preparation.displayedStateUnchanged)panels[3][2].push(['本次方向状态','未接受流场更新；云图保持原检查点，拒绝诊断单独列出']);
   if(preparation.kind==='temperature-comparison-QR-LSQR-v1')panels[3][2].push(['界面温度方向约束',preparation.constraintApplied?'原内部面19207：一阶温度变化约为0；未增加物理温度边界，非线性候选仍须原温域检查':'同态自由方向对照；未施加界面温度投影']);
   if(preparation.kind==='current-band-active-temperature-v1'){
    const count=v=>Array.isArray(v)?v.length:v==null?0:1;
    const lookahead=preparation.lookaheadSelection;
    const selectionLabel=lookahead?'原10 K带内 '+count(lookahead.originalBandFaceIDs)+' 面；'+(preparation.repreview?'累计暴露越温 ':'预测越温 ')+count(lookahead.predictedViolatingCellIDs)+' 单元；实际选择 '+count(preparation.bandFaceIDs)+' 面，激活 '+count(preparation.activeFaceIDs)+' 行':'当前两端10 K带内 '+count(preparation.bandFaceIDs)+' 面；实际激活 '+count(preparation.activeFaceIDs)+' 行';
    panels[3][2].push(['当前温域方向选择',selectionLabel],
     ['方向与边界身份','每步重建当前导数；有限步必要条件可强制激活。未增加物理温度边界，实际候选仍须原NS和温域检查']);
    if(preparation.repreview){
     // 模块用途：投影后末轮与初始自由预览分开；内部温域通过不代替完整NS、开放边界或稳态。
     const r=preparation.repreview,h=Array.isArray(r.history)?r.history:[r.history],last=h.at(-1);
     const reasons={'specified-target-candidate-temperature-checked-with-all-current-rows-feasible':'指定步长内部温域通过；原方向门槛另验','repreview-has-no-new-actionable-original-condition':'没有可新增条件，停止受限重算','bounded-repreview-budget-exhausted':'已用尽本次受限重算预算'};
     panels[3][2].push(['投影后温域复查','目标步长 '+num(r.targetAlpha)+'；预览 '+r.previewCount+' 次；受限重算 '+preparation.activeSolveCount+' 次'],
      ['末轮内部温域',!last||!last.previewCellStateEvaluable?'未执行完整内部检查：候选单元状态无法评估':
       r.targetCandidateThermallyFeasible?'指定步长全部单元与内部界面温域通过；完整NS另验':
       '未通过：越温单元 '+count(last.predictedViolatingCellIDs)+'；可检查越温界面 '+count(last.predictedViolatingInterfaceFaceIDs)+'；阻挡面 '+last.blockedInterfaceCount],
      ['原方向精度门槛',typeof preparation.bandDirectionFeasible!=='boolean'?'未知：缺少原方向检查结果':preparation.bandDirectionFeasible?'通过原128eps点积界及全部必要条件检查':'未通过；内部温度合法也不能放行方向'],
      ['复查停止原因',reasons[preparation.repreviewStopReason]||preparation.repreviewStopReason||'未知'],
      ['复查适用范围','仅指定步长的全部单元和内部Riemann温域；不包含开放边界、完整NS、稳态或物理验证']);
    }else if(lookahead&&lookahead.interfacePreview&&!Array.isArray(lookahead.interfacePreview)){
     const preview=lookahead.interfacePreview;
     panels[3][2].push(['Riemann界面前瞻','预测越温 '+count(preview.predictedViolatingInterfaceFaceIDs)+' 面；已检查 '+count(preview.assessedFaceIDs)+' 面；因单元超域阻挡 '+count(preview.blockedFaceIDs)+' 面'],
      ['前瞻适用范围',preview.allInteriorFacesAssessed?'本次自由候选的全部内部面已检查；投影后候选仍须完整NS验收':'仅检查候选输入单元合法的内部面；阻挡面尚未检查，不能宣称候选合法']);
    }else if(lookahead&&['original-and-predicted-cell-and-interface-faces','original-and-iteratively-predicted-cell-and-interface-faces'].includes(lookahead.selectionPolicy)){
     panels[3][2].push(['Riemann界面前瞻','未执行：目标候选的单元状态无法评估；没有界面已检查结论']);
    }
    if(preparation.freeCheckpointReplay)panels[3][2].push(['检查点接续','同一当前状态和方程已核对；复用自由方向 '+preparation.reusedFreeNativeIterations+' 步，新算受限方向 '+preparation.newProjectedNativeIterations+' 步']);
   }
  }
  // 当前接受态与方向比较分别绑定，避免把并列分支拼成连续迭代。
  const temperatureStudy=d.temperatureProjectionResearch||d.temperatureActiveResearch;
  let trendNote='';
  if(trend){
   panels[3][2].push(['本段RMS变化',num(100*trend.rmsRelativeChange)+' %'],['本段最大残差变化',num(100*trend.maximumRelativeChange)+' %']);
   if(trend.rmsDecreased&&trend.maximumIncreased)trendNote='本段RMS下降，但最大局部残差升高；一次接受不能认定稳定收敛。';
   else if(!trend.steadyNumericalTargetsMet)trendNote='当前未满足稳态数值门槛；接受更新不能代表区域验证通过。';
   if(trendNote){panels[3][2].push(['收敛判断',trendNote]);panels[4][2].push(['收敛判断',trendNote]);}
  }
  // 有限步诊断来自独立冻结记录，不将内部小域状态映射成整箭物理解。
  const finite=d.finiteStepResearch;
  if(finite){
   panels[3][2].push(['有限步诊断起点',finite.sourceLabel],['主要机制',finite.mechanismLabel],
    ['同预算六步小复现',finite.patchComparisonLabel],['整箭修正对照',finite.fullComparisonLabel],
    ['当前资格','未收敛研究；新增正式求解面积0']);
   panels[4][2].push(['有限步诊断',finite.scopeLabel]);
  }
  const p=panels[step];$('title').textContent=p[0];$('description').textContent=p[1];$('details').replaceChildren();
  p[2].forEach(([key,value])=>{const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=key;dd.textContent=String(value);$('details').append(dt,dd);});
  if(finite&&(step===3||step===4)){
   const dt=document.createElement('dt'),dd=document.createElement('dd'),a=document.createElement('a');
   dt.textContent='诊断与连续步记录';a.textContent='查看对照图、各步残差及适用限制';a.setAttribute('href',finite.reportURL);
   a.addEventListener('click',event=>{event.preventDefault();$('viewer').setAttribute('src',finite.reportURL);});dd.append(a);$('details').append(dt,dd);
  }
  if(temperatureStudy&&temperatureStudy.selectedCase===s.resultEvidenceFolder&&(step===3||step===4)){
   const dt=document.createElement('dt'),dd=document.createElement('dd'),detail=document.createElement('details'),summary=document.createElement('summary'),a=document.createElement('a');
   dt.textContent='方向对照与温域限制';summary.textContent='展开定量比较与适用限制';detail.append(summary);
   [temperatureStudy.comparisonLabel,temperatureStudy.scopeLabel,temperatureStudy.constraintLabel].forEach(text=>{const p=document.createElement('p');p.textContent=text;detail.append(p);});
   a.textContent='在当前工作台查看定量对照';a.setAttribute('href',temperatureStudy.reportURL);
   a.addEventListener('click',event=>{event.preventDefault();$('viewer').setAttribute('src',temperatureStudy.reportURL);});detail.append(a);dd.append(detail);$('details').append(dt,dd);
  }
  buttons.forEach((b,i)=>b.setAttribute('aria-selected',String(i===step)));
  const url=step===4?c.resultsURL:step===1?d.domainURL:d.geometryURL;if($('viewer').getAttribute('src')!==url)$('viewer').setAttribute('src',url);
  $('historyRows').replaceChildren();const h=c.history;
  if(h&&Array.isArray(h.iteration))h.iteration.forEach((it,j)=>{const tr=document.createElement('tr');[it,num(h.pseudoDt[j]),num(h.alpha[j]),num(h.maximumResidualPerSecond[j]),num(h.rmsMerit[j]),num(h.spatialConservationRelative[j]),num(h.steadyBoundaryFluxBalance[j]),h.residualEvaluations[j]].forEach(x=>{const td=document.createElement('td');td.textContent=String(x);tr.append(td);});$('historyRows').append(tr);});
  // 按需展示本段真实失败检查；不将线预条件器的微小缺陷当成完整NS线性精度。
  const diag=c.solverDiagnostics,rows=x=>Array.isArray(x)?x:x&&Object.keys(x).length?[x]:[];
  $('linearRows').replaceChildren();$('trialRows').replaceChildren();$('layerRows').replaceChildren();
  const layerFractions=diag&&diag.linearRemainingLayerFractions,equationFractions=diag&&diag.linearRemainingEquationFractions;
  if(Array.isArray(layerFractions)&&layerFractions.length){
   layerFractions.forEach((fraction,i)=>{const tr=document.createElement('tr');[i+1,num(100*fraction)].forEach(value=>{const td=document.createElement('td');td.textContent=String(value);tr.append(td);});$('layerRows').append(tr);});
   const names=['质量','X动量','Y动量','Z动量','总能量'];
   $('layerNote').textContent='最后一次原NS组合方向复核的线性剩余，按原缩放五量平方范数分解；不是压力、热流误差或物理解。方程占比：'+names.map((name,i)=>name+' '+num(100*equationFractions[i])+'%').join('，')+'。';
  }else $('layerNote').textContent='当前状态没有可用的非零原NS方向剩余分布。';
  if(diag&&diag.localSegmentOnly){
   const lineMethod=diag.normalEulerLinearization==='exact-branch'?'精确当前波系分支导数':'冻结Rusanov导数';
   const equation=diag.linearEquation==='steady-newton'?'实际线性方程为(-J)δ=R，τ仅用于右预条件器；下降斜率只表示数值模型。':'实际线性方程为(I/τ-J)δ=R。';
   const graphMode=['lateral-colored-line','two-level-minres'].includes(diag.rightPreconditioner);
   let strategy=graphMode?'侧向邻居Euler也采用精确当前波系分支导数；固定 '+num(diag.graphSweeps)+' 次对称着色法向线扫描，图残差仅为预条件诊断，不是精确图解。':'独立法向线精确求解近似预条件系统。';
   if(diag.wlsNormalCorrection)strategy+='本段显式保留同一原面法向线内的WLS五带耦合，包含应力、导热及两项黏性功；其他原面之间仍用既有图近似。新增持久块存储保守上界 '+num(diag.wlsNormalCorrection.additionalPersistentBlockUpperBoundBytes/1048576)+' MiB（不是进程峰值）。这是预条件修正，原完整NS和验收门槛保持。';
   if(diag.rightPreconditioner==='two-level-minres')strategy+='随后用 '+num(diag.coarseUnknowns)+' 个固定粗方向最小化细层模型残差，粗像存储 '+num(diag.coarseImageBytes/1048576)+' MiB；原物理网格不变。粗层代数残差不能替代完整NS缺陷。本段构建/复用 '+num(diag.coarseOperatorBuilds)+' / '+num(diag.coarseOperatorReuses)+'；构建、预条件、复用检查分别增加 '+num(diag.coarseConstructionModelApplications)+' / '+num(diag.coarsePreconditionerModelApplications)+' / '+num(diag.coarseReuseVerificationModelApplications)+' 次模型作用，均不计NS。';
   if(diag.coarseBasis==='conservative-targeted-spatial')strategy+='本段采用鼻端近壁局部空间细分，父空间 '+num(diag.parentCoarseUnknowns)+' 个未知量，目标 '+num(diag.coarseTargetCellCount)+' 个原控制体；新子组保持连通和真实体积，不是物理网格加密。';
   if(diag.coarseBasis==='mixed-conservative-primitive')strategy+='混合粗空间保留 '+num(diag.coarseBaseConservativeColumns)+' 个保守底空间列，附加 '+num(diag.coarseAdditionalPrimitiveColumns)+' 个原始变量耦合方向；逐列提升身份固定，右端变化不重新筛选，不是物理网格加密。';
   if(diag.coarseBasis==='local-conservative-enriched-mixed')strategy+='局部增强保留 '+num(diag.coarseRetainedMixedColumns)+' 个混合列，增加 '+num(diag.coarseLocalConservativeColumns)+' 个保守方向，作用于 '+num(diag.coarseLocalDistinctCells)+' 个原控制体；全部 '+num(diag.coarseBaseConservativeColumns)+' 个保守底列保留。逐列提升身份固定，原CV与方程不删除，右端变化不重新筛选，不是物理网格加密。';
   if(diag.steadyInitialGuessModelApplications>0)strategy+='同一基态保存方向作为初猜，重新执行 '+num(diag.steadyInitialGuessModelApplications)+' 次当前模型作用形成剩余；不复用旧导数。原NS复核的是初猜与增量相加的总方向，分母仍是原始右端。GMRES报告与容差以增量右端为分母，不能与实际NS缺陷混同。'+rows(diag.linearAttempts).filter(a=>a.initialGuessUsed).map(a=>'初猜模型缺陷 '+num(a.initialGuessModelRelativeDefect)+'，增量容差 '+num(a.effectiveCorrectionTolerance)+(a.gmresSkipped?'；模型阈值已满足，跳过GMRES但仍检查原NS。':'。')).join('');
   if(diag.steadyInitialGuessReport&&diag.steadyInitialGuessReport.sourceVariant==='ns-probe')strategy+='本段初猜检查实际调用当前完整NS，已计入NS工作量；完整几何独立绑定，不加载旧粗层像。';
   if(diag.krylovMethod==='SVD-selected-augmented-restarted-right-GMRES')strategy+='本段显式采用SVD多方向增强重启，实际Arnoldi步数 '+num(diag.recycledKrylovArnoldiSteps)+'，其完整模型/预条件作用 '+num(diag.recycledKrylovModelApplications)+' / '+num(diag.recycledKrylovPreconditionerApplications)+'；构造保留像与周期复核均计入模型预算。'+rows(diag.linearAttempts).filter(a=>a.recycled&&a.recycled.method).map(a=>'每周期最多 '+num(a.recycled.options.restartLength)+' 步，最多 '+num(a.recycled.options.maximumCycles)+' 周期，保留最多 '+num(a.recycled.options.retainedDirections)+' 个方向；基数据保守上界 '+num(a.recycled.basisStorageUpperBoundBytes/1048576)+' MiB（不是进程峰值）。实际各周期保留 '+rows(a.recycled.cycles).map(c=>c.retainedBefore+'→'+c.retainedAfter).join(' / ')+'。').join('')+'表中预条件残差取最后一次Arnoldi预条件应用；完整NS总方向缺陷仍独立验收。';
   const action=diag.jacobianAction==='selected-analytic'?'本段使用选中分支的完整WLS/特征边界解析方向模型（共 '+num(diag.selectedAnalyticApplications)+' 次，含可能的初猜检查）；组合方向已由原NS差分重新检查。解析模型作用不计作NS调用。':'GMRES逐方向调用原NS差分。';
   $('diagnosticNote').textContent='仅本段 '+diag.sourceFolder+'，承接此前 '+diag.priorAcceptedCorrections+' 次修正；累计接受历史见上表。'+equation+action+'法向线Euler使用'+lineMethod+'，仅为预条件，正式NS不改。'+strategy+'完整NS线性缺陷须≤'+num(diag.linearRelativeTolerance)+'，预条件器残差另列。接受历史记录当时调用数；累计最后记录状态后另有 '+(diag.residualEvaluationsAfterLastAcceptedHistory??0)+' 次NS调用，均计入累计总数；本段调用数单列。';
   rows(diag.linearAttempts).forEach(a=>{const tr=document.createElement('tr');[a.attempt,num(a.fullSystemRelativeDefect),num(a.target),a.fullSystemTargetMet?'通过':'拒绝',num(a.gmresReportedResidual),num(a.preconditionerRelativeDefect),num(a.normalizedMeritDerivative),num(a.inexactNewtonUpperBound),num(a.coarseAlgebraicRelativeResidual)].forEach(x=>{const td=document.createElement('td');td.textContent=String(x);tr.append(td);});$('linearRows').append(tr);});
   rows(diag.nonlinearTrials).forEach(a=>{const tr=document.createElement('tr');[a.iteration,a.retry,num(a.pseudoDt),num(a.alpha),a.accepted?'接受':'拒绝',a.errorIdentifier||'无护栏错误',num(a.candidateRmsMerit)].forEach(x=>{const td=document.createElement('td');td.textContent=String(x);tr.append(td);});$('trialRows').append(tr);});
  }else $('diagnosticNote').textContent='当前为冻结源状态或同Q基线，没有在此状态重新执行非线性修正。';
  // 后续内部求逆研究与所选物理场分开；拒绝方向不能改写云图或接受历史。
  $('linearResearchRows').replaceChildren();
  if(d.linearResearch&&diag&&diag.sourceFolder===d.linearResearch.sourceCase){
   $('linearResearchNote').textContent='以下为同一基态的后续内部线性研究，未更新当前云图或接受历史。C剩余属于近似矩阵，不能替代原NS缺陷、热流误差或面积覆盖。GMRES的正交与组合关系门槛均为1e−8；LSQR单列伴随检查，Arnoldi指标标为不适用。未执行指标显示未知，分解拒绝不能算作已解结果。';
   rows(d.linearResearch.studies).forEach(a=>{const tr=document.createElement('tr');
    const lsqr=[d.linearResearch.fullFLeastSquaresDiagnostics,d.linearResearch.fullFImageScaleDiagnostics,d.linearResearch.parentLeastSquaresDiagnostics].some(q=>q&&a.caseID===q.caseID);
    const iterationMetric=v=>lsqr&&v==null?'不适用（LSQR）':num(v);
    [a.label,num(a.innerRelativeCResidual),iterationMetric(a.basisOrthogonalityError),iterationMetric(a.modelRelationRelativeDifference),a.displayState||(a.candidateUsable?'内部方向通过':'拒绝'),a.newActualNSEvaluations,a.note||'内部方向检查'].forEach(x=>{const td=document.createElement('td');td.textContent=String(x);tr.append(td);});$('linearResearchRows').append(tr);
   });
  }else $('linearResearchNote').textContent='当前所选状态没有绑定后续内部线性研究。';
  $('wallNote').textContent='固定初猜的同Q受控检查：负机械耗散墙面 '+d.integrationReport.oldMechanicalCoolingWalls+' → '+d.integrationReport.newMechanicalCoolingWalls+'，负静止墙熵配对 '+d.integrationReport.oldNegativeWallEntropyFaces+' → '+d.integrationReport.newNegativeWallEntropyFaces+'；墙压/Fourier、内部通量、质量/能量残差保持。不是全NS熵稳定或物理资格证明。';
  $('stateNote').textContent='当前 '+c.label+'；温度变化率不是壁温历程，黏性内能转换来自邻壁速度代理，不是静止壁外功。所有研究有限值不计正式合格面积。'+trendNote;
 }
 buttons.forEach((b,i)=>b.addEventListener('click',()=>{step=i;panel();}));$('caseChoice').addEventListener('change',()=>{selected=Number($('caseChoice').value);panel();});panel();
})();
