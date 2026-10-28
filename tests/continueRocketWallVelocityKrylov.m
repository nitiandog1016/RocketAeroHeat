function report=continueRocketWallVelocityKrylov(folder,previousFolder,options)
%CONTINUEROCKETWALLVELOCITYKRYLOV 接续同一静止壁约束方程，逐字段核对上一真实端点。
% 模块一：只合并同一新方程历史，原v82源方法历史不拼接。
% 模块二：仅在上一完整步被接受时按原增长规则设置下一伪步，不新增物理时间。
% 模块三：保存完整更新/拒绝及各段/累计调用；未稳态不增加区域合格面积。
if nargin<3,options=struct;end
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'WallVelocityContinue:Overwrite','不得覆盖真实接续。');mkdir(folder);
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');
s=load(fullfile(previousFolder,'wall_velocity_krylov_evidence.mat'));
assert(strcmp(cutCellWallVelocityPolicy(s.boundary.solidWall),'dirichlet-tangent'));
options.residualScale=s.result.residualScale;
if ~isfield(options,'initialPseudoDt')
    tau=s.result.history.pseudoDt(end);
    if tau==0,tau=s.result.options.initialPseudoDt;
    elseif s.result.history.alpha(end)==1,tau=min(s.result.options.maximumPseudoDt,tau*s.result.options.pseudoDtGrowth);end
    options.initialPseudoDt=tau;
end
if ~isfield(options,'maximumIterations'),options.maximumIterations=3;end
if ~isfield(options,'maximumKrylovIterations'),options.maximumKrylovIterations=6;end
% 仅保留接续核对所需数据；旧图矩阵/失败方向/旧候选不必同时驻留求解内存。
previousQ=s.result.Q;previousRate=s.result.rate;previousInfo=s.result.info;
if isfield(s,'displayHistory'),previousHistory=s.displayHistory;else,previousHistory=s.result.history;end
priorEvaluations=s.result.residualEvaluations;
if isfield(s.report,'totalNewEquationResidualEvaluations'),priorEvaluations=s.report.totalNewEquationResidualEvaluations;end
priorReport=s.report;lineage=s.lineage;gas=s.gas;boundary=s.boundary;clear s;
timer=tic;result=solveRocketShellSteadyKrylov(previousQ,d.domain,gas,boundary,options);elapsed=toc(timer);
assert(isequaln(result.initialQ,previousQ),'WallVelocityContinue:Identity','接续初态不同。');
replayReport=assertRocketNSReplay(result.initialRate,result.initialInfo,previousRate,previousInfo,d.domain.cells.volume);
[displayHistory,work]=mergeRocketKrylovHistory(previousHistory,priorEvaluations,result.history,result.residualEvaluations);
priorAccepted=work.priorAcceptedIterations;
mapped=mapCutCellWallFields(d.domain,result.info,result.steadyConverged);assert(~mapped.eligibleForSurfaceCoverage);
report=struct('status','bounded-real-20deg-same-static-wall-equation-Krylov-continuation', ...
    'sourceFaces',size(d.domain.inputGeometry.F,1),'cellCount',size(result.Q,1),'acceptedNonlinearIterations',result.acceptedNonlinearIterations, ...
    'priorAcceptedNonlinearIterations',priorAccepted,'totalNewEquationAcceptedIterations',priorAccepted+result.acceptedNonlinearIterations, ...
    'candidateTrials',result.candidateTrials,'residualEvaluations',result.residualEvaluations, ...
    'totalNewEquationResidualEvaluations',priorEvaluations+result.residualEvaluations,'operatorBuilds',result.operatorBuilds, ...
    'initialMaximumResidual',result.initialMaximumResidual,'maximumResidual',result.maximumResidual, ...
    'initialRmsMerit',result.initialRmsMerit,'rmsMerit',result.rmsMerit,'elapsedSeconds',elapsed, ...
    'matrixVectorProducts',result.matrixVectorProducts,'jacobianResidualEvaluations',result.jacobianResidualEvaluations, ...
    'preconditionerApplications',result.preconditionerApplications,'lineage',lineage,'wallVelocityGradient','dirichlet-tangent', ...
    'sameEquationHistoryContinued',true,'steadyConverged',result.steadyConverged,'stopReason',result.stopReason, ...
    'spatialConservationRelative',result.spatialConservationRelative,'steadyBoundaryFluxBalance',result.steadyBoundaryFluxBalance, ...
    'signedHeatPowerW',mapped.summedKnownHeatPower,'wallHeatRangeWm2',[min(mapped.heatFluxIntoWall.value),max(mapped.heatFluxIntoWall.value)], ...
    'physicalTimeAdvanced',0,'BETimeStepsClaimed',0,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
report.priorActualResidualEvaluations=priorEvaluations;report.replayReport=replayReport;
% 新储库问题接续时保留其边界身份，只累计同一边界问题；旧特征历史仍不拼接。
if isfield(priorReport,'openBoundaryDiscretization')
    report.openBoundaryDiscretization=priorReport.openBoundaryDiscretization;
    report.sourceHistoryAppended=false;
end
report.residualEvaluationsAfterLastAcceptedHistory=work.residualEvaluationsAfterLastAcceptedHistory;
if isfield(result,'normalEulerLinearization')
    report.normalEulerLinearization=result.normalEulerLinearization;
    report.physicalSpatialEquationChanged=false;
end
if isfield(result,'linearEquation')
    report.linearEquation=result.linearEquation;report.preconditionerTimeStepOnly=result.preconditionerTimeStepOnly;
end
if isfield(result,'rightPreconditioner')
    report.rightPreconditioner=result.rightPreconditioner;report.graphSweeps=result.options.graphSweeps;
    report.fixedLinearPreconditioner=result.fixedLinearPreconditioner;
end
if isfield(result,'jacobianAction')
    report.jacobianAction=result.jacobianAction;report.selectedAnalyticApplications=result.selectedAnalyticApplications;
    report.analyticCacheBuildSeconds=result.analyticCacheBuildSeconds;
    report.originalNSCombinedDirectionCheckPreserved=result.originalNSCombinedDirectionCheckPreserved;
end
if isfield(result,'steadyInitialGuessReport')
    report.steadyInitialGuessReport=result.steadyInitialGuessReport;
    report.steadyInitialGuessModelApplications=result.steadyInitialGuessModelApplications;
end
if isfield(result,'krylovMethod')
    for key={'krylovMethod','recycledKrylovModelApplications','recycledKrylovPreconditionerApplications','recycledKrylovArnoldiSteps'},report.(key{1})=result.(key{1});end
end
if isfield(result,'wlsNormalCorrection'),report.wlsNormalCorrection=result.wlsNormalCorrection;end
if isfield(result,'coarseUnknowns')
    for key={'coarseUnknowns','coarseImageBytes','coarseOperatorBuilds','coarseOperatorReuses', ...
            'coarseConstructionModelApplications','coarsePreconditionerModelApplications', ...
            'coarseReuseVerificationModelApplications','coarseReuseReport','coarseSetupSeconds'}
        report.(key{1})=result.(key{1});
    end
    report.coarseCorrectionChangesPhysicalMesh=false;
    if isfield(result,'coarseBasis')
        for key={'coarseBasis','coarseTargetCellCount','parentCoarseUnknowns','coarseBaseConservativeColumns','coarseAdditionalPrimitiveColumns', ...
                'coarseRetainedMixedColumns','coarseLocalConservativeColumns','coarseLocalDistinctCells'}
            if isfield(result,key{1}),report.(key{1})=result.(key{1});end
        end
    end
end
save(fullfile(folder,'wall_velocity_krylov_evidence.mat'),'result','report','mapped','gas','boundary','lineage','displayHistory','priorReport','-v7');
fid=fopen(fullfile(folder,'wall_velocity_krylov_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
