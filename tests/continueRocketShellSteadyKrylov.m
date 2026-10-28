function report=continueRocketShellSteadyKrylov(folder,previousFolder,options)
%CONTINUEROCKETSHELLSTEADYKRYLOV 从真实保存端点接续有界完整三维Krylov，不伪造时间/历史。
% 模块一：重新调用原完整初始算子，并逐字段核对前一真实末态；不沿用仅摘要。
% 模块二：本次与此前非线性计数分开，合并历史只用于展示，原证据不改。
% 模块三：所有空间/线性/状态门槛保持，仍不取得物理面积资格。
if nargin<3,options=struct;end
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'ShellKrylovContinue:Overwrite','不得覆盖真实接续记录。');mkdir(folder);
domainPath=fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat');
if ~exist(domainPath,'file'),domainPath=fullfile(root,'validation','v80','domain02','shell_domain_evidence.mat');end
d=load(domainPath,'domain');s=load(fullfile(previousFolder,'steady_krylov_evidence.mat'),'result','report','gas','boundary');
options.residualScale=s.result.residualScale;
if ~isfield(options,'initialPseudoDt'),options.initialPseudoDt=s.result.history.pseudoDt(end);end
if ~isfield(options,'maximumIterations'),options.maximumIterations=3;end
t=tic;result=solveRocketShellSteadyKrylov(s.result.Q,d.domain,s.gas,s.boundary,options);elapsed=toc(t);
assert(isequaln(result.initialQ,s.result.Q)&&isequaln(result.initialRate,s.result.rate)&&isequaln(result.initialInfo,s.result.info), ...
    'ShellKrylovContinue:Source','接续初态完整算子没有复现上一末态。');
mapped=mapCutCellWallFields(d.domain,result.info,result.steadyConverged);priorAccepted=s.result.acceptedNonlinearIterations;
if isfield(s.report,'totalAcceptedNonlinearIterations'),priorAccepted=s.report.totalAcceptedNonlinearIterations;end
report=struct('status','bounded-real-20deg-full-Krylov-continuation-from-exact-previous-state', ...
    'sourceFaces',19752,'cellCount',numel(d.domain.cells.volume),'acceptedNonlinearIterations',result.acceptedNonlinearIterations, ...
    'priorAcceptedNonlinearIterations',priorAccepted,'totalAcceptedNonlinearIterations',priorAccepted+result.acceptedNonlinearIterations, ...
    'candidateTrials',result.candidateTrials,'residualEvaluations',result.residualEvaluations,'operatorBuilds',result.operatorBuilds, ...
    'initialMaximumResidual',result.initialMaximumResidual,'maximumResidual',result.maximumResidual, ...
    'initialRmsMerit',result.initialRmsMerit,'rmsMerit',result.rmsMerit,'elapsedSeconds',elapsed, ...
    'matrixVectorProducts',result.matrixVectorProducts,'jacobianResidualEvaluations',result.jacobianResidualEvaluations, ...
    'preconditionerApplications',result.preconditionerApplications,'viscousWorkVelocity',cutCellViscousWorkPolicy(s.boundary), ...
    'steadyConverged',result.steadyConverged,'stopReason',result.stopReason,'spatialConservationRelative',result.spatialConservationRelative, ...
    'steadyBoundaryFluxBalance',result.steadyBoundaryFluxBalance,'physicalTimeAdvanced',0,'BETimeStepsClaimed',0, ...
    'newProductionSolvedAreaM2',0,'physicalSolution',false,'signedHeatPowerW',mapped.summedKnownHeatPower, ...
    'wallHeatRangeWm2',[min(mapped.heatFluxIntoWall.value),max(mapped.heatFluxIntoWall.value)],'engineProcessID',feature('getpid'));
displayHistory=s.result.history;start=result.history;
for name=fieldnames(start)'
    key=name{1};values=start.(key)(2:end);
    if strcmp(key,'iteration'),values=values+priorAccepted;end
    if strcmp(key,'residualEvaluations'),values=values+displayHistory.residualEvaluations(end);end
    displayHistory.(key)=[displayHistory.(key),values];
end
gas=s.gas;boundary=s.boundary;priorReport=s.report;
save(fullfile(folder,'steady_krylov_evidence.mat'),'result','mapped','report','displayHistory','priorReport','gas','boundary','-v7');
fid=fopen(fullfile(folder,'steady_krylov_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
