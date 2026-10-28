function report=testRocketShellCoarsePreconditioner(folder)
%TESTROCKETSHELLCOARSEPRECONDITIONER 实际NS两层接入、缓存拒绝与旧解析路径完整回放。
% 模块一：复用已核对小域粗像，组合方向和实际候选仍交原NS验收。
% 模块二：受限一方向预算不能绕过实际线性目标；原缺省及解析图路径逐字段保持。
% 模块三：错误基态、边界、源码、像列、因子、体积空间及预算必须拒绝复用。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'CoarseSolverTest:Overwrite','不得覆盖两层接入验证。');mkdir(folder);
s=load(fullfile(root,'outputs','v83_validation','selectedsolver01','steady_newton_unit_evidence.mat'),'result','cut','gas','boundary');
cut=s.cut;gas=s.gas;boundary=s.boundary;old=s.result;Q=old.initialQ;options=old.options;
baseline=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
assert(isequaln(untimed(baseline),untimed(old)),'CoarseSolverTest:Compatibility','旧解析图完整输出未复现。');
sourceFile=fullfile(root,'outputs','v83_validation','coarselsunit01','least_squares_coarse_evidence.mat');
d=load(sourceFile,'coarse');space=d.coarse.space;cache=buildCutCellSelectedDirectionCache(Q,old.initialInfo,cut,gas,boundary,old.residualScale);
[loaded,reuse]=loadRocketLeastSquaresCoarseOperator(sourceFile,cache,space,256*1024^2);
assert(isequaln(loaded,d.coarse)&&reuse.reuseVerificationModelApplications==1);
rejected=0;
for kind=1:8
    coarse=d.coarse;limit=256*1024^2;
    switch kind
        case 1,coarse.baseQ(1,1)=1.01*coarse.baseQ(1,1);
        case 2,coarse.boundary.solidWall.temperature=301;
        case 3,coarse.sourceIdentity.applyCutCellSelectedDirection='wrong';
        case 4,coarse.image(:,1)=1.01*coarse.image(:,1);
        case 5,coarse.upper(1,1)=1.01*coarse.upper(1,1);
        case 6,coarse.space.groupVolume(1)=1.01*coarse.space.groupVolume(1);
        case 7,limit=1;
        case 8,coarse.columnLiftScale(1)=1.01*coarse.columnLiftScale(1);
    end
    badFile=fullfile(folder,sprintf('rejected_cache_%02d.mat',kind));save(badFile,'coarse','-v7');
    caught=false;try,loadRocketLeastSquaresCoarseOperator(badFile,cache,space,limit);catch ex,caught=startsWith(ex.identifier,'RocketCoarse:');end
    assert(caught,'CoarseSolverTest:Rejection','错误粗像没有拒绝，kind=%d。',kind);rejected=rejected+1;
end
options.rightPreconditioner='two-level-minres';options.coarseOperatorFile=sourceFile;options.coarseSpaceOptions=space.options;
result=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
assert(result.acceptedNonlinearIterations==1&&result.candidateTrials==1&&result.residualEvaluations==3&& ...
    result.jacobianResidualEvaluations==1&&result.originalNSCombinedDirectionCheckPreserved&&result.fixedLinearPreconditioner);
assert(result.coarseOperatorBuilds==0&&result.coarseOperatorReuses==1&&result.coarseConstructionModelApplications==0&& ...
    result.coarseReuseVerificationModelApplications==1&&result.coarsePreconditionerModelApplications==result.preconditionerApplications);
a=result.linearLog{1};assert(a.fullSystemTargetMet&&a.coarseAlgebraicRelativeResidual<1e-8&&a.fineModelResidualNonincreaseVerified);
limitedOptions=options;limitedOptions.maximumKrylovIterations=1;limited=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,limitedOptions);
assert(limited.candidateTrials==0&&limited.acceptedNonlinearIterations==0&&limited.residualEvaluations==2&& ...
    limited.jacobianResidualEvaluations==1&&isequaln(limited.Q,Q)&&strcmp(limited.stopReason,'steady-newton-linear-target-not-met'));
% 模块四：无保存粗像时从当前Q新建；数值结果相同，仅构建/复用身份与耗时不同。
freshOptions=options;freshOptions=rmfield(freshOptions,'coarseOperatorFile');
fresh=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,freshOptions);
assert(fresh.coarseOperatorBuilds==1&&fresh.coarseOperatorReuses==0&&fresh.coarseConstructionModelApplications==space.coarseUnknowns&& ...
    fresh.coarseReuseVerificationModelApplications==0&&isequaln(fresh.Q,result.Q)&&isequaln(fresh.rate,result.rate));
for kind=1:3
    bad=options;
    switch kind
        case 1,bad.jacobianAction='original-difference';
        case 2,bad.linearEquation='pseudo-transient';
        case 3,bad.rightPreconditioner='lateral-colored-line';
    end
    caught=false;try,solveRocketShellSteadyKrylov(Q,cut,gas,boundary,bad);catch ex,caught=startsWith(ex.identifier,'ShellKrylov:');end
    assert(caught);rejected=rejected+1;
end
report=struct('status','passed-two-level-original-NS-acceptance-rejection-reuse-build-and-old-analytic-result-replay', ...
    'oldAnalyticGraphCompleteResultExactlyReplayed',true,'invalidCacheOrPolicyRejections',rejected, ...
    'coarseUnknowns',result.coarseUnknowns,'actualAcceptedIterations',result.acceptedNonlinearIterations, ...
    'actualSteadyLinearDefect',a.fullSystemRelativeDefect,'selectedAnalyticApplications',result.selectedAnalyticApplications, ...
    'coarsePreconditionerModelApplications',result.coarsePreconditionerModelApplications, ...
    'coarseReuseVerificationModelApplications',result.coarseReuseVerificationModelApplications,'actualNSCalls',result.residualEvaluations, ...
    'limitedActualNSCalls',limited.residualEvaluations,'freshBuildModelApplications',fresh.coarseConstructionModelApplications, ...
    'freshAndReusedActualEndpointIdentical',true,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'steady_newton_unit_evidence.mat'),'report','result','limited','fresh','cut','gas','boundary','reuse','-v7');
fid=fopen(fullfile(folder,'coarse_solver_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function out=untimed(out)
out.history=rmfield(out.history,{'operatorBuildSeconds','linearSolveSeconds'});out=rmfield(out,'analyticCacheBuildSeconds');
end
