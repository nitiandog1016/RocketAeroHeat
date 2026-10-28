function report=testRocketTargetedCoarseIntegration(folder)
%TESTROCKETTARGETEDCOARSEINTEGRATION 真实NS验收、旧完整粗路径回放及局部细分缓存/新Q重建。
% 模块一：旧90/80类体积粗校正完整输出保持；新选项只在显式两层模式启用。
% 模块二：原NS检查方向和候选，严格单位预算拒绝，不以粗层小残差替代。
% 模块三：新Q重建粗像；错空间/预算必须拒绝，缓存与新构造取得一致端点。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'TargetIntegration:Overwrite','不得覆盖局部细分接入测试。');mkdir(folder);
s=load(fullfile(root,'outputs','v83_validation','coarsesolver01','steady_newton_unit_evidence.mat'),'result','cut','gas','boundary');
old=s.result;cut=s.cut;gas=s.gas;boundary=s.boundary;Q=old.initialQ;options=old.options;
baseline=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
assert(isequaln(untimed(baseline),untimed(old)),'TargetIntegration:Compatibility','旧完整体积粗路径没有复现。');
sourceFile=fullfile(root,'outputs','v83_validation','basiscandidatesunit01','spatial','coarse_candidate_evidence.mat');d=load(sourceFile,'coarse');
options.coarseSpaceOptions=d.coarse.space.options;options.coarseTargetRefinement=d.coarse.space.targetRefinement.options;options.coarseOperatorFile=sourceFile;
result=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
save(fullfile(folder,'targeted_first_result_evidence.mat'),'result','cut','gas','boundary','-v7');
assert(result.acceptedNonlinearIterations==1&&result.residualEvaluations==3&&result.jacobianResidualEvaluations==1&&result.candidateTrials==3);
assert(all(~[result.lineLog(1:2).accepted])&&all(strcmp({result.lineLog(1:2).errorIdentifier},'CutNS:State'))&& ...
    result.lineLog(3).accepted&&result.lineLog(3).alpha==.25,'TargetIntegration:Safeguard','原状态拒绝/缩步过程没有正确保留。');
assert(strcmp(result.coarseBasis,'conservative-targeted-spatial')&&result.parentCoarseUnknowns==120&&result.coarseUnknowns==140&&result.coarseTargetCellCount>0);
assert(result.coarseOperatorReuses==1&&result.coarseOperatorBuilds==0&&result.coarseConstructionModelApplications==0&&result.coarsePreconditionerModelApplications==result.preconditionerApplications);
limitedOptions=options;limitedOptions.maximumKrylovIterations=1;limitedOptions.linearRelativeTolerance=1e-6;
limited=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,limitedOptions);
assert(limited.acceptedNonlinearIterations==0&&limited.candidateTrials==0&&limited.residualEvaluations==2&&isequaln(limited.Q,Q));
rebuildOptions=options;rebuildOptions.maximumIterations=2;rebuilt=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,rebuildOptions);
assert(rebuilt.coarseOperatorReuses==1&&rebuilt.coarseOperatorBuilds==1&&rebuilt.coarseConstructionModelApplications==140&& ...
    rebuilt.analyticCacheBuilds==2&&rebuilt.acceptedNonlinearIterations>=1,'TargetIntegration:Rebuild','接受新Q后没有重建当前粗像。');
rejected=0;bad=options;bad.coarseTargetRefinement.maximumCoarseUnknowns=5;
try,solveRocketShellSteadyKrylov(Q,cut,gas,boundary,bad);catch ex,assert(strcmp(ex.identifier,'TargetCoarse:Budget'));rejected=rejected+1;end
bad=options;bad.coarseOperatorFile=old.options.coarseOperatorFile;
try,solveRocketShellSteadyKrylov(Q,cut,gas,boundary,bad);catch ex,assert(strcmp(ex.identifier,'RocketCoarse:Cache'));rejected=rejected+1;end
assert(rejected==2);
report=struct('status','passed-targeted-coarse-original-NS-acceptance-strict-budget-rejection-new-Q-rebuild-and-old-full-output-replay', ...
    'oldCompleteCoarseResultReplayed',true,'coarseUnknowns',result.coarseUnknowns,'acceptedIterations',result.acceptedNonlinearIterations, ...
    'actualNSLinearDefect',result.linearLog{1}.fullSystemRelativeDefect,'actualNSCalls',result.residualEvaluations, ...
    'candidateTrials',result.candidateTrials,'invalidStateCandidatesRejected',2,'acceptedAlpha',result.lineLog(end).alpha, ...
    'limitedActualNSCalls',limited.residualEvaluations,'limitedUnitLinearTolerance',limitedOptions.linearRelativeTolerance, ...
    'newQRebuiltModelApplications',rebuilt.coarseConstructionModelApplications,'newQAcceptedIterations',rebuilt.acceptedNonlinearIterations, ...
    'invalidBudgetOrSpaceRejections',rejected,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'steady_newton_unit_evidence.mat'),'report','result','limited','rebuilt','cut','gas','boundary','-v7');
fid=fopen(fullfile(folder,'targeted_integration_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function out=untimed(out)
out.history=rmfield(out.history,{'operatorBuildSeconds','linearSolveSeconds'});out=rmfield(out,{'analyticCacheBuildSeconds','coarseBuildSeconds','coarseSetupSeconds'});
end
