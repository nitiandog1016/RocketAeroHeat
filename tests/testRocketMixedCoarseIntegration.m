function report=testRocketMixedCoarseIntegration(folder)
%TESTROCKETMIXEDCOARSEINTEGRATION 旧路径完整回放、混合原NS实际验收、来源拒绝和新Q像重建。
% 模块一：混合列仅在显式配置启用；原保守/局部体积粗路径逐字段保持。
% 模块二：缓存与新构造都经原NS组合检查及实际候选保护，不假设完整步必可接受。
% 模块三：错误Q、逐列身份、提升来源和预算明确拒绝；温域和物理网格保持。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'MixedIntegration:Overwrite','不得覆盖混合接入验证。');mkdir(folder);
s=load(fullfile(root,'outputs','v83_validation','targetintegration02','steady_newton_unit_evidence.mat'),'result','cut','gas','boundary');
old=s.result;cut=s.cut;gas=s.gas;boundary=s.boundary;Q=old.initialQ;options=old.options;
baseline=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
assert(isequaln(untimed(baseline),untimed(old)),'MixedIntegration:Compatibility','旧完整局部粗空间结果没有复现。');
sourceFile=fullfile(root,'outputs','v83_validation','mixedunit01','mixed_coarse_evidence.mat');d=load(sourceFile,'coarse');coarse=d.coarse;
options.coarseOperatorFile=sourceFile;options.coarseMixedColumns=coarse.selectedPrimitiveColumns;
result=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
save(fullfile(folder,'mixed_first_result_evidence.mat'),'result','cut','gas','boundary','-v7');
assert(result.acceptedNonlinearIterations==1&&result.jacobianResidualEvaluations==1&&result.originalNSCombinedDirectionCheckPreserved&& ...
    result.linearLog{1}.fullSystemTargetMet&&result.coarseOperatorReuses==1&&result.coarseOperatorBuilds==0);
assert(strcmp(result.coarseBasis,'mixed-conservative-primitive')&&result.coarseUnknowns==147&&result.coarseAdditionalPrimitiveColumns==7&&result.coarseBaseConservativeColumns==140);
limitedOptions=options;limitedOptions.maximumKrylovIterations=1;limitedOptions.linearRelativeTolerance=1e-7;
limited=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,limitedOptions);
assert(limited.acceptedNonlinearIterations==0&&limited.candidateTrials==0&&limited.residualEvaluations==2&&isequaln(limited.Q,Q));
freshOptions=rmfield(options,'coarseOperatorFile');fresh=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,freshOptions);
assert(fresh.coarseOperatorBuilds==1&&fresh.coarseOperatorReuses==0&&fresh.coarseConstructionModelApplications==147&&fresh.acceptedNonlinearIterations==1);
% 新Q采用很小且合法的确定性扰动，单独核对新像；不复用旧Q的数值列。
newQ=Q;newQ(:,5)=newQ(:,5).*(1+1e-5*sin((1:size(Q,1))'));
[~,newInfo]=cutCellNavierStokesResidual(newQ,cut,gas,boundary);
cache=buildCutCellSelectedDirectionCache(newQ,newInfo,cut,gas,boundary,result.residualScale);
parent=buildRocketVolumeCoarseSpace(cut,options.coarseSpaceOptions);space=buildRocketTargetedCoarseSpace(parent,options.coarseTargetRefinement);
rebuilt=rebuildRocketMixedCoarseOperator(cache,space,parent,options.coarseMixedColumns,options.coarseMaximumImageBytes);
weights=sin((1:size(rebuilt.image,2))'*.29);direction=liftRocketMixedCoarse(rebuilt,weights);actual=-applyCutCellSelectedDirection(cache,direction);
rebuildError=norm(actual(:)-rebuilt.image*weights)/max(1,norm(actual(:)));assert(rebuildError<1e-9&&isequaln(rebuilt.baseQ,newQ));
rejected=0;
for kind=1:4
    bad=coarse;limit=options.coarseMaximumImageBytes;ids=options.coarseMixedColumns;
    if kind==1,bad.baseQ(1,1)=bad.baseQ(1,1)*1.01;
    elseif kind==2,bad.liftSourceIdentity.liftRocketMixedCoarse='wrong';
    elseif kind==3,ids=flipud(ids(:));
    else,limit=bad.imageBytes-1;end
    file=fullfile(folder,sprintf('bad_mixed_%02d.mat',kind));coarse=bad;save(file,'coarse','-v7');coarse=d.coarse;
    currentCache=buildCutCellSelectedDirectionCache(Q,old.initialInfo,cut,gas,boundary,result.residualScale);
    caught=false;try,loadRocketMixedCoarseOperator(file,currentCache,space,parent,ids,limit);catch ex,caught=startsWith(ex.identifier,'MixedCoarse:');end
    assert(caught);rejected=rejected+1;
end
report=struct('status','passed-mixed-original-NS-acceptance-strict-rejection-fresh-build-state-rebuild-source-guards-and-old-complete-replay', ...
    'oldTargetedCompleteResultReplayed',true,'coarseUnknowns',result.coarseUnknowns,'actualNSCalls',result.residualEvaluations, ...
    'actualNSLinearDefect',result.linearLog{1}.fullSystemRelativeDefect,'acceptedIterations',result.acceptedNonlinearIterations, ...
    'candidateTrials',result.candidateTrials,'limitedActualNSCalls',limited.residualEvaluations,'limitedUnitLinearTolerance',limitedOptions.linearRelativeTolerance, ...
    'freshConstructionModelApplications',fresh.coarseConstructionModelApplications,'newQConstructionModelApplications',rebuilt.constructionModelApplications, ...
    'newQProbeRelativeDifference',rebuildError,'invalidCacheOrLiftRejections',rejected, ...
    'separateNewStateValidationNSEvaluations',1,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'steady_newton_unit_evidence.mat'),'report','result','limited','fresh','rebuilt','newQ','newInfo','cut','gas','boundary','-v7');
fid=fopen(fullfile(folder,'mixed_integration_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function out=untimed(out)
out.history=rmfield(out.history,{'operatorBuildSeconds','linearSolveSeconds'});out=rmfield(out,{'analyticCacheBuildSeconds','coarseBuildSeconds','coarseSetupSeconds'});
end
