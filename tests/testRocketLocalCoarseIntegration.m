function report=testRocketLocalCoarseIntegration(folder)
%TESTROCKETLOCALCOARSEINTEGRATION 验证显式局部策略、原NS接受/拒绝、新Q重建及旧完整输出。
% 模块一：旧混合路径逐字段回放，新增配置只在显式选择时生效。
% 模块二：缓存与新构造均核对原NS组合方向和实际候选，不假定完整步必被接受。
% 模块三：嵌套Q、原CV身份、面映射、来源、因子及存储预算分别拒绝。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'LocalIntegration:Overwrite','不得覆盖局部接入验证。');mkdir(folder);
s=load(fullfile(root,'outputs','v83_validation','mixedintegration01','steady_newton_unit_evidence.mat'),'result','cut','gas','boundary');
old=s.result;cut=s.cut;gas=s.gas;boundary=s.boundary;Q=old.initialQ;options=old.options;
baseline=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
assert(isequaln(untimed(baseline),untimed(old)),'LocalIntegration:Compatibility','旧混合空间完整输出未逐字段复现。');
sourceFile=fullfile(root,'outputs','v83_validation','localunit02','local_coarse_evidence.mat');s=load(sourceFile,'coarse');coarse=s.coarse;
columns=struct('retainedMixedColumns',coarse.retainedMixedColumns,'cellIDs',int64(coarse.localCellIDs),'equationIDs',int64(coarse.localEquationIDs));
options.coarseOperatorFile=sourceFile;options.coarseLocalColumns=columns;
result=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
save(fullfile(folder,'local_first_result_evidence.mat'),'result','cut','gas','boundary','-v7');
assert(result.acceptedNonlinearIterations==1&&result.coarseOperatorReuses==1&&result.coarseOperatorBuilds==0&& ...
    result.coarseConstructionModelApplications==0&&result.originalNSCombinedDirectionCheckPreserved&&result.linearLog{1}.fullSystemTargetMet);
assert(strcmp(result.coarseBasis,'local-conservative-enriched-mixed')&&result.coarseUnknowns==149&& ...
    result.coarseRetainedMixedColumns==147&&result.coarseLocalConservativeColumns==2);
limitedOptions=options;limitedOptions.maximumKrylovIterations=1;limitedOptions.linearRelativeTolerance=1e-7;
limited=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,limitedOptions);
assert(limited.acceptedNonlinearIterations==0&&limited.candidateTrials==0&&limited.residualEvaluations==2&&isequaln(limited.Q,Q));
freshOptions=rmfield(options,'coarseOperatorFile');fresh=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,freshOptions);
assert(fresh.coarseOperatorBuilds==1&&fresh.coarseOperatorReuses==0&&fresh.coarseConstructionModelApplications==149&&fresh.acceptedNonlinearIterations==1);
parent=buildRocketVolumeCoarseSpace(cut,options.coarseSpaceOptions);space=buildRocketTargetedCoarseSpace(parent,options.coarseTargetRefinement);
newQ=Q;newQ(:,5)=newQ(:,5).*(1+1e-5*sin((1:size(Q,1))'));
[~,newInfo]=cutCellNavierStokesResidual(newQ,cut,gas,boundary);newCache=buildCutCellSelectedDirectionCache(newQ,newInfo,cut,gas,boundary,result.residualScale);
rebuilt=rebuildRocketLocalCoarseOperator(newCache,space,parent,options.coarseMixedColumns,columns,options.coarseMaximumImageBytes);
weights=sin((1:size(rebuilt.image,2))'*.29);actual=-applyCutCellSelectedDirection(newCache,liftRocketLocalCoarse(rebuilt,weights));
rebuildError=norm(actual(:)-rebuilt.image*weights)/max(1,norm(actual(:)));assert(rebuildError<1e-9&&isequaln(rebuilt.baseQ,newQ));
currentCache=buildCutCellSelectedDirectionCache(Q,old.initialInfo,cut,gas,boundary,result.residualScale);original=coarse;rejected=0;
for kind=1:8
    coarse=original;limit=options.coarseMaximumImageBytes;ids=columns;
    if kind==1,coarse.baseQ(1,1)=coarse.baseQ(1,1)*1.01;
    elseif kind==2,coarse.mixedLift.baseQ(1,1)=coarse.mixedLift.baseQ(1,1)*1.01;
    elseif kind==3,ids.cellIDs=flipud(ids.cellIDs);ids.equationIDs=flipud(ids.equationIDs);
    elseif kind==4,coarse.localOriginalFaceIDs(1)=coarse.localOriginalFaceIDs(1)+1;
    elseif kind==5,coarse.liftSourceIdentity.liftRocketLocalCoarse='wrong';
    elseif kind==6,coarse.upper(1,1)=coarse.upper(1,1)*1.01;
    elseif kind==7,limit=coarse.imageBytes-1;
    else,ids.cellIDs(1)=size(Q,1)+1;end
    file=fullfile(folder,sprintf('bad_local_%02d.mat',kind));save(file,'coarse','-v7');caught=false;
    try,loadRocketLocalCoarseOperator(file,currentCache,space,parent,options.coarseMixedColumns,ids,limit);
    catch ex,caught=startsWith(ex.identifier,'LocalCoarse:');end
    assert(caught,'LocalIntegration:Guard','第%d项非法缓存没有明确拒绝。',kind);rejected=rejected+1;
end
coarse=original;
report=struct('status','passed-local-cache-fresh-current-state-rebuild-original-NS-gates-and-old-mixed-complete-replay', ...
    'oldMixedCompleteResultReplayed',true,'coarseUnknowns',result.coarseUnknowns,'actualNSCalls',result.residualEvaluations, ...
    'actualNSLinearDefect',result.linearLog{1}.fullSystemRelativeDefect,'acceptedIterations',result.acceptedNonlinearIterations, ...
    'candidateTrials',result.candidateTrials,'limitedActualNSCalls',limited.residualEvaluations,'limitedUnitLinearTolerance',1e-7, ...
    'freshConstructionModelApplications',fresh.coarseConstructionModelApplications,'newQConstructionModelApplications',rebuilt.constructionModelApplications, ...
    'newQProbeRelativeDifference',rebuildError,'invalidCacheOrLiftRejections',rejected,'integerConfigurationInputTested',true, ...
    'separateNewStateValidationNSEvaluations',1,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'steady_newton_unit_evidence.mat'),'report','result','limited','fresh','rebuilt','newQ','newInfo','coarse','cut','gas','boundary','-v7');
fid=fopen(fullfile(folder,'local_integration_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function out=untimed(out)
out.history=rmfield(out.history,{'operatorBuildSeconds','linearSolveSeconds'});out=rmfield(out,{'analyticCacheBuildSeconds','coarseBuildSeconds','coarseSetupSeconds'});
end
