function report=testRocketSteadyDirectionGuess(folder)
%TESTROCKETSTEADYDIRECTIONGUESS 校核同Q初猜、修正/跳过路径、总方向原NS及旧完整输出。
% 模块一：原局部策略逐字段复现；四分之一保存方向强制走增量求解。
% 模块二：完整保存方向可跳过GMRES，但仍原NS复核；严格单位目标不放宽。
% 模块三：错误基态、尺度、边界、原面几何、接受身份及配置分别拒绝。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'SteadyGuessTest:Overwrite','不得覆盖方向初猜验证。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','localintegration01','steady_newton_unit_evidence.mat');
s=load(sourceFile,'result','limited','gas','boundary','cut');old=s.result;Q=old.initialQ;cut=s.cut;gas=s.gas;boundary=s.boundary;options=old.options;
baseline=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
assert(isequaln(untimed(baseline),untimed(old)),'SteadyGuessTest:Compatibility','零初猜旧局部路径完整输出没有复现。');
options.steadyInitialGuess=struct('file',sourceFile,'variant','limited','scale',.25);
result=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
save(fullfile(folder,'guess_first_result_evidence.mat'),'result','gas','boundary','cut','-v7');
a=result.linearLog{1};assert(result.acceptedNonlinearIterations==1&&a.initialGuessUsed&&~a.gmresSkipped&&a.initialGuessModelRelativeDefect>.1);
assert(a.effectiveCorrectionTolerance>.1&&a.fullSystemRelativeDefect<=.1&&a.target==.1&&a.meritDescentVerified);
assert(isequaln(a.steadyDirectionEvidence.directionU,a.initialGuessEvidence.directionU+a.initialGuessEvidence.correctionDirectionU));
readyOptions=options;readyOptions.steadyInitialGuess.scale=1;ready=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,readyOptions);
assert(ready.linearLog{1}.gmresSkipped&&ready.acceptedNonlinearIterations==1&&ready.jacobianResidualEvaluations==1&&ready.steadyInitialGuessModelApplications==1);
limitedOptions=options;limitedOptions.maximumKrylovIterations=1;limitedOptions.linearRelativeTolerance=1e-7;
limited=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,limitedOptions);
assert(limited.acceptedNonlinearIterations==0&&limited.candidateTrials==0&&limited.residualEvaluations==2&&isequaln(limited.Q,Q));
rejected=0;original=s.limited;originalBoundary=boundary;strictLimited=limited;
for kind=1:8
    limitedSource=original;badBoundary=originalBoundary;config=options.steadyInitialGuess;
    if kind==1,limitedSource.initialQ(1,1)=limitedSource.initialQ(1,1)*1.01;
    elseif kind==2,limitedSource.residualScale(1)=limitedSource.residualScale(1)*1.01;
    elseif kind==3,badBoundary.solidWall.temperature=badBoundary.solidWall.temperature+1;
    elseif kind==4,limitedSource.linearLog{end}.steadyDirectionEvidence.baseAcceptedIteration=1;
    elseif kind==5
        t=load(limitedSource.options.coarseOperatorFile,'coarse');coarse=t.coarse;coarse.space.geometry.cells.centroid(1,1)=coarse.space.geometry.cells.centroid(1,1)+.01;
        badCoarseFile=fullfile(folder,'bad_guess_geometry.mat');save(badCoarseFile,'coarse','-v7');limitedSource.options.coarseOperatorFile=badCoarseFile;
    elseif kind==6,config.scale=0;
    elseif kind==7,limitedSource.acceptedNonlinearIterations=1;
    else,config.variant='wrong';end
    file=fullfile(folder,sprintf('bad_guess_%02d.mat',kind));limited=limitedSource;boundary=badBoundary;save(file,'limited','gas','boundary','-v7');boundary=originalBoundary;config.file=file;
    caught=false;try,loadRocketSteadyDirectionGuess(config,Q,cut,gas,boundary,result.residualScale);
    catch ex,caught=startsWith(ex.identifier,'SteadyGuess:')||(kind==6&&startsWith(ex.identifier,'MATLAB:'));end
    assert(caught,'SteadyGuessTest:Guard','第%d项非法方向初猜未拒绝。',kind);rejected=rejected+1;
end
limited=strictLimited;
report=struct('status','passed-same-Q-guess-correction-skip-total-NS-gates-source-guards-and-old-full-replay', ...
    'oldLocalCompleteResultReplayed',true,'correctedActualNSDefect',a.fullSystemRelativeDefect,'correctedActualNSCalls',result.residualEvaluations, ...
    'initialGuessModelDefect',a.initialGuessModelRelativeDefect,'effectiveCorrectionTolerance',a.effectiveCorrectionTolerance,'originalTarget',a.target, ...
    'correctedAcceptedIterations',result.acceptedNonlinearIterations,'readySkippedGMRES',ready.linearLog{1}.gmresSkipped, ...
    'readyActualNSCalls',ready.residualEvaluations,'strictLimitedActualNSCalls',limited.residualEvaluations, ...
    'invalidSourceRejections',rejected, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'steady_newton_unit_evidence.mat'),'report','result','ready','limited','gas','boundary','cut','-v7');
fid=fopen(fullfile(folder,'steady_guess_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function out=untimed(out)
out.history=rmfield(out.history,{'operatorBuildSeconds','linearSolveSeconds'});out=rmfield(out,{'analyticCacheBuildSeconds','coarseBuildSeconds','coarseSetupSeconds'});
end
