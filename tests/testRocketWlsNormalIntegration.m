function report=testRocketWlsNormalIntegration(folder)
%TESTROCKETWLSNORMALINTEGRATION 核对显式五带接入、旧完整输出及原NS实际接受/拒绝。
% 模块一：旧同Q初猜路径逐字段回放，五带配置只在显式选择时生效。
% 模块二：总方向继续用原NS及原norm(R)门槛，不能只验增量。
% 模块三：预算/模式拒绝与单位严格目标拒绝分列，不加入真实主链。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'WlsNormalIntegration:Overwrite','不得覆盖五带接入验证。');mkdir(folder);
s=load(fullfile(root,'outputs','v83_validation','guessintegration01','steady_newton_unit_evidence.mat'),'result','gas','boundary','cut');
old=s.result;Q=old.initialQ;cut=s.cut;gas=s.gas;boundary=s.boundary;options=old.options;
baseline=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
assert(isequaln(untimed(old),untimed(baseline)),'WlsNormalIntegration:Compatibility','旧同Q初猜完整输出没有复现。');
options.maximumKrylovIterations=24;options.wlsNormalCorrection=struct('maximumAdditionalBytes',32*1024^2);
result=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
save(fullfile(folder,'wls_normal_first_result.mat'),'result','cut','gas','boundary','-v7');
assert(result.acceptedNonlinearIterations==1&&result.originalNSCombinedDirectionCheckPreserved&&result.linearLog{1}.fullSystemTargetMet);
assert(result.wlsNormalCorrection.stressHeatAndBothWorkTermsIncluded&&~result.wlsNormalCorrection.crossOriginalFacetWlsBlocksIncluded);
limitedOptions=options;limitedOptions.maximumKrylovIterations=1;limitedOptions.linearRelativeTolerance=1e-7;
limited=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,limitedOptions);
assert(limited.acceptedNonlinearIterations==0&&limited.candidateTrials==0&&limited.residualEvaluations==2&&isequaln(limited.Q,Q));
bad=options;bad.wlsNormalCorrection.maximumAdditionalBytes=1;rejected=false;
try,solveRocketShellSteadyKrylov(Q,cut,gas,boundary,bad);catch ex,rejected=strcmp(ex.identifier,'WlsNormal:Budget');end;assert(rejected);
bad=options;bad.rightPreconditioner='lateral-colored-line';rejectedMode=false;
try,solveRocketShellSteadyKrylov(Q,cut,gas,boundary,bad);catch ex,rejectedMode=startsWith(ex.identifier,'ShellKrylov:')||strcmp(ex.identifier,'WlsNormal:Mode');end;assert(rejectedMode);
report=struct('status','passed-explicit-WLS-normal-five-band-main-integration-old-guess-full-replay-and-original-NS-gates', ...
    'oldSingleGuessCompleteOutputReplayed',true,'actualNSDefect',result.linearLog{1}.fullSystemRelativeDefect,'actualNSCalls',result.residualEvaluations, ...
    'acceptedIterations',result.acceptedNonlinearIterations,'limitedActualNSCalls',limited.residualEvaluations, ...
    'additionalBlockBudgetRejected',rejected,'invalidModeRejected',rejectedMode,'budgetFailureBaselineNSCalls',1, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'steady_newton_unit_evidence.mat'),'report','result','limited','gas','boundary','cut','-v7');
fid=fopen(fullfile(folder,'wls_normal_integration_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function out=untimed(out)
out.history=rmfield(out.history,{'operatorBuildSeconds','linearSolveSeconds'});out=rmfield(out,{'analyticCacheBuildSeconds','coarseBuildSeconds','coarseSetupSeconds'});
end
