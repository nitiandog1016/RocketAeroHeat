function report=testRocketRecycledIntegration(folder)
%TESTROCKETRECYCLEDINTEGRATION 显式多方向增强接入、旧单方向完整回放及原NS真实接受/拒绝。
% 模块一：保留既有初猜/粗像合同；默认路径不增加新字段或改变数值。
% 模块二：新模式直接输出物理增量，总方向由原NS重算，严格目标不绕过。
% 模块三：基存储预算拒绝独立检查，单位结果不拼入真实主链。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'RecycledIntegration:Overwrite','不得覆盖增强接入测试。');mkdir(folder);
s=load(fullfile(root,'outputs','v83_validation','guessintegration01','steady_newton_unit_evidence.mat'),'result','gas','boundary','cut');
old=s.result;Q=old.initialQ;cut=s.cut;gas=s.gas;boundary=s.boundary;options=old.options;
baseline=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
assert(isequaln(untimed(old),untimed(baseline)),'RecycledIntegration:Compatibility','旧单方向初猜完整输出没有复现。');
options.maximumKrylovIterations=24;options.recycledKrylov=struct('restartLength',8,'maximumCycles',3,'retainedDirections',3, ...
    'maximumBasisBytes',64*1024^2,'saveBasisEvidence',true);
result=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
save(fullfile(folder,'recycled_first_result_evidence.mat'),'result','gas','boundary','cut','-v7');
assert(result.acceptedNonlinearIterations==1&&result.originalNSCombinedDirectionCheckPreserved&&result.linearLog{1}.fullSystemTargetMet);
assert(result.recycledKrylovArnoldiSteps>0&&result.recycledKrylovModelApplications==result.linearLog{1}.recycledKrylovReport.modelApplications);
limitedOptions=options;limitedOptions.maximumKrylovIterations=1;limitedOptions.linearRelativeTolerance=1e-7;
limited=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,limitedOptions);
assert(limited.acceptedNonlinearIterations==0&&limited.candidateTrials==0&&limited.residualEvaluations==2&&isequaln(limited.Q,Q));
bad=options;bad.recycledKrylov.maximumBasisBytes=1;rejected=false;
try,solveRocketShellSteadyKrylov(Q,cut,gas,boundary,bad);catch ex,rejected=strcmp(ex.identifier,'RecycledGMRES:Budget');end
assert(rejected,'RecycledIntegration:Budget','主链未明确拒绝基存储预算。');
report=struct('status','passed-explicit-SVD-recycled-main-integration-original-NS-gates-and-old-guess-full-replay', ...
    'oldSingleGuessCompleteResultReplayed',true,'actualNSDefect',result.linearLog{1}.fullSystemRelativeDefect, ...
    'actualNSCalls',result.residualEvaluations,'acceptedIterations',result.acceptedNonlinearIterations, ...
    'arnoldiSteps',result.recycledKrylovArnoldiSteps,'recycledModelApplications',result.recycledKrylovModelApplications, ...
    'limitedActualNSCalls',limited.residualEvaluations,'basisBudgetRejected',rejected,'failedBudgetBaselineNSCalls',1, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'steady_newton_unit_evidence.mat'),'report','result','limited','gas','boundary','cut','-v7');
fid=fopen(fullfile(folder,'recycled_integration_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function out=untimed(out)
out.history=rmfield(out.history,{'operatorBuildSeconds','linearSolveSeconds'});out=rmfield(out,{'analyticCacheBuildSeconds','coarseBuildSeconds','coarseSetupSeconds'});
end
