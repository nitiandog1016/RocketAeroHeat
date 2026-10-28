function report=testRocketShellSelectedDirection(folder)
%TESTROCKETSHELLSELECTEDDIRECTION 解析GMRES模型接入后仍用原NS复核组合方向与候选。
% 模块一：先逐字段复现旧直接稳态缺省，再实际测试新模式与受限预算拒绝。
% 模块二：解析应用次数和原NS调用分开，不能将廉价模型作用记作原NS调用。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'SelectedDirectionTest:Overwrite','不得覆盖解析作用接入测试。');mkdir(folder);
s=load(fullfile(root,'outputs','v83_validation','newtonunit03','steady_newton_unit_evidence.mat'),'result','cut','gas','boundary');
Q=s.result.initialQ;cut=s.cut;gas=s.gas;boundary=s.boundary;options=s.result.options;
options.jacobianAction='original-difference';normal=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
assert(isequaln(untimed(normal),untimed(s.result)),'SelectedDirectionTest:Compatibility','显式原差分没有复现完整旧直接稳态结果。');
options.jacobianAction='selected-analytic';options.rightPreconditioner='lateral-colored-line';options.graphSweeps=1;
result=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
assert(result.acceptedNonlinearIterations==1&&result.selectedAnalyticApplications>0&&result.originalNSCombinedDirectionCheckPreserved);
e=result.linearLog{1}.steadyDirectionEvidence;
assert(strcmp(e.scheme,'forward')&&e.actualOriginalResidualEvaluations==1&&result.linearLog{1}.actualCombinedDirectionCheckedWithOriginalNS);
assert(result.jacobianResidualEvaluations==1&&result.residualEvaluations==3&&result.candidateTrials==1);
limitedOptions=options;limitedOptions.maximumKrylovIterations=1;limited=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,limitedOptions);
assert(limited.acceptedNonlinearIterations==0&&limited.candidateTrials==0&&isequaln(limited.Q,Q)&& ...
    strcmp(limited.stopReason,'steady-newton-linear-target-not-met')&&limited.jacobianResidualEvaluations==1&&limited.residualEvaluations==2);
report=struct('status','passed-selected-analytic-GMRES-model-original-NS-combined-check-and-budget-rejection', ...
    'acceptedIterations',result.acceptedNonlinearIterations,'selectedAnalyticApplications',result.selectedAnalyticApplications, ...
    'actualOriginalResidualEvaluations',result.residualEvaluations,'originalDirectionalRechecks',result.jacobianResidualEvaluations, ...
    'actualSteadyLinearDefect',result.linearLog{1}.fullSystemRelativeDefect,'limitedOriginalResidualEvaluations',limited.residualEvaluations, ...
    'originalDifferenceCompleteResultReplayed',true,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'steady_newton_unit_evidence.mat'),'report','result','limited','cut','gas','boundary','-v7');
fid=fopen(fullfile(folder,'selected_direction_solver_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function out=untimed(out)
out.history=rmfield(out.history,{'operatorBuildSeconds','linearSolveSeconds'});
end
