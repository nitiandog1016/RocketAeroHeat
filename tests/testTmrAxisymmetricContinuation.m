function report=testTmrAxisymmetricContinuation(folder,source)
%TESTTMRAXISYMMETRICCONTINUATION 检查真实检查点续算合同，不重跑初始长段。
% 模块一：保留父文件和尺度，拒绝覆盖、验收放宽和已有累计账损坏。
% 模块二：已稳态检查点仍需新段连续三步；累计账与伪时间分别核对。
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'tests'));
if nargin<1,folder=fullfile(root,'outputs','v63_validation','continuation_contract');end
if nargin<2
    source=fullfile(root,'outputs','v63_validation','tmr_hlle_continue1','tmr_implicit_budget.mat');
    if ~exist(source,'file'),source=fullfile(root,'validation','v63','tmr_hlle_continue1','tmr_implicit_budget.mat');end
end
if ~exist(folder,'dir'),mkdir(folder);end
saved=load(source);assert(saved.implicit.steadyConverged&&isfield(saved,'cumulative'));
refusals=0;
reject(@()studyTmrAxisymmetricContinuation(source,fileparts(source),3),'TmrContinuation:Overwrite');
reject(@()studyTmrAxisymmetricContinuation(source,fullfile(folder,'bad_controls'),3, ...
    struct('steadyResidualTolerance',1)),'TmrContinuation:Controls');
bad=saved;bad.implicit.options.residualScale(1)=2*bad.implicit.residualScale(1);
badFile=fullfile(folder,'bad_scale.mat');save(badFile,'-struct','bad');
reject(@()studyTmrAxisymmetricContinuation(badFile,fullfile(folder,'bad_scale_output'),3),'TmrContinuation:Scale');
bad=saved;bad.cumulative.boundaryIntegral(1)=bad.cumulative.boundaryIntegral(1)+1;
badFile=fullfile(folder,'bad_ledger.mat');save(badFile,'-struct','bad');
reject(@()studyTmrAxisymmetricContinuation(badFile,fullfile(folder,'bad_ledger_output'),3),'TmrContinuation:ParentConservation');
result=studyTmrAxisymmetricContinuation(source,fullfile(folder,'restart'),3);
assert(result.steadyConverged&&result.segmentAcceptedSteps==3&&result.cumulative.balanceRelative<5e-10&& ...
    result.cumulative.acceptedSteps==saved.cumulative.acceptedSteps+3&& ...
    result.diagnosticReplayResidualEvaluations==1&&~result.pseudoTimeIsPhysical&&result.physicalTimeAdvanced==0);
report=struct('status','passed-real-checkpoint-restart-contract-and-refusals','refusals',refusals, ...
    'restartReport',result,'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'continuation_contract_evidence.mat'),'report');
fid=fopen(fullfile(folder,'continuation_contract_summary.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s',jsonencode(report));fclose(fid);
fprintf('Actual TMR continuation contract passed: %d refusals and three-step current-source restart.\n',refusals);
    function reject(action,identifier)
        try,action();catch exception
            assert(strcmp(exception.identifier,identifier),'ContinuationTest:Unexpected','实际拒绝%s，预期%s。',exception.identifier,identifier);
            refusals=refusals+1;return;
        end
        error('ContinuationTest:MissingRefusal','无效检查点合同未拒绝。');
    end
end
