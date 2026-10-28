function report=testRocketLocalTrialReplay(folder,source)
%TESTROCKETLOCALTRIALREPLAY 从冻结同输入重新计算真实完整BE，逐值复现全部Newton候选。
% 模块一：仅使用传入完整工况/Q/选项，不把保存的末态传给求解器。
% 模块二：完整result与候选逐字段一致；时间/物理精度不由此认证。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v74_validation','trial_replay');end
if nargin<2,source=fullfile(root,'outputs','v74_validation','local_difference_2em06','local_difference_evidence.mat');end
assert(~exist(folder,'dir'),'LocalTrialReplay:Overwrite','不能覆盖真实复算证据。');
old=load(source,'c','Q','options','trial','diagnostics');c=old.c;
[trial,diagnostics]=advanceAxisymmetricSteadyLocal(old.Q,c.mesh,c.gas,c.boundary,old.options, ...
    struct('enabled',true,'maxRecords',1024));
assert(isequaln(trial,old.trial)&&isequaln(diagnostics,old.diagnostics), ...
    'LocalTrialReplay:Different','同输入完整BE结果或真实候选无法逐字段复现。');
assert(trial.acceptedSteps==1&&trial.rejectedAttempts==0&&diagnostics.complete);
report=struct('status','passed-actual-local-full-BE-and-all-Newton-candidates-exact-replay', ...
    'matlabEngineProcessID',feature('getpid'),'acceptedSteps',trial.acceptedSteps, ...
    'candidateCount',diagnostics.seenRecords,'residualEvaluations',trial.residualEvaluations, ...
    'acceptedDtS',trial.history.dt(end),'residual',trial.residualNorm,'fullResultExactlyEqual',true, ...
    'fullDiagnosticsExactlyEqual',true,'physicalSolution',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'trial_replay_evidence.mat'),'report','trial','diagnostics','-v7');
fid=fopen(fullfile(folder,'trial_replay_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;fprintf('%s\n',jsonencode(report));
end
