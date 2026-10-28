function report=testRocketSteadyRootReplay(folder,sourceFolder)
%TESTROCKETSTEADYROOTREPLAY 从原Q6重新求根，完整复现Newton记录及独立三步BE确认。
% 模块一：不把保存的根传给Newton；输入为同一真实原始Q6与完整工况。
% 模块二：两算术的完整解/回溯/成本/BE记录逐字段比较，防止以元数据代替求解。
% 模块三：单迭代预算必须返回未收敛，不能把中途停止标成稳态。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v75_validation','root_replay');end
if nargin<2,sourceFolder=fullfile(root,'outputs','v75_validation');end
assert(~exist(folder,'dir'),'RootReplay:Overwrite','不能覆盖真实复算。');rows=struct([]);fresh={};
for name={'double','two_double_euler'}
    old=load(fullfile(sourceFolder,['spatial_root_' name{1}],'spatial_root_evidence.mat'));c=old.c;
    solution=solveAxisymmetricSteadyNewton(old.solution.initialQ,c.mesh,c.gas,c.boundary,old.options);
    assert(isequaln(solution,old.solution)&&solution.spatialRootConverged,'RootReplay:Solution','同输入Newton完整记录无法复现。');
    confirmation=advanceAxisymmetricSteadyLocal(solution.Q,c.mesh,c.gas,c.boundary,old.confirmation.options);
    assert(isequaln(confirmation,old.confirmation)&&confirmation.acceptedSteps==3&&confirmation.steadyConverged, ...
        'RootReplay:Confirmation','原BE后验完整记录无法复现。');
    row=struct('arithmetic',name{1},'newtonIterations',solution.acceptedNewtonIterations,'rootResidual',solution.residualNorm, ...
        'boundaryBalanceRelative',solution.boundaryFluxBalanceRelative,'BEAcceptedSteps',confirmation.acceptedSteps, ...
        'allNewtonAndBEFieldsExactlyReplayed',true);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
    fresh{end+1}=struct('solution',solution,'confirmation',confirmation); %#ok<AGROW>
    if strcmp(name{1},'double')
        limited=old.options;limited.maximumIterations=1;
        stopped=solveAxisymmetricSteadyNewton(old.solution.initialQ,c.mesh,c.gas,c.boundary,limited);
        assert(stopped.acceptedNewtonIterations==1&&~stopped.spatialRootConverged&&strcmp(stopped.stopReason,'maximum-iterations'));
    end
end
report=struct('status','passed-actual-Q6-direct-Newton-and-original-three-BE-exact-replay', ...
    'matlabEngineProcessID',feature('getpid'),'cases',rows,'budgetStopNotPromoted',true, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'root_replay_evidence.mat'),'report','fresh','stopped','-v7');
fid=fopen(fullfile(folder,'root_replay_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;fprintf('%s\n',jsonencode(report));
end
