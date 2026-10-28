function report=runRocketFullMinimumStateComparison(folder)
%RUNROCKETFULLMINIMUMSTATECOMPARISON 完整工程包的真实Q24整箭对照入口。
% 从随包真实检查点重建当前系统，固定三段128步和三个非线性方向预算。
% 运行可能较久；禁止复用旧矩阵资格或在运行途中重启同一目录。
root=fileparts(mfilename('fullpath'));addpath(root,fullfile(root,'tests'),fullfile(root,'internalFlow'));setupCutCellFlow;
if nargin<1,folder=fullfile(root,'outputs',['full_minimum_state_' datestr(now,'yyyymmdd_HHMMSS')]);end
assert(~exist(folder,'dir'),'MinimumState:Overwrite','请使用新目录。');mkdir(folder);
stateFolder=fullfile(root,'outputs/v83_validation/krylov24');geometry=fullfile(root,'outputs/v80_validation/domain02/shell_domain_evidence.mat');
systemFolder=fullfile(folder,'system');prepareRocketCurrentPeakSystem(systemFolder,fullfile(stateFolder,'wall_velocity_krylov_evidence.mat'),geometry,'domain');
preparation=fullfile(systemFolder,'peak_system.json');candidate='';
for k=1:3
    next=fullfile(folder,sprintf('direction%d',k));solveRocketPreparedPeakDirection(next,preparation,struct,candidate);candidate=fullfile(next,'peak_zero_candidate_evidence.mat');
end
report=studyRocketMinimumStateFull(fullfile(folder,'comparison'),preparation,candidate,stateFolder,3);
end
