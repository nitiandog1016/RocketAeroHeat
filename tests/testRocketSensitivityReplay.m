function report=testRocketSensitivityReplay(folder,sourceFolder)
%TESTROCKETSENSITIVITYREPLAY 从三档实际插值初猜重新求根，逐字段复现完整Newton记录。
% 模块一：不把已保存新根传回求解器，真实重算各新网格五量方程。
% 模块二：源工况、解、回溯/成本完全一致，不能以重读摘要充当实际复算。
% 模块三：这里不是第二套RANS或实验验证；BE不可变链另外逐步验算。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v76_validation','root_replay');end
if nargin<2,sourceFolder=fullfile(root,'outputs','v76_validation');end
assert(~exist(folder,'dir'),'SensitivityReplay:Overwrite','禁止覆盖真实复算。');rows=struct([]);
for kinds={'radial48','nearwall_half','farfield1p5'}
    kind=kinds{1};old=load(fullfile(sourceFolder,kind,'spatial_root_evidence.mat'));
    c=createRocketAxisymmetricSensitivityCase(old.c.options);assert(isequaln(c,old.c));
    fresh=solveAxisymmetricSteadyNewton(c.Q,c.mesh,c.gas,c.boundary,old.options);
    assert(isequaln(fresh,old.solution)&&fresh.spatialRootConverged, ...
        'SensitivityReplay:Solution','三档各自完整Newton解/候选/成本不能逐字段复现。');
    row=struct('kind',kind,'cellCount',c.mesh.cellCount,'NewtonIterations',fresh.acceptedNewtonIterations, ...
        'residualEvaluations',fresh.residualEvaluations,'originalResidual',fresh.residualNorm,'completeSolutionExactlyReplayed',true);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
report=struct('status','passed-three-actual-interpolated-initial-guesses-new-grid-Newton-exact-replay', ...
    'cases',rows,'engineProcessID',feature('getpid'),'physicalSolution',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'sensitivity_replay_evidence.mat'),'report');
fid=fopen(fullfile(folder,'sensitivity_replay_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
fprintf('%s\n',jsonencode(report));
end
