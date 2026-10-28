function report=testRocketWallResolutionReplay(folder,sourceFolder)
%TESTROCKETWALLRESOLUTIONREPLAY 从两档真实初猜重新求根，完整复现Newton并真实调用原三步BE。
% 模块一：未把保存的新根送进Newton；144层父72胶囊必须与本次重算根逐值一致。
% 模块二：准备Newton与BE确认分别计数，原残差/通量账/BE门槛保持。
% 模块三：保存真实复算摘要，不将包内重跑等同于独立物理或实验验证。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v77_validation','root_replay');end
if nargin<2,sourceFolder=fullfile(root,'outputs','v77_validation');end
assert(~exist(folder,'dir'),'RocketWallReplay:Overwrite','禁止覆盖真实复算。');rows=struct([]);previous=[];
for kinds={'collar72','nested144'}
    kind=kinds{1};old=load(fullfile(sourceFolder,kind,'spatial_root_evidence.mat'));
    c=createRocketAxisymmetricWallResolutionCase(old.c.options);assert(isequaln(c,old.c));
    fresh=solveAxisymmetricSteadyNewton(c.Q,c.mesh,c.gas,c.boundary,old.options);
    assert(isequaln(fresh,old.solution)&&fresh.spatialRootConverged,'RocketWallReplay:Root','完整Newton解/候选/成本不能复现。');
    frozen=load(fullfile(sourceFolder,kind,'confirmed','rocket_axisymmetric_evidence.mat'),'implicit');o=frozen.implicit.options;o.maxSteps=3;
    confirmation=advanceAxisymmetricSteadyLocal(fresh.Q,c.mesh,c.gas,c.boundary,o);
    assert(confirmation.acceptedSteps==3&&confirmation.steadyConverged&&all(confirmation.history.dt(2:end)>0)&& ...
        isequal(confirmation.Q,frozen.implicit.Q)&&isequal(confirmation.history.residual,frozen.implicit.history.residual)&& ...
        isequal(confirmation.history.dt,frozen.implicit.history.dt)&& ...
        isequal(confirmation.history.newtonTarget,frozen.implicit.history.newtonTarget)&& ...
        isequal(confirmation.history.newtonDefect,frozen.implicit.history.newtonDefect), ...
        'RocketWallReplay:BE','真实原三步BE确认不能复现冻结端点与原目标/缺陷。');
    if strcmp(kind,'nested144')
        parent=load(fullfile(sourceFolder,'collar72','spatial_root_evidence.mat'),'solution');
        assert(isequal(previous,parent.solution.Q),'RocketWallReplay:Parent','本次72层重算根与实际144层内部供体不同。');
    end
    previous=fresh.Q;
    row=struct('kind',kind,'cellCount',c.mesh.cellCount,'NewtonIterations',fresh.acceptedNewtonIterations, ...
        'residualEvaluations',fresh.residualEvaluations,'originalResidual',fresh.residualNorm, ...
        'completeSolutionExactlyReplayed',true,'actualOriginalBEConfirmationSteps',confirmation.acceptedSteps);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
report=struct('status','passed-two-real-initial-guesses-complete-Newton-and-original-BE-exact-replay', ...
    'cases',rows,'engineProcessID',feature('getpid'),'physicalSolution',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'wall_resolution_replay_evidence.mat'),'report');
fid=fopen(fullfile(folder,'wall_resolution_replay_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
