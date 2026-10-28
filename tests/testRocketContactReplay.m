function report=testRocketContactReplay(folder,sourceFolder)
%TESTROCKETCONTACTREPLAY 真实重算HLLC粗根及粗根插值后的精根，再真实调用原三步BE。
% 模块一：不把保存新根当Newton输入；完整解/回溯/成本逐字段复现。
% 模块二：细网格供体与本次实际重算粗根一致；失败直接初猜只保留，不改判成功。
% 模块三：Newton和BE分开计数，温域、SA、星态、原残差及守恒标准不改。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v78_validation','root_replay');end
if nargin<2,sourceFolder=fullfile(root,'outputs','v78_validation');end
assert(~exist(folder,'dir'),'RocketContactReplay:Overwrite','不能覆盖真实重算。');rows=struct([]);coarse=[];
kinds={'collar72','nested144_from_hllc72'};
for k=1:2
    old=load(fullfile(sourceFolder,kinds{k},'spatial_root_evidence.mat'));
    if k==1,c=createRocketAxisymmetricContactCase(old.c.options);else,c=createRocketAxisymmetricContactRefinedCase(old.c.options);end
    assert(isequaln(c,old.c));fresh=solveAxisymmetricSteadyNewton(c.Q,c.mesh,c.gas,c.boundary,old.options);
    assert(isequaln(fresh,old.solution)&&fresh.spatialRootConverged,'RocketContactReplay:Root','实际新HLLC方程完整求根不能复现。');
    f=load(fullfile(sourceFolder,kinds{k},'confirmed','rocket_axisymmetric_evidence.mat'),'implicit');o=f.implicit.options;o.maxSteps=3;
    confirmation=advanceAxisymmetricSteadyLocal(fresh.Q,c.mesh,c.gas,c.boundary,o);
    assert(confirmation.acceptedSteps==3&&confirmation.steadyConverged&& ...
        isequal(confirmation.Q,f.implicit.Q)&&isequal(confirmation.history.residual,f.implicit.history.residual)&& ...
        isequal(confirmation.history.dt,f.implicit.history.dt)&&isequal(confirmation.history.newtonTarget,f.implicit.history.newtonTarget)&& ...
        isequal(confirmation.history.newtonDefect,f.implicit.history.newtonDefect),'RocketContactReplay:BE','原三步BE的真实状态/目标/缺陷不能复现。');
    if k==2
        p=load(fullfile(sourceFolder,'collar72','spatial_root_evidence.mat'),'solution');assert(isequal(coarse,p.solution.Q));
    end
    coarse=fresh.Q;
    row=struct('kind',kinds{k},'cellCount',c.mesh.cellCount,'NewtonIterations',fresh.acceptedNewtonIterations, ...
        'residualEvaluations',fresh.residualEvaluations,'originalResidual',fresh.residualNorm, ...
        'completeSolutionExactlyReplayed',true,'actualOriginalBEConfirmationSteps',confirmation.acceptedSteps);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
failed=load(fullfile(sourceFolder,'nested144','spatial_root_evidence.mat'),'solution');
assert(~failed.solution.spatialRootConverged&&strcmp(failed.solution.stopReason,'line-search-failed'));
report=struct('status','passed-two-HLLC-real-initial-guesses-complete-Newton-and-original-BE-replay', ...
    'cases',rows,'failedDirectTrialRecomputed',false,'failedDirectTrialRetainedAsFailure',true, ...
    'engineProcessID',feature('getpid'),'physicalSolution',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'contact_replay_evidence.mat'),'report');
fid=fopen(fullfile(folder,'contact_replay_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
