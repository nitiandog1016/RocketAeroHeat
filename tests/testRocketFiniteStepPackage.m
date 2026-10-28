function report=testRocketFiniteStepPackage(folder)
%TESTROCKETFINITESTEPPACKAGE 在独立工程目录验证本包入口、夹具与整箭原NS基态。
% 不把一次全域基态回放冒充全部整箭对照的隔离重算或物理验证。
root=fileparts(fileparts(mfilename('fullpath')));assert(contains(root,'isolated'),'FiniteStep:Isolation','必须在隔离目录运行。');
restoredefaultpath;addpath(root,fullfile(root,'tests'),fullfile(root,'internalFlow'));setupCutCellFlow;
assert(startsWith(which('cutCellNavierStokesResidual'),root));
local=runRocketFiniteStepResearch(folder);timer=tic;
s=load(fullfile(root,'outputs/v83_validation/krylov24/wall_velocity_krylov_evidence.mat'),'result','gas','boundary');
d=load(fullfile(root,'outputs/v80_validation/domain02/shell_domain_evidence.mat'),'domain');
[rate,info]=cutCellNavierStokesResidual(s.result.Q,d.domain,s.gas,s.boundary);
replay=assertRocketNSReplay(rate,info,s.result.rate,s.result.info,d.domain.cells.volume);
report=struct('status','passed-isolated-package-local-six-step-and-complete-Q24-NS-replay', ...
    'local',local,'completeQ24Replay',replay,'completeQ24ReplaySeconds',toc(timer), ...
    'fullMultistepComparisonRecomputedInIsolation',false,'resolvedNSFile',which('cutCellNavierStokesResidual'), ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
fid=fopen(fullfile(folder,'isolated_package.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(report));fclose(fid);
end
