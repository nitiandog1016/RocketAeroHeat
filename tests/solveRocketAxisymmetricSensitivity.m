function report=solveRocketAxisymmetricSensitivity(kind,folder)
%SOLVEROCKETAXISYMMETRICSENSITIVITY 新网格上的完整原空间定常方程实算。
% 模块一：完整保存插值初猜、几何和来源；单次运行不覆盖旧结果。
% 模块二：原2e-6残差/5e-10通量账，不裁状态、不放宽阈值。
% 模块三：成功求根后仍须独立BE确认；失败原样保存，不输出合格云图。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,kind='radial48';end
if nargin<2,folder=fullfile(root,'outputs','v76_validation',kind);end
assert(~exist(folder,'dir'),'RocketSensitivity:Overwrite','禁止覆盖真实敏感性计算。');
c=createRocketAxisymmetricSensitivityCase(struct('kind',kind,'rootEvidenceSHA256',''));
sourceBindings=rocketSensitivitySourceBindings;options=struct('residualScale',c.residualScale, ...
    'steadyResidualTolerance',2e-6,'conservationTolerance',5e-10,'differenceStep',2e-6, ...
    'maximumIterations',24,'maximumLineSearch',14);
mkdir(folder);save(fullfile(folder,'initial_case.mat'),'c','options','sourceBindings','-v7');
clock=tic;solution=solveAxisymmetricSteadyNewton(c.Q,c.mesh,c.gas,c.boundary,options);elapsed=toc(clock);
assert(isequaln(sourceBindings,rocketSensitivitySourceBindings),'RocketSensitivity:Source','运行中求解源码变化。');
report=struct('kind',kind,'status','actual-new-grid-spatial-root-attempt', ...
    'cellCount',c.mesh.cellCount,'originalResidual',solution.residualNorm, ...
    'spatialRootConverged',solution.spatialRootConverged,'stopReason',solution.stopReason, ...
    'acceptedNewtonIterations',solution.acceptedNewtonIterations,'residualEvaluations',solution.residualEvaluations, ...
    'boundaryFluxBalanceRelative',solution.boundaryFluxBalanceRelative,'elapsedSeconds',elapsed, ...
    'engineProcessID',feature('getpid'),'processRecord',getenv('ROCKET_AEROHEAT_PROCESS_RECORD'), ...
    'BEConfirmationPerformed',false,'physicalSolution',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'spatial_root_evidence.mat'),'c','solution','options','sourceBindings','report','-v7');
fid=fopen(fullfile(folder,'spatial_root_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
fprintf('%s\n',jsonencode(report));
end
