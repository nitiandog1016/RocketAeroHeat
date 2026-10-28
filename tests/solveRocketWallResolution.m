function report=solveRocketWallResolution(kind,parent72Hash)
%SOLVEROCKETWALLRESOLUTION 有界近壁网格上的真实五量定常Newton实算。
% 模块一：新工况/插值与源SHA前锁；各运行独立目录不覆盖旧解。
% 模块二：原2e-6与5e-10平衡门槛保持；实际失败完整保留，不补状态。
% 模块三：返回真实根胶囊SHA供下一档内部初猜，不把数值根当边界或BE步。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,kind='collar72';end;if nargin<2,parent72Hash='';end
folder=fullfile(root,'outputs','v77_validation',kind);assert(~exist(folder,'dir'),'RocketWallResolution:Overwrite','禁止覆盖真实求解。');
c=createRocketAxisymmetricWallResolutionCase(struct('kind',kind,'parent72SHA256',parent72Hash,'rootEvidenceSHA256',''));
sourceBindings=rocketWallResolutionSourceBindings;options=struct('residualScale',c.residualScale, ...
    'steadyResidualTolerance',2e-6,'conservationTolerance',5e-10,'differenceStep',2e-6, ...
    'maximumIterations',24,'maximumLineSearch',14);
mkdir(folder);save(fullfile(folder,'initial_case.mat'),'c','sourceBindings','options','-v7');
clock=tic;solution=solveAxisymmetricSteadyNewton(c.Q,c.mesh,c.gas,c.boundary,options);elapsed=toc(clock);
assert(isequaln(sourceBindings,rocketWallResolutionSourceBindings),'RocketWallResolution:Source','实际求解源码变化。');
report=struct('kind',kind,'cellCount',c.mesh.cellCount,'spatialRootConverged',solution.spatialRootConverged, ...
    'originalResidual',solution.residualNorm,'stopReason',solution.stopReason, ...
    'acceptedNewtonIterations',solution.acceptedNewtonIterations,'residualEvaluations',solution.residualEvaluations, ...
    'boundaryFluxBalanceRelative',solution.boundaryFluxBalanceRelative,'elapsedSeconds',elapsed, ...
    'engineProcessID',feature('getpid'),'processRecord',getenv('ROCKET_AEROHEAT_PROCESS_RECORD'), ...
    'BEConfirmationPerformed',false,'physicalSolution',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
target=fullfile(folder,'spatial_root_evidence.mat');save(target,'c','solution','options','sourceBindings','report','-v7');
fid=fopen(target,'rb');guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);report.rootEvidenceSHA256=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
fid=fopen(fullfile(folder,'spatial_root_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
fprintf('%s\n',jsonencode(report));
end
