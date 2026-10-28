function report=solveRocketContactStudy(kind)
%SOLVEROCKETCONTACTSTUDY 从同网格有效HLLE根真实求原HLLC黏性定常方程。
% 模块一：同一物理输入，只更改已有内部通量；原星态护栏不删，不隐式回退。
% 模块二：原2e-6空间、5e-10通量账及完整五量温域/SA约束保持。
% 模块三：成功/未收敛均保存完整初末、候选、成本及实际源码，不冒称物理资格。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,kind='collar72';end
folder=fullfile(root,'outputs','v78_validation',kind);assert(~exist(folder,'dir'),'RocketContact:Overwrite','不能覆盖真实受控试验。');
c=createRocketAxisymmetricContactCase(struct('kind',kind,'rootEvidenceSHA256',''));
sourceBindings=rocketContactSourceBindings;options=struct('residualScale',c.residualScale, ...
    'steadyResidualTolerance',2e-6,'conservationTolerance',5e-10,'differenceStep',2e-6,'maximumIterations',24,'maximumLineSearch',14);
mkdir(folder);save(fullfile(folder,'initial_case.mat'),'c','sourceBindings','options','-v7');
clock=tic;solution=solveAxisymmetricSteadyNewton(c.Q,c.mesh,c.gas,c.boundary,options);elapsed=toc(clock);
assert(isequaln(sourceBindings,rocketContactSourceBindings),'RocketContact:Source','实际求解源码变化。');
report=struct('kind',kind,'cellCount',c.mesh.cellCount,'interiorConvectiveFlux','hllc','spatialRootConverged',solution.spatialRootConverged, ...
    'originalResidual',solution.residualNorm,'stopReason',solution.stopReason,'acceptedNewtonIterations',solution.acceptedNewtonIterations, ...
    'residualEvaluations',solution.residualEvaluations,'boundaryFluxBalanceRelative',solution.boundaryFluxBalanceRelative, ...
    'elapsedSeconds',elapsed,'engineProcessID',feature('getpid'),'processRecord',getenv('ROCKET_AEROHEAT_PROCESS_RECORD'), ...
    'BEConfirmationPerformed',false,'physicalSolution',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
target=fullfile(folder,'spatial_root_evidence.mat');save(target,'c','solution','options','sourceBindings','report','-v7');
fid=fopen(target,'rb');guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);report.rootEvidenceSHA256=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
fid=fopen(fullfile(folder,'spatial_root_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
fprintf('%s\n',jsonencode(report));
end
