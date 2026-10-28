function report=solveRocketContactRefinedStudy
%SOLVEROCKETCONTACTREFINEDSTUDY 以已收敛同方程粗根作初猜的独立144层受控试验。
% 保留原失败试验，原24次预算、2e-6/5e-10门槛及温域/SA/星态护栏保持。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
kind='nested144_from_hllc72';folder=fullfile(root,'outputs','v78_validation',kind);
assert(~exist(folder,'dir'),'RocketContactRefined:Overwrite','不能覆盖真实初猜对照。');
c=createRocketAxisymmetricContactRefinedCase;sourceBindings=rocketContactRefinedSourceBindings;
options=struct('residualScale',c.residualScale,'steadyResidualTolerance',2e-6,'conservationTolerance',5e-10, ...
    'differenceStep',2e-6,'maximumIterations',24,'maximumLineSearch',14);
mkdir(folder);save(fullfile(folder,'initial_case.mat'),'c','sourceBindings','options','-v7');
clock=tic;solution=solveAxisymmetricSteadyNewton(c.Q,c.mesh,c.gas,c.boundary,options);elapsed=toc(clock);
assert(isequaln(sourceBindings,rocketContactRefinedSourceBindings),'RocketContactRefined:Source','实际执行源变化。');
report=struct('kind',kind,'cellCount',c.mesh.cellCount,'interiorConvectiveFlux','hllc', ...
    'spatialRootConverged',solution.spatialRootConverged,'originalResidual',solution.residualNorm,'stopReason',solution.stopReason, ...
    'acceptedNewtonIterations',solution.acceptedNewtonIterations,'residualEvaluations',solution.residualEvaluations, ...
    'boundaryFluxBalanceRelative',solution.boundaryFluxBalanceRelative,'elapsedSeconds',elapsed, ...
    'engineProcessID',feature('getpid'),'processRecord',getenv('ROCKET_AEROHEAT_PROCESS_RECORD'), ...
    'BEConfirmationPerformed',false,'physicalSolution',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
target=fullfile(folder,'spatial_root_evidence.mat');save(target,'c','solution','options','sourceBindings','report','-v7');
fid=fopen(target,'rb');guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);report.rootEvidenceSHA256=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
fid=fopen(fullfile(folder,'spatial_root_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
