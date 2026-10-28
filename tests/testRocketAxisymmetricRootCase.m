function report=testRocketAxisymmetricRootCase(folder)
%TESTROCKETAXISYMMETRICROOTCASE 核对内部定常解胶囊、物理输入和新BE链身份。
% 模块一：普通/双double两胶囊均为真实根，不能作为给定入口或物理精度认证。
% 模块二：三次Newton准备与本链零BE步分开，根源码/输入不匹配必须拒绝。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v75_validation','root_case');end
assert(~exist(folder,'dir'),'RootCaseTest:Overwrite','不能覆盖测试。');
ordinary=createRocketAxisymmetricRootCase;accurate=createRocketAxisymmetricRootCase(struct('rootArithmetic','two-double-euler'));
assert(isequaln(ordinary,createRocketAxisymmetricRootCase)&&isequaln(ordinary.mesh,accurate.mesh)&&isequaln(ordinary.gas,accurate.gas));
assert(~isfield(ordinary.boundary,'interiorFluxArithmetic')&&strcmp(accurate.boundary.interiorFluxArithmetic,'two-double-euler'));
for c={ordinary,accurate}
    item=c{1};seed=item.initialization.seed;
    assert(strcmp(seed.kind,'v75-spatial-root')&&seed.priorNewtonIterations==3&&seed.priorAcceptedSteps==0&& ...
        seed.parentAcceptedSteps==6&&seed.parentSeedPriorAcceptedSteps==23&&~seed.initialGuessIsBoundaryCondition&& ...
        ~seed.externalCFDImported&&~seed.physicalSolution&&seed.spatialRootResidualNorm<=2e-6&& ...
        item.solverConfiguration.initialPseudoCFL==5&&item.solverConfiguration.maximumPseudoCFL==1e6&& ...
        ~item.initialization.allFiveComponentsUniform&&item.initialization.timeStepsTaken==0);
    [rate,info]=axisymmetricRANSResidual(item.Q,item.mesh,item.gas,item.boundary);
    assert(max(abs(rate)./item.residualScale,[],'all')==seed.spatialRootResidualNorm&&~info.physicalSolution);
end
refusals=0;reject(@()createRocketAxisymmetricRootCase(struct('rootArithmetic','unknown')),'RocketRootCase:Arithmetic');
reject(@()createRocketAxisymmetricRootCase(struct('initialGuess','uniform')),'RocketRootCase:Guess');
reject(@()createRocketAxisymmetricRootCase(struct('flight',struct('Twall',301))),'RocketRootCase:Inputs');
report=struct('status','passed-bound-internal-spatial-root-factory-and-separate-BE-identity', ...
    'matlabEngineProcessID',feature('getpid'),'sourceRootCapsules',{{ordinary.initialization.seed.sourceEvidenceSHA256,accurate.initialization.seed.sourceEvidenceSHA256}}, ...
    'NewtonPreparationIterations',3,'newBEChainInitialAcceptedSteps',0,'refusals',refusals, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'root_case_evidence.mat'),'report','ordinary','accurate','-v7');
fid=fopen(fullfile(folder,'root_case_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));clear guard;
fprintf('%s\n',jsonencode(report));
    function reject(action,id)
        try,action();catch exception,assert(strcmp(exception.identifier,id));refusals=refusals+1;return;end;error('RootCaseTest:Refusal','非法输入未拒绝。');
    end
end
