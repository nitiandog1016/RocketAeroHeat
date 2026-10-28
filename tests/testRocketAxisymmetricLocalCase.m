function report=testRocketAxisymmetricLocalCase(folder)
%TESTROCKETAXISYMMETRICLOCALCASE 核对新链初猜、数值选项与物理边界分离。
% 模块一：冻结Q23初猜逐值接入，保留真实来源SHA和此前23步，当前链仍从0开始。
% 模块二：均匀选择仍等于原来流；只改求解配置不改几何、气体或给定边界。
% 模块三：不支持的初猜/非法CFL/相对扰动明确拒绝，不静默换初猜或放宽门槛。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v74_validation','local_case');end
assert(~exist(folder,'dir'),'LocalCaseTest:Overwrite','不能覆盖证据。');
c=createRocketAxisymmetricLocalCase;again=createRocketAxisymmetricLocalCase;assert(isequaln(c,again));
base=createRocketAxisymmetricFluxCase;uniform=createRocketAxisymmetricLocalCase(struct('initialGuess','uniform'));
assert(isequaln(c.mesh,base.mesh)&&isequaln(c.gas,base.gas)&&isequaln(c.boundary,base.boundary)&& ...
    isequal(c.residualScale,base.residualScale)&&isequal(uniform.Q,base.Q)&&~isequal(c.Q,base.Q));
seed=c.initialization.seed;assert(strcmp(seed.kind,'v73-uniform23')&&seed.priorAcceptedSteps==23&& ...
    ~seed.initialGuessIsBoundaryCondition&&~seed.externalCFDImported&&~seed.physicalSolution&& ...
    ~c.initialization.allFiveComponentsUniform&&c.initialization.timeStepsTaken==0&& ...
    ~c.initialization.initialStateIsSolution&&~c.inputContract.newChainIncludesPriorAcceptedSteps);
assert(uniform.initialization.allFiveComponentsUniform&&uniform.initialization.seed.priorAcceptedSteps==0&& ...
    uniform.solverConfiguration.initialPseudoCFL==.5&&c.solverConfiguration.initialPseudoCFL==1e8&& ...
    c.solverConfiguration.differenceStep==2e-6&&c.solverConfiguration.steadyResidualTolerance==2e-6&& ...
    c.solverConfiguration.conservationTolerance==5e-10);
changed=createRocketAxisymmetricLocalCase(struct('solverDifferenceStep',2e-7,'solverInitialPseudoCFL',1e7));
assert(isequal(changed.Q,c.Q)&&isequaln(changed.boundary,c.boundary)&&changed.solverConfiguration.differenceStep==2e-7);
refusals=0;reject(@()createRocketAxisymmetricLocalCase(struct('initialGuess','external-CFD')));
reject(@()createRocketAxisymmetricLocalCase(struct('solverDifferenceStep',0)));
reject(@()createRocketAxisymmetricLocalCase(struct('solverInitialPseudoCFL',1e9)));
report=struct('status','passed-explicit-local-solver-case-internal-initial-guess-contract', ...
    'matlabEngineProcessID',feature('getpid'),'sourceEvidenceSHA256',seed.sourceEvidenceSHA256, ...
    'priorAcceptedSteps',23,'newChainAcceptedStepsInitially',0,'physicalBoundaryUnchanged',true, ...
    'uniformOptionPreserved',true,'refusals',refusals,'physicalSolution',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'local_case_evidence.mat'),'report','c','uniform','changed','-v7');
fid=fopen(fullfile(folder,'local_case_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;fprintf('%s\n',jsonencode(report));
    function reject(action)
        try,action();catch,refusals=refusals+1;return;end;error('LocalCaseTest:Refusal','非法输入未拒绝。');
    end
end
