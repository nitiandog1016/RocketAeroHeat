function [value,report]=solveRocketCoupledPatchKrylov(packed,factor,rhs,options)
%SOLVEROCKETCOUPLEDPATCHKRYLOV 用保留五量耦合的局部块求解完整原C，不作动量Schur分裂。
% 模块一：完整C作用含全部原块；局部逆仅为固定右预条件。
% 模块二：一次有界GMRES，原正交/组合1e-8门槛不变，内部结果不冒称原NS解。
if nargin<4,options=struct;end
defaults=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2,'saveBasisEvidence',false,'failureEvidenceFile','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'CoupledPatchKrylov:Options','未知五量内部选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
assert(packed.cellCount==factor.cellCount&&isequal(size(rhs),[packed.cellCount 5])&&all(isfinite(rhs),'all'),'CoupledPatchKrylov:Input','原C与右端不匹配。');
configuration=struct('restartLength',options.maximumIterations,'maximumCycles',1,'retainedDirections',0,'maximumArnoldiSteps',options.maximumIterations, ...
    'absoluteTolerance',max(realmin,options.relativeTolerance*norm(rhs,'fro')),'maximumBasisBytes',options.maximumBasisBytes,'saveBasisEvidence',options.saveBasisEvidence,'robustArnoldi',true);
if ~isempty(options.failureEvidenceFile),configuration.failureEvidenceFile=options.failureEvidenceFile;end
timer=tic;[x,krylov]=solveRocketRecycledGMRES(@action,@precondition,rhs(:),configuration);value=reshape(x,size(rhs));
report=struct('method','full-five-conservation-C-GMRES-with-bounded-strong-connection-local-inverses','options',options,'krylov',krylov, ...
    'relativeCResidual',krylov.relativeResidual,'CApplications',krylov.modelApplications,'localBlockApplications',krylov.preconditionerApplications,'elapsedSeconds',toc(timer), ...
    'baseLocalMapFixedLinear',true,'innerSolveFixedLinearMapClaimed',false,'requiresFlexibleOuterMethod',true,'originalNSAcceptanceClaimed',false,'physicalSolution',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
    function out=action(v),a=applyRocketPackedGraph(packed,reshape(v,size(rhs)));out=a(:);end
    function out=precondition(v),a=solveRocketCoupledPatches(factor,reshape(v,size(rhs)));out=a(:);end
end
