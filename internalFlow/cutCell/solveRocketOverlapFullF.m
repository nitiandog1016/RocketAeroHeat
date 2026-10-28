function [value,report]=solveRocketOverlapFullF(graph,factor,rhs,options)
%SOLVEROCKETOVERLAPFULLF 使用受限核心行重叠逆求完整A_FF，全部原跨组系数仍参与。
% 模块一：只作一次有界GMRES，原正交/组合1e-8门槛保持。
% 模块二：F自身右端范数与父原C另列，不改变旧近似S、原NS或物理场。
if nargin<4,options=struct;end
defaults=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2,'saveBasisEvidence',false,'failureEvidenceFile','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'OverlapF:Options','未知重叠F内部选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumIterations,{'numeric'},{'scalar','integer','positive'});validateattributes(options.relativeTolerance,{'numeric'},{'scalar','real','finite','positive','<',1});
assert(isfield(graph,'systemScope')&&strcmp(graph.systemScope,'exact-spatial-F-submatrix-A_FF-not-original-full-C')&& ...
    graph.cellCount==factor.cellCount&&isequal(size(rhs),[graph.cellCount 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'OverlapF:Input','完整F图、重叠因子或右端不符。');
configuration=struct('restartLength',options.maximumIterations,'maximumCycles',1,'retainedDirections',0,'maximumArnoldiSteps',options.maximumIterations, ...
    'absoluteTolerance',max(realmin,options.relativeTolerance*norm(rhs,'fro')),'maximumBasisBytes',options.maximumBasisBytes,'saveBasisEvidence',options.saveBasisEvidence,'robustArnoldi',true);
if ~isempty(options.failureEvidenceFile),configuration.failureEvidenceFile=options.failureEvidenceFile;end
timer=tic;[x,krylov]=solveRocketRecycledGMRES(@action,@precondition,rhs(:),configuration);value=reshape(x,size(rhs));
report=struct('method','complete-A_FF-GMRES-with-certified-extended-inverses-returning-core-rows-only','options',options,'krylov',krylov, ...
    'relativeFullFResidual',krylov.relativeResidual,'FApplications',krylov.modelApplications,'localFApplications',krylov.preconditionerApplications, ...
    'elapsedSeconds',toc(timer),'fullFOperatorIncludesCrossPatchBlocks',true,'baseLocalFMapFixedLinear',true, ...
    'innerSolveFixedLinearMapClaimed',false,'oldApproximateSchurIdentityReplaced',false,'originalCOrNSConvergenceClaimed',false, ...
    'restrictedCoreOutput',true,'physicalSolution',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
    function out=action(v),a=applyRocketPackedGraph(graph,reshape(v,size(rhs)));out=a(:);end
    function out=precondition(v),a=solveRocketRestrictedOverlap(factor,reshape(v,size(rhs)));out=a(:);end
end
