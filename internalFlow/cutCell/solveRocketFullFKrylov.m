function [value,report]=solveRocketFullFKrylov(graph,factor,rhs,options)
%SOLVEROCKETFULLFKRYLOV 用已有固定局部F逆迭代完整A_FF，保留全部跨组F连接。
% 模块一：子系统残差按自身右端归一化，与原C和近似S分别记录。
% 模块二：只作一个有界GMRES周期，原Arnoldi1e-8护栏保持。
% 模块三：内部解随右端变化，不把它直接替换为旧近似S的固定M_F。
if nargin<4,options=struct;end
defaults=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2,'saveBasisEvidence',false,'failureEvidenceFile','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'FullFKrylov:Options','未知完整F内部选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumIterations,{'numeric'},{'scalar','integer','positive'});validateattributes(options.relativeTolerance,{'numeric'},{'scalar','real','finite','positive','<',1});
assert(isfield(graph,'systemScope')&&strcmp(graph.systemScope,'exact-spatial-F-submatrix-A_FF-not-original-full-C')&& ...
    graph.cellCount==factor.cellCount&&isequal(size(rhs),[graph.cellCount 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'FullFKrylov:Input','完整F图、局部F因子或右端不符。');
configuration=struct('restartLength',options.maximumIterations,'maximumCycles',1,'retainedDirections',0,'maximumArnoldiSteps',options.maximumIterations, ...
    'absoluteTolerance',max(realmin,options.relativeTolerance*norm(rhs,'fro')),'maximumBasisBytes',options.maximumBasisBytes,'saveBasisEvidence',options.saveBasisEvidence,'robustArnoldi',true);
if ~isempty(options.failureEvidenceFile),configuration.failureEvidenceFile=options.failureEvidenceFile;end
timer=tic;[x,krylov]=solveRocketRecycledGMRES(@action,@precondition,rhs(:),configuration);value=reshape(x,size(rhs));
report=struct('method','complete-spatial-A_FF-GMRES-with-existing-fixed-within-patch-F-map','options',options,'krylov',krylov, ...
    'relativeFullFResidual',krylov.relativeResidual,'FApplications',krylov.modelApplications,'localFApplications',krylov.preconditionerApplications, ...
    'elapsedSeconds',toc(timer),'fullFOperatorIncludesCrossPatchBlocks',true,'baseLocalFMapFixedLinear',true, ...
    'innerSolveFixedLinearMapClaimed',false,'oldApproximateSchurIdentityReplaced',false,'originalCOrNSConvergenceClaimed',false, ...
    'physicalSolution',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
    function out=action(v),a=applyRocketPackedGraph(graph,reshape(v,size(rhs)));out=a(:);end
    function out=precondition(v),a=solveRocketSpatialFBlocks(factor,reshape(v,size(rhs)));out=a(:);end
end
