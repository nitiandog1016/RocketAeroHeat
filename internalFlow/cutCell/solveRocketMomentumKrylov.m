function [value,report]=solveRocketMomentumKrylov(matrix,factor,rhs,options)
%SOLVEROCKETMOMENTUMKRYLOV 以固定法向列逆预条件，有界求解完整原三动量A。
% 模块一：完整A保留跨列作用；一个GMRES周期核对真实A剩余。
% 模块二：内层求解随右端改变，不能作为固定线性Schur作用，后续需柔性外层。
% 模块三：不调用原NS或更新任何物理场，原Arnoldi正交/组合门槛保持。
if nargin<4,options=struct;end
defaults=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2,'saveBasisEvidence',true,'failureEvidenceFile','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'MomentumKrylov:Options','未知三动量内部选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumIterations,{'numeric'},{'scalar','integer','positive'});validateattributes(options.relativeTolerance,{'numeric'},{'scalar','real','finite','positive','<',1});
assert(matrix.cellCount==factor.cellCount&&isequal(size(rhs),[matrix.cellCount 3])&&all(isfinite(rhs),'all'),'MomentumKrylov:Input','原A、因子或右端不一致。');
matrix.owners=uint32(repelem((1:matrix.cellCount)',diff(double(matrix.rowPointer))));
configuration=struct('restartLength',options.maximumIterations,'maximumCycles',1,'retainedDirections',0,'maximumArnoldiSteps',options.maximumIterations, ...
    'absoluteTolerance',max(realmin,options.relativeTolerance*norm(rhs,'fro')),'maximumBasisBytes',options.maximumBasisBytes,'saveBasisEvidence',options.saveBasisEvidence,'robustArnoldi',true);
if ~isempty(options.failureEvidenceFile),configuration.failureEvidenceFile=options.failureEvidenceFile;end
timer=tic;[x,krylov]=solveRocketRecycledGMRES(@action,@precondition,rhs(:),configuration);value=reshape(x,size(rhs));
report=struct('method','one-cycle-full-original-momentum-A-right-GMRES-with-fixed-normal-line-inverses','options',options,'krylov',krylov, ...
    'relativeAResidual',krylov.relativeResidual,'AApplications',krylov.modelApplications,'lineInverseApplications',krylov.preconditionerApplications,'elapsedSeconds',toc(timer), ...
    'baseLineInverseFixedLinear',true,'innerSolveFixedLinearMapClaimed',false,'requiresFlexibleOuterMethod',true, ...
    'newSelectedJApplications',0,'newActualNSEvaluations',0,'originalNSAcceptanceClaimed',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
    function out=action(v),y=applyRocketMomentumGraph(matrix,reshape(v,size(rhs)));out=y(:);end
    function out=precondition(v),y=solveRocketMomentumLines(factor,reshape(v,size(rhs)));out=y(:);end
end
