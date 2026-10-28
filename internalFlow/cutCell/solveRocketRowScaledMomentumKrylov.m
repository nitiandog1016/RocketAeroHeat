function [value,report,rowScale,originalImage]=solveRocketRowScaledMomentumKrylov(matrix,factor,rhs,options)
%SOLVEROCKETROWSCALEDMOMENTUMKRYLOV 对等价R*A求解并在原单位A中单独验收。
% 模块一：右预条件必须为P*R^{-1}，保持法向预条件器与行变换的回映一致。
% 模块二：原单位残差与尺度残差分别报告；尺度目标通过不等于原A目标通过。
% 模块三：只复用已存原A/固定因子，不调用原NS或改变任何物理场。
if nargin<4,options=struct;end
defaults=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2,'saveBasisEvidence',true,'failureEvidenceFile','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'RowMomentum:Options','未知等价行尺度选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumIterations,{'numeric'},{'scalar','real','finite','integer','positive'});validateattributes(options.relativeTolerance,{'numeric'},{'scalar','real','finite','positive','<',1});
assert(matrix.cellCount==factor.cellCount&&isequal(size(rhs),[matrix.cellCount 3])&&isreal(rhs)&&all(isfinite(rhs),'all'),'RowMomentum:Input','原A、因子或右端不一致。');
[rowScale,scaleReport]=rocketMomentumRowScales(matrix);matrix.owners=uint32(repelem((1:matrix.cellCount)',diff(double(matrix.rowPointer))));balancedRhs=rowScale.*rhs;
configuration=struct('restartLength',options.maximumIterations,'maximumCycles',1,'retainedDirections',0,'maximumArnoldiSteps',options.maximumIterations, ...
    'absoluteTolerance',max(realmin,options.relativeTolerance*norm(balancedRhs,'fro')),'maximumBasisBytes',options.maximumBasisBytes,'saveBasisEvidence',options.saveBasisEvidence,'robustArnoldi',true);
if ~isempty(options.failureEvidenceFile),configuration.failureEvidenceFile=options.failureEvidenceFile;end
timer=tic;[x,krylov]=solveRocketRecycledGMRES(@action,@precondition,balancedRhs(:),configuration);value=reshape(x,size(rhs));
originalImage=applyRocketMomentumGraph(matrix,value);remaining=rhs-originalImage;originalResidual=norm(remaining,'fro')/max(realmin,norm(rhs,'fro'));scaledResidual=norm(rowScale.*remaining,'fro')/max(realmin,norm(balancedRhs,'fro'));
assert(abs(scaledResidual-krylov.relativeResidual)<1e-8,'RowMomentum:ResidualScope','尺度剩余与保存的当前作用不一致。');
report=struct('method','one-cycle-equivalent-row-scaled-original-momentum-A-with-P-R-inverse-right-preconditioning', ...
    'options',options,'krylov',krylov,'relativeAResidual',originalResidual,'relativeRowScaledResidual',scaledResidual, ...
    'originalUnitTargetMet',originalResidual<=options.relativeTolerance,'rowScaledTargetMet',krylov.flag==0, ...
    'rowScalingReport',scaleReport,'AApplications',krylov.modelApplications+1,'solverAApplications',krylov.modelApplications,'originalUnitDiagnosticAApplications',1, ...
    'lineInverseApplications',krylov.preconditionerApplications,'elapsedSeconds',toc(timer),'baseLineInverseFixedLinear',true, ...
    'innerSolveFixedLinearMapClaimed',false,'requiresFlexibleOuterMethod',true,'newSelectedJApplications',0,'newActualNSEvaluations',0, ...
    'originalNSAcceptanceClaimed',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
    function out=action(v),y=rowScale.*applyRocketMomentumGraph(matrix,reshape(v,size(rhs)));out=y(:);end
    function out=precondition(v),y=solveRocketMomentumLines(factor,reshape(v,size(rhs))./rowScale);out=y(:);end
end
