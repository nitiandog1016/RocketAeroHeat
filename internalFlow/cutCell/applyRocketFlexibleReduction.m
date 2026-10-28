function [value,report]=applyRocketFlexibleReduction(reduction,coarseFactor,Fgraph,columnScale,rhs,options)
%APPLYROCKETFLEXIBLEREDUCTION 用有界全F内解和旧粗近似构成显式可变右预条件。
% 模块一：先近似解F，再把E耦合传到全部C点；旧S局部逆仅作近似预条件。
% 模块二：C校正通过B回到全部F点，再次有界解F；保留每个原CV的全部五量。
% 模块三：LSQR随右端变化且可未达F目标，不能称固定线性映射、精确Schur或物理解。
if nargin<6,options=struct;end
defaults=struct('maximumFIterations',32,'FRelativeTolerance',1e-3,'evidenceFile','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'FlexibleReduction:Options','未知可变空间回代选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumFIterations,{'numeric'},{'scalar','finite','integer','positive'});
assert(isscalar(options.FRelativeTolerance)&&isfinite(options.FRelativeTolerance)&&options.FRelativeTolerance>0&&options.FRelativeTolerance<=1e-3,'FlexibleReduction:Options','原F目标不能放宽超过1e-3。');
split=reduction.split;
assert(strcmp(reduction.kind,'five-field-spatial-local-F-approximate-Schur-reduction')&& ...
    isequal(size(rhs),[reduction.cellCount 5])&&isreal(rhs)&&all(isfinite(rhs),'all')&& ...
    coarseFactor.cellCount==split.coarseCells&&Fgraph.cellCount==split.fineCells&& ...
    isequal(Fgraph.originalPhysicalCVs,split.fineOriginalCVs),'FlexibleReduction:Input','原C/F映射、完整F图、旧粗对象或五量右端不符。');
if ~isempty(options.evidenceFile),assert(exist(options.evidenceFile,'file')~=2,'FlexibleReduction:Overwrite','不得覆盖内层作用证据。');end
bF=rhs(split.fineOriginalCVs,:);bC=rhs(split.coarseOriginalCVs,:);
inner=struct('maximumIterations',options.maximumFIterations,'relativeTolerance',options.FRelativeTolerance);
[firstF,firstReport]=solveRocketFullFLeastSquares(Fgraph,reduction.Ffactor,bF,zeros(size(bF)),inner,columnScale);
coarseRHS=bC-applyRocketRectangularCoupling(reduction.E,firstF);
xC=solveRocketCoupledPatches(coarseFactor,coarseRHS);
secondRHS=bF-applyRocketRectangularCoupling(reduction.B,xC);
[secondF,secondReport]=solveRocketFullFLeastSquares(Fgraph,reduction.Ffactor,secondRHS,zeros(size(bF)),inner,columnScale);
value=zeros(size(rhs));value(split.fineOriginalCVs,:)=secondF;value(split.coarseOriginalCVs,:)=xC;
assert(all(isfinite(value),'all'),'FlexibleReduction:Nonfinite','完整可变块回代非有限。');
counts=firstReport.operatorApplications;
for key=fieldnames(counts)',counts.(key{1})=counts.(key{1})+secondReport.operatorApplications.(key{1});end
report=struct('method','variable-bounded-scaled-F-solves-old-S-local-approximation-and-full-block-back-substitution','options',options, ...
    'firstFReport',firstReport,'secondFReport',secondReport,'operatorApplications',counts,'coarseLocalApplications',1,'rectangularCouplingApplications',2, ...
    'FInnerTargetsNotMet',double(~firstReport.fullFInternalTargetMet)+double(~secondReport.fullFInternalTargetMet), ...
    'fixedLinearMapClaimed',false,'oldSUsedOnlyAsApproximatePreconditioningObject',true,'oldSchurIdentityClaimedForNewMap',false, ...
    'fullFFOrCoarseInverseClaimed',false,'physicalSolution',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
if ~isempty(options.evidenceFile)
    save(options.evidenceFile,'rhs','firstF','coarseRHS','xC','secondRHS','secondF','firstReport','secondReport','report','-v7');
end
end
