function [value,report]=solveRocketSpatialReductionKrylov(packed,reduction,coarseFactor,rhs,options)
%SOLVEROCKETSPATIALREDUCTIONKRYLOV 用固定空间块回代预条件完整原C，保留全部方程和护栏。
% 模块一：32步等预算对照，完整原C含未在M_F中精确求逆的跨组F-F项。
% 模块二：原C通过后仍须选中J和原NS总方向验收，不把代数结果当物理解。
if nargin<5,options=struct;end
defaults=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2,'saveBasisEvidence',false,'failureEvidenceFile','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'SpatialReduction:Options','未知空间内部Krylov选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumIterations,{'numeric'},{'scalar','integer','positive'});validateattributes(options.relativeTolerance,{'numeric'},{'scalar','real','finite','positive','<',1});
assert(packed.cellCount==reduction.cellCount&&isequal(size(rhs),[packed.cellCount 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'SpatialReduction:Input','原C和五量右端不符。');
configuration=struct('restartLength',options.maximumIterations,'maximumCycles',1,'retainedDirections',0,'maximumArnoldiSteps',options.maximumIterations, ...
    'absoluteTolerance',max(realmin,options.relativeTolerance*norm(rhs,'fro')),'maximumBasisBytes',options.maximumBasisBytes,'saveBasisEvidence',options.saveBasisEvidence,'robustArnoldi',true);
if ~isempty(options.failureEvidenceFile),configuration.failureEvidenceFile=options.failureEvidenceFile;end
fApps=0;cApps=0;rectApps=0;timer=tic;[x,krylov]=solveRocketRecycledGMRES(@action,@precondition,rhs(:),configuration);value=reshape(x,size(rhs));
report=struct('method','original-full-C-GMRES-with-fixed-five-field-spatial-local-F-Schur-block-reduction','options',options,'krylov',krylov, ...
    'relativeCResidual',krylov.relativeResidual,'CApplications',krylov.modelApplications,'spatialReductionApplications',krylov.preconditionerApplications, ...
    'FLocalApplications',fApps,'coarseLocalApplications',cApps,'rectangularCouplingApplications',rectApps,'elapsedSeconds',toc(timer), ...
    'baseSpatialMapFixedLinear',true,'innerSolveFixedLinearMapClaimed',false,'originalNSAcceptanceClaimed',false,'physicalSolution',false, ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
    function out=action(v),a=applyRocketPackedGraph(packed,reshape(v,size(rhs)));out=a(:);end
    function out=precondition(v)
        [a,r]=applyRocketSpatialReduction(reduction,coarseFactor,reshape(v,size(rhs)));fApps=fApps+r.FLocalApplications;
        cApps=cApps+r.coarseLocalApplications;rectApps=rectApps+r.rectangularCouplingApplications;out=a(:);
    end
end
