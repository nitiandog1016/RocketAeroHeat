function [value,report]=solveRocketFlexibleReductionKrylov(packed,reduction,coarseFactor,Fgraph,columnScale,rhs,options)
%SOLVEROCKETFLEXIBLEREDUCTIONKRYLOV 固定完整父C作用，显式可变右预条件，单周期无回收。
% 模块一：内层LSQR的工作量/右端可变，不把它当固定Schur算子。
% 模块二：复用保存每次实际Z及A*Z的稳健Arnoldi内核，最终用完整父C核对组合。
% 模块三：原完整C的1e-3与Arnoldi正交/组合1e-8门槛保持；单列所有内层实际调用。
if nargin<7,options=struct;end
defaults=struct('maximumIterations',8,'relativeTolerance',1e-3,'innerIterationSequence',32,'innerRelativeTolerance',1e-3, ...
    'maximumInnerOperatorApplications',2048,'maximumBasisBytes',768*1024^2,'saveBasisEvidence',false,'preconditionerEvidenceFolder','','failureEvidenceFile','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'FlexibleReduction:Options','未知完整父C可变预条件选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumIterations,{'numeric'},{'scalar','finite','integer','positive'});
sequence=options.innerIterationSequence(:);assert(~isempty(sequence)&&isreal(sequence)&&all(isfinite(sequence)&sequence>=1&sequence==floor(sequence)),'FlexibleReduction:Options','内层工作量序列必须为正整数。');
assert(isscalar(options.relativeTolerance)&&options.relativeTolerance>0&&options.relativeTolerance<=1e-3&& ...
    isscalar(options.innerRelativeTolerance)&&options.innerRelativeTolerance>0&&options.innerRelativeTolerance<=1e-3,'FlexibleReduction:Options','原父C/F目标不能放宽。');
assert(packed.cellCount==reduction.cellCount&&isequal(size(rhs),[packed.cellCount 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'FlexibleReduction:Input','完整父C或五量右端非法。');
upper=options.maximumIterations*2*(2*max(sequence)+16);
assert(isscalar(options.maximumInnerOperatorApplications)&&isfinite(options.maximumInnerOperatorApplications)&&upper<=options.maximumInnerOperatorApplications,'FlexibleReduction:WorkBudget','可变内层显式最坏调用预算不足。');
if ~isempty(options.preconditionerEvidenceFolder),assert(~exist(options.preconditionerEvidenceFolder,'dir'),'FlexibleReduction:Overwrite','不得覆盖预条件证据目录。');mkdir(options.preconditionerEvidenceFolder);end
configuration=struct('restartLength',options.maximumIterations,'maximumCycles',1,'retainedDirections',0,'maximumArnoldiSteps',options.maximumIterations, ...
    'absoluteTolerance',max(realmin,options.relativeTolerance*norm(rhs,'fro')),'maximumBasisBytes',options.maximumBasisBytes,'saveBasisEvidence',options.saveBasisEvidence,'robustArnoldi',true);
if ~isempty(options.failureEvidenceFile),configuration.failureEvidenceFile=options.failureEvidenceFile;end
application=0;innerWork=0;preconditionerReports={};counts=struct('F',0,'FT',0,'M',0,'MT',0);coarseApps=0;rectApps=0;unmet=0;timer=tic;
[x,krylov]=solveRocketRecycledGMRES(@action,@precondition,rhs(:),configuration);value=reshape(x,size(rhs));
report=struct('method','single-cycle-flexible-right-preconditioned-complete-original-C','options',options,'krylov',krylov, ...
    'relativeCResidual',krylov.relativeResidual,'internalCTargetMet',krylov.relativeResidual<=options.relativeTolerance, ...
    'CApplications',krylov.modelApplications,'flexiblePreconditionerApplications',application,'innerOperatorApplications',counts, ...
    'coarseLocalApplications',coarseApps,'rectangularCouplingApplications',rectApps,'FInnerTargetsNotMet',unmet, ...
    'preconditionerReports',{preconditionerReports},'innerFAndFTWorkUpperBound',upper,'actualInnerFAndFTApplications',innerWork, ...
    'elapsedSeconds',toc(timer),'completeParentOperatorFixed',true,'preconditionerFixedLinearMapClaimed',false, ...
    'actualStoredZUsedForCombination',true,'oldSchurIdentityClaimedForNewMap',false, ...
    'originalNSAcceptanceClaimed',false,'physicalSolution',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
    function out=action(v),a=applyRocketPackedGraph(packed,reshape(v,size(rhs)));out=a(:);end
    function out=precondition(v)
        application=application+1;budget=sequence(1+mod(application-1,numel(sequence)));
        assert(innerWork+2*(2*budget+16)<=options.maximumInnerOperatorApplications,'FlexibleReduction:WorkBudget','剩余内层最坏调用预算不足，拒绝启动。');
        p=struct('maximumFIterations',budget,'FRelativeTolerance',options.innerRelativeTolerance,'evidenceFile','');
        if ~isempty(options.preconditionerEvidenceFolder),p.evidenceFile=fullfile(options.preconditionerEvidenceFolder,sprintf('application_%03d.mat',application));end
        [a,r]=applyRocketFlexibleReduction(reduction,coarseFactor,Fgraph,columnScale,reshape(v,size(rhs)),p);
        for key=fieldnames(counts)',counts.(key{1})=counts.(key{1})+r.operatorApplications.(key{1});end
        innerWork=innerWork+r.operatorApplications.F+r.operatorApplications.FT;
        assert(innerWork<=options.maximumInnerOperatorApplications,'FlexibleReduction:WorkBudget','实际内层调用超预算，不能称可用方向。');
        coarseApps=coarseApps+1;rectApps=rectApps+2;unmet=unmet+r.FInnerTargetsNotMet;
        preconditionerReports{end+1}=r;out=a(:);
    end
end
