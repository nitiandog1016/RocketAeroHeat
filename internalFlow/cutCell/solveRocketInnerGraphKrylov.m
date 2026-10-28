function [update,report]=solveRocketInnerGraphKrylov(operator,rhs,options)
%SOLVEROCKETINNERGRAPHKRYLOV 有界求解当前近似系统B*x=rhs，分开评价图求逆误差。
% 模块一：B包含既有五带WLS同线块；固定一次着色扫描仅作为内部右预条件。
% 模块二：复用已验证两遍正交内核，只用一周期、零保留，不调用完整J或原NS。
% 模块三：内层GMRES映射随右端/停止条件变化，不能冒充固定线性外预条件器。
% 若以后用于外层，必须选用保存物理Z方向的柔性方法并重新通过原NS总方向验收。
if nargin<3,options=struct;end
defaults=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2,'saveBasisEvidence',false);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'InnerGraph:Options','近似系统内迭代选项无效。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumIterations,{'numeric'},{'real','finite','scalar','integer','positive'});
validateattributes(options.relativeTolerance,{'numeric'},{'real','finite','scalar','positive','<',1});
validateattributes(options.maximumBasisBytes,{'numeric'},{'real','finite','scalar','positive'});
assert(isfield(operator,'wlsNormalFactor')&&operator.fixedLinearPreconditioner&&operator.fixedSymmetricSweeps==1, ...
    'InnerGraph:Operator','本研究要求已验的五带算子与固定一次对称着色扫描。');
assert(isequal(size(rhs),[operator.sourceFaces*operator.layers 5])&&all(isfinite(rhs),'all'),'InnerGraph:RHS','五量右端非法。');
configuration=struct('restartLength',options.maximumIterations,'maximumCycles',1,'retainedDirections',0, ...
    'maximumArnoldiSteps',options.maximumIterations,'absoluteTolerance',max(realmin,options.relativeTolerance*norm(rhs,'fro')), ...
    'maximumBasisBytes',options.maximumBasisBytes,'saveBasisEvidence',options.saveBasisEvidence);
lastScan=struct;timer=tic;[vector,krylov]=solveRocketRecycledGMRES(@graphAction,@linePrecondition,rhs(:),configuration);
update=reshape(vector,size(rhs));
report=struct('method','bounded-one-cycle-right-GMRES-for-approximate-five-band-graph-B', ...
    'linearSystem','B*x=rhs','options',options,'krylov',krylov,'lastLineScanReport',lastScan, ...
    'approximateBApplications',krylov.modelApplications,'linePreconditionerApplications',krylov.preconditionerApplications, ...
    'lineScanDiagnosticBApplications',krylov.preconditionerApplications, ...
    'totalBApplicationsIncludingScanDiagnostics',krylov.modelApplications+krylov.preconditionerApplications, ...
    'relativeBResidual',krylov.relativeResidual,'elapsedSeconds',toc(timer), ...
    'fullSelectedJApplications',0,'newActualNSEvaluations',0,'outerFixedLinearMapClaimed',false, ...
    'requiresFlexibleOuterMethod',true,'steadyNewtonEquationSolved',false,'originalNSAcceptanceClaimed',false, ...
    'physicalEquationChanged',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
    function out=graphAction(vector)
        value=applyRocketWlsNormalOperator(operator,reshape(vector,size(rhs)));out=value(:);
    end
    function out=linePrecondition(vector)
        [value,lastScan]=solveRocketWlsNormalPreconditioner(operator,reshape(vector,size(rhs)));out=value(:);
    end
end
