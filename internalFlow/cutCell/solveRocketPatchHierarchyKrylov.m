function [value,report]=solveRocketPatchHierarchyKrylov(hierarchy,rhs,options)
%SOLVEROCKETPATCHHIERARCHYKRYLOV 用固定五量代数V循环预条件完整原C，保留原验收门槛。
% 模块一：外层完整C与每个预条件内部C作用分别记账，不计作原NS。
% 模块二：一次有界GMRES，完整C通过后仍需原NS总方向验收。
if nargin<3,options=struct;end
defaults=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2,'saveBasisEvidence',false,'failureEvidenceFile','');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'PatchHierarchy:Options','未知内部Krylov选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.maximumIterations,{'numeric'},{'scalar','integer','positive'});
validateattributes(options.relativeTolerance,{'numeric'},{'scalar','real','finite','positive','<',1});
assert(isequal(size(rhs),[hierarchy.cellCount 5])&&all(isfinite(rhs),'all'),'PatchHierarchy:Input','完整五量右端不符。');
configuration=struct('restartLength',options.maximumIterations,'maximumCycles',1,'retainedDirections',0,'maximumArnoldiSteps',options.maximumIterations, ...
    'absoluteTolerance',max(realmin,options.relativeTolerance*norm(rhs,'fro')),'maximumBasisBytes',options.maximumBasisBytes,'saveBasisEvidence',options.saveBasisEvidence,'robustArnoldi',true);
if ~isempty(options.failureEvidenceFile),configuration.failureEvidenceFile=options.failureEvidenceFile;end
nl=numel(hierarchy.levels);insideActions=zeros(nl,1);insideLocal=insideActions;insideTerminal=insideActions;timer=tic;
[x,krylov]=solveRocketRecycledGMRES(@action,@precondition,rhs(:),configuration);value=reshape(x,size(rhs));
report=struct('method','original-full-C-GMRES-with-fixed-five-field-strong-connection-V-cycle','options',options,'krylov',krylov, ...
    'relativeCResidual',krylov.relativeResidual,'solverCApplications',krylov.modelApplications,'VcycleApplications',krylov.preconditionerApplications, ...
    'preconditionerCApplicationsByLevel',insideActions,'preconditionerLocalApplicationsByLevel',insideLocal, ...
    'preconditionerTerminalApplicationsByLevel',insideTerminal,'totalRootCApplications',krylov.modelApplications+insideActions(1), ...
    'elapsedSeconds',toc(timer),'baseHierarchyFixedLinear',true,'innerSolveFixedLinearMapClaimed',false, ...
    'originalNSAcceptanceClaimed',false,'physicalSolution',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
    function out=action(v),a=applyRocketPackedGraph(hierarchy.levels{1}.packed,reshape(v,size(rhs)));out=a(:);end
    function out=precondition(v)
        [a,r]=applyRocketPatchHierarchy(hierarchy,reshape(v,size(rhs)));insideActions=insideActions+r.CApplicationsByLevel;
        insideLocal=insideLocal+r.localMapApplicationsByLevel;insideTerminal=insideTerminal+r.terminalApplicationsByLevel;out=a(:);
    end
end
