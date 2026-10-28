function report=probeRocketShellGraphPreconditioner(folder,sourceCase,sweeps,pseudoDt)
%PROBEROCKETSHELLGRAPHPRECONDITIONER 在真实保存端点测侧向图成本与固定线性作用。
% 模块一：不调用原NS，不更新Q；直接使用已保存的完整基态和真实方向。
% 模块二：测缓存大小/进程内存观察值、构建与扫描时间，不能称完整求解加速。
% 模块三：同方向对已独立审计的图作用；另作叠加检验，不放宽任何物理门槛。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
if nargin<4,pseudoDt=1e-4;end
validateattributes(pseudoDt,{'numeric'},{'real','finite','positive','scalar'});
assert(~exist(folder,'dir'),'GraphProbe:Overwrite','不得覆盖真实图预条件探针。');
assert(strcmp(sourceCase,'krylov06'),'GraphProbe:Source','当前成本对照明确绑定已保存的krylov06方向。');
val=fullfile(root,'outputs','v83_validation');sourceFile=fullfile(val,sourceCase,'wall_velocity_krylov_evidence.mat');
s=load(sourceFile,'result','gas','boundary');Q=s.result.initialQ;scale=s.result.residualScale;gas=s.gas;boundary=s.boundary;
rhs=s.result.initialRate./scale;directionU=s.result.linearLog{1}.steadyDirectionEvidence.directionU;clear s;
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');cut=d.domain;clear d;
before=memory;t=tic;operator=buildRocketShellGraphPreconditioner(Q,cut,gas,boundary,scale,pseudoDt,sweeps);buildSeconds=toc(t);afterBuild=memory;
buffers=whos('operator');operatorBytes=buffers.bytes;t=tic;[update,scan]=solveRocketShellGraphPreconditioner(operator,rhs);applySeconds=toc(t);afterApply=memory;
graphAction=directionU/pseudoDt-applyRocketShellGraphOperator(operator,directionU);
prior=load(fullfile(val,'lateralrocket01','lateral_coupling_evidence.mat'),'exactSideGraphAction');
priorActionError=norm(graphAction(:)-prior.exactSideGraphAction(:))/norm(prior.exactSideGraphAction(:));
assert(priorActionError<2e-11,'GraphProbe:Action','实际图装配没有复现此前独立核对的全侧面作用。');
second=rhs.*reshape(cos((1:numel(rhs))*.11),size(rhs));other=solveRocketShellGraphPreconditioner(operator,second);
combined=solveRocketShellGraphPreconditioner(operator,.7*rhs-1.3*second);expected=.7*update-1.3*other;
linearity=norm(combined(:)-expected(:))/max(1,norm(expected(:)));assert(linearity<1e-11,'GraphProbe:Linearity','真实图扫描叠加关系未通过。');
color=operator.color;
report=struct('status','measured-real-graph-cache-fixed-scan-cost-and-prior-independent-action', ...
    'sourceCase',sourceCase,'cellCount',size(Q,1),'sideFaceCount',numel(operator.sideMinus),'colorCount',numel(operator.colorGroups), ...
    'fixedSymmetricSweeps',sweeps,'preconditionerPseudoDt',pseudoDt,'operatorCacheBytes',operatorBytes,'buildSeconds',buildSeconds,'applySeconds',applySeconds, ...
    'matlabBytesBeforeBuild',before.MemUsedMATLAB,'matlabBytesAfterBuild',afterBuild.MemUsedMATLAB,'matlabBytesAfterApply',afterApply.MemUsedMATLAB, ...
    'peakMemoryMeasured',false,'fullSparseGraphAllocated',false,'graphPreconditionerRelativeResidual',scan.graphRelativeResidual, ...
    'maximumCachedLineFactorDefect',operator.maximumFactorDefect,'priorIndependentActionRelativeDifference',priorActionError, ...
    'linearSuperpositionRelativeDifference',linearity,'newOriginalNSEvaluations',0,'newAcceptedIterations',0, ...
    'fullOriginalNSJacobian',false,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
mkdir(folder);save(fullfile(folder,'graph_probe_evidence.mat'),'report','sourceFile','scale','rhs','update','scan','graphAction','second','other','combined','color','-v7');
fid=fopen(fullfile(folder,'graph_probe_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
