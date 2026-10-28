function report=testRocketVolumeCoarse(folder)
%TESTROCKETVOLUMECOARSE 用原NS粗方向差分核对R*(-J)*P与两层校正。
% 模块一：真实体积R/P恒等式、连通分组和粗层直接LU，不使用人工替代残差。
% 模块二：逐粗未知量中心差分原完整NS，独立于缓存方向构造小矩阵。
% 模块三：固定两层映射的叠加及校正后粗残差，粗层不当作物理解。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'RocketCoarseTest:Overwrite','不得覆盖粗层验证。');mkdir(folder);
s=load(fullfile(root,'outputs','v83_validation','wlsunit01','wls_direction_unit_evidence.mat'),'Q','baseInfo','cut','gas','boundary','scale');
Q=s.Q;baseInfo=s.baseInfo;cut=s.cut;gas=s.gas;boundary=s.boundary;scale=s.scale;clear s;
cache=buildCutCellSelectedDirectionCache(Q,baseInfo,cut,gas,boundary,scale);
space=buildRocketVolumeCoarseSpace(cut,struct('cartesianBins',[2 2 2],'layerBreaks',2,'maximumCoarseUnknowns',160));
coarse=buildRocketGalerkinCoarseOperator(cache,space);ng=space.groupCount;n=size(Q,1);dimension=5*ng;
fdMatrix=zeros(dimension);calls=0;h=1e-6;
for b=1:5,for g=1:ng
    direction=zeros(n,5);direction(space.cellToGroup==g,b)=1;
    rp=cutCellNavierStokesResidual(Q+h*direction.*scale,cut,gas,boundary);
    rm=cutCellNavierStokesResidual(Q-h*direction.*scale,cut,gas,boundary);calls=calls+2;
    column=-space.restriction*((rp-rm)./scale/(2*h));fdMatrix(:,(b-1)*ng+g)=column(:);
end,end
fdError=norm(coarse.matrix-fdMatrix,'fro')/norm(fdMatrix,'fro');assert(fdError<2e-6,'RocketCoarseTest:FD','粗矩阵不符合原NS逐粗未知量差分。');
coarseRhs=reshape(sin((1:dimension)*.37),ng,5);[coarseSolution,coarseReport]=solveRocketGalerkinCoarseSystem(coarse,coarseRhs);
reference=reshape(coarse.matrix\coarseRhs(:),ng,5);solveError=norm(coarseSolution-reference,'fro')/max(1,norm(reference,'fro'));
assert(solveError<1e-11);
rhs=reshape(cos((1:5*n)*.23),n,5);graph=buildRocketShellGraphPreconditioner(Q,cut,gas,boundary,scale,1e-7,1);
[update,linear]=solveRocketTwoLevelPreconditioner(graph,coarse,cache,rhs);
remaining=space.restriction*(rhs+applyCutCellSelectedDirection(cache,update));
coarseRemaining=norm(remaining,'fro')/max(1,norm(space.restriction*rhs,'fro'));assert(coarseRemaining<1e-8,'RocketCoarseTest:Correction','校正后真实完整模型的粗残差未被消除。');
second=reshape(sin((1:5*n)*.47),n,5);other=solveRocketTwoLevelPreconditioner(graph,coarse,cache,second);
combined=solveRocketTwoLevelPreconditioner(graph,coarse,cache,.7*rhs-1.3*second);expected=.7*update-1.3*other;
linearity=norm(combined-expected,'fro')/max(1,norm(expected,'fro'));assert(linearity<1e-11,'RocketCoarseTest:Linearity','固定两层预条件叠加检查失败。');
rejected=0;badOptions=space.options;badOptions.maximumCoarseUnknowns=5;
try,buildRocketVolumeCoarseSpace(cut,badOptions);catch ex,assert(strcmp(ex.identifier,'RocketCoarse:Budget'));rejected=rejected+1;end
bad=coarse;bad.baseQ(1,1)=1.01*bad.baseQ(1,1);
try,solveRocketTwoLevelPreconditioner(graph,bad,cache,rhs);catch ex,assert(strcmp(ex.identifier,'RocketCoarse:Cache'));rejected=rejected+1;end
assert(rejected==2);
report=struct('status','passed-actual-NS-coarse-FD-volume-transfer-LU-and-fixed-two-level-correction', ...
    'cellCount',n,'coarseGroups',ng,'coarseUnknowns',dimension,'originalNSCoarseFDRelativeDifference',fdError, ...
    'coarseLUReferenceRelativeDifference',solveError,'correctedModelCoarseRelativeResidual',coarseRemaining, ...
    'fixedTwoLevelLinearityRelativeDifference',linearity,'restrictionProlongationIdentityError',space.restrictionProlongationIdentityError, ...
    'totalVolumeRelativeDifference',space.totalVolumeRelativeDifference,'reciprocalConditionEstimate',coarse.reciprocalConditionEstimate, ...
    'actualOriginalNSEvaluations',calls,'constructionModelApplications',coarse.constructionModelApplications, ...
    'invalidBudgetOrCacheRejections',rejected,'physicalMeshChanged',false,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'volume_coarse_unit_evidence.mat'),'report','Q','baseInfo','cut','gas','boundary','scale','space','coarse','fdMatrix', ...
    'coarseRhs','coarseSolution','reference','coarseReport','rhs','graph','update','linear','remaining','second','other','combined','-v7');
fid=fopen(fullfile(folder,'volume_coarse_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
