function report=studyRocketVolumeCoarse(folder)
%STUDYROCKETVOLUMECOARSE 在真实冻结端点检查完整空间Galerkin粗层及两层校正成本。
% 模块一：原Q/rate只读，释放旧图后建有限维粗矩阵，不进行原NS或非线性更新。
% 模块二：用随机粗方向独立核对RAP，粗层投影消除与细层残差变化分别报告。
% 模块三：不将粗层代数小残差当真实流场收敛或新增求解面积。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'RocketCoarseStudy:Overwrite','不得覆盖真实粗层研究。');
sourceFile=fullfile(root,'outputs','v83_validation','krylov08','wall_velocity_krylov_evidence.mat');
s=load(sourceFile,'result','gas','boundary');Q=s.result.Q;rate=s.result.rate;scale=s.result.residualScale;gas=s.gas;boundary=s.boundary;
baseInfo=struct;for key={'primitive','gradient','viscousFaceFlux'},baseInfo.(key{1})=s.result.info.(key{1});end;clear s;
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');cut=d.domain;clear d;
options=struct('cartesianBins',[2 2 2],'layerBreaks',4,'maximumCoarseUnknowns',200);space=buildRocketVolumeCoarseSpace(cut,options);
before=memory;t=tic;cache=buildCutCellSelectedDirectionCache(Q,baseInfo,cut,gas,boundary,scale);directionCacheBuildSeconds=toc(t);
coarse=buildRocketGalerkinCoarseOperator(cache,space);afterCoarse=memory;
ng=space.groupCount;coarseVector=reshape(sin((1:5*ng)*.37),ng,5);fineVector=space.prolongation*coarseVector;
modelCoarse=space.restriction*(-applyCutCellSelectedDirection(cache,fineVector));matrixCoarse=reshape(coarse.matrix*coarseVector(:),ng,5);
projectionError=norm(modelCoarse-matrixCoarse,'fro')/max(1,norm(matrixCoarse,'fro'));
assert(projectionError<1e-10,'RocketCoarseStudy:Projection','真实粗矩阵没有复现完整空间RAP作用。');
graph=buildRocketShellGraphPreconditioner(Q,cut,gas,boundary,scale,1e-7,1);rhs=rate./scale;
t=tic;[initial,graphReport]=solveRocketShellGraphPreconditioner(graph,rhs);initialSeconds=toc(t);
beforeCorrection=rhs+applyCutCellSelectedDirection(cache,initial);
t=tic;[update,linear]=solveRocketTwoLevelPreconditioner(graph,coarse,cache,rhs);applySeconds=toc(t);
afterCorrection=rhs+applyCutCellSelectedDirection(cache,update);restricted=space.restriction*afterCorrection;
coarseRemaining=norm(restricted,'fro')/max(1,norm(space.restriction*rhs,'fro'));
assert(coarseRemaining<1e-8,'RocketCoarseStudy:Correction','真实完整模型的校正后粗残差未达粗层代数门槛。');
afterApply=memory;
report=struct('status','measured-real-volume-Galerkin-projection-and-fixed-two-level-model-correction', ...
    'fineCells',size(Q,1),'coarseGroups',ng,'coarseUnknowns',5*ng,'options',options, ...
    'projectionRelativeDifference',projectionError,'coarseReciprocalConditionEstimate',coarse.reciprocalConditionEstimate, ...
    'correctedModelCoarseRelativeResidual',coarseRemaining,'initialFineModelRelativeResidual',norm(beforeCorrection,'fro')/norm(rhs,'fro'), ...
    'correctedFineModelRelativeResidual',norm(afterCorrection,'fro')/norm(rhs,'fro'), ...
    'directionCacheBuildSeconds',directionCacheBuildSeconds,'coarseBuildSeconds',coarse.buildSeconds,'initialGraphApplySeconds',initialSeconds,'twoLevelApplySeconds',applySeconds, ...
    'constructionModelApplications',coarse.constructionModelApplications,'otherValidationModelApplications',4, ...
    'matlabBytesBeforeCache',before.MemUsedMATLAB,'matlabBytesAfterCoarse',afterCoarse.MemUsedMATLAB,'matlabBytesAfterApply',afterApply.MemUsedMATLAB, ...
    'peakMemoryMeasured',false,'fineJacobianMatrixAllocated',false,'allFineCoarseColumnsStored',false, ...
    'actualOriginalNSEvaluations',0,'newAcceptedIterations',0,'physicalMeshChanged',false,'physicalSolution',false, ...
    'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
mkdir(folder);save(fullfile(folder,'volume_coarse_real_evidence.mat'),'report','sourceFile','space','coarse','coarseVector','modelCoarse','matrixCoarse', ...
    'rhs','initial','update','graphReport','linear','beforeCorrection','afterCorrection','restricted','-v7');
fid=fopen(fullfile(folder,'volume_coarse_real_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
