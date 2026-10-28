function report=studyRocketLeastSquaresCoarse(folder,kind)
%STUDYROCKETLEASTSQUARESCOARSE 验证细层最小二乘粗校正，不重算物理解或原NS。
% 模块一：复用已测连通体积聚合，稀疏像列另存以避免以后重复构建。
% 模块二：R_V*A*P与原Galerkin矩阵逐项对照；单位还对独立原NS差分矩阵核对。
% 模块三：检查完整模型细残差不增、固定映射叠加和单位稠密QR参考。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'CoarseLSTest:Overwrite','不得覆盖最小残差粗层研究。');
if strcmp(kind,'unit')
    sourceFile=fullfile(root,'outputs','v83_validation','coarseunit01','volume_coarse_unit_evidence.mat');
    s=load(sourceFile,'coarse','Q','baseInfo','gas','boundary','scale','rhs','fdMatrix');
    Q=s.Q;info=s.baseInfo;gas=s.gas;boundary=s.boundary;scale=s.scale;rhs=s.rhs;fdMatrix=s.fdMatrix;
else
    assert(strcmp(kind,'real'));sourceFile=fullfile(root,'outputs','v83_validation','coarserocket01','volume_coarse_real_evidence.mat');
    s=load(sourceFile,'coarse','rhs');rhs=s.rhs;
    d=load(fullfile(root,'outputs','v83_validation','krylov08','wall_velocity_krylov_evidence.mat'),'result','gas','boundary');
    Q=d.result.Q;gas=d.gas;boundary=d.boundary;scale=d.result.residualScale;info=struct;
    for key={'primitive','gradient','viscousFaceFlux'},info.(key{1})=d.result.info.(key{1});end;clear d;fdMatrix=[];
end
space=s.coarse.space;priorMatrix=s.coarse.matrix;assert(isequaln(Q,s.coarse.baseQ));clear s;cut=space.geometry;
before=memory;cache=buildCutCellSelectedDirectionCache(Q,info,cut,gas,boundary,scale);
coarse=buildRocketLeastSquaresCoarseOperator(cache,space,256*1024^2);afterCoarse=memory;ng=space.groupCount;
projected=full(kron(speye(5),space.restriction)*coarse.image*spdiags(1./coarse.columnLiftScale,0,5*ng,5*ng));
projectionError=norm(projected-priorMatrix,'fro')/norm(priorMatrix,'fro');assert(projectionError<1e-10,'CoarseLSTest:Projection','稀疏像未复现原完整空间Galerkin矩阵。');
fdError=NaN;if ~isempty(fdMatrix),fdError=norm(projected-fdMatrix,'fro')/norm(fdMatrix,'fro');assert(fdError<2e-6);end
graph=buildRocketShellGraphPreconditioner(Q,cut,gas,boundary,scale,1e-7,1);modelApplications=0;
initial=solveRocketShellGraphPreconditioner(graph,rhs);remaining=rhs+applyCutCellSelectedDirection(cache,initial);modelApplications=modelApplications+1;
t=tic;[update,linear]=solveRocketMinResidualCoarsePreconditioner(graph,coarse,cache,rhs);applySeconds=toc(t);modelApplications=modelApplications+linear.coarseCorrectionModelApplications;
after=rhs+applyCutCellSelectedDirection(cache,update);modelApplications=modelApplications+1;
assert(norm(after,'fro')<=norm(remaining,'fro')+1e-9*max(1,norm(remaining,'fro')),'CoarseLSTest:Nonincrease','完整细层模型残差不增检查失败。');
second=rhs.*reshape(cos((1:numel(rhs))*.47),size(rhs));[other,linear2]=solveRocketMinResidualCoarsePreconditioner(graph,coarse,cache,second);modelApplications=modelApplications+linear2.coarseCorrectionModelApplications;
[combined,linear3]=solveRocketMinResidualCoarsePreconditioner(graph,coarse,cache,.7*rhs-1.3*second);modelApplications=modelApplications+linear3.coarseCorrectionModelApplications;
expected=.7*update-1.3*other;linearity=norm(combined-expected,'fro')/max(1,norm(expected,'fro'));assert(linearity<1e-9,'CoarseLSTest:Linearity','固定最小残差粗层叠加检查失败。');
qrError=NaN;
if strcmp(kind,'unit')
    reference=full(coarse.image)\remaining(:);referenceUpdate=initial+space.prolongation*reshape(reference.*coarse.columnLiftScale,ng,5);
    qrError=norm(update-referenceUpdate,'fro')/max(1,norm(referenceUpdate,'fro'));assert(qrError<2e-7,'CoarseLSTest:QR','小Gram粗解不符合独立稠密QR参考。');
end
report=struct('status','passed-sparse-fine-L2-coarse-projection-model-nonincreasing-and-fixed-linearity','kind',kind, ...
    'fineCells',size(Q,1),'coarseGroups',ng,'coarseUnknowns',5*ng,'imageNonzeros',coarse.imageNonzeros,'imageBytes',coarse.imageBytes, ...
    'maximumImageBytes',coarse.maximumImageBytes,'coarseReciprocalConditionEstimate',coarse.reciprocalConditionEstimate, ...
    'volumeGalerkinProjectionRelativeDifference',projectionError,'originalNSCoarseFDRelativeDifference',fdError, ...
    'initialFineModelRelativeResidual',norm(remaining,'fro')/norm(rhs,'fro'),'correctedFineModelRelativeResidual',norm(after,'fro')/norm(rhs,'fro'), ...
    'fixedLinearityRelativeDifference',linearity,'denseQRReferenceRelativeDifference',qrError,'buildSeconds',coarse.buildSeconds,'twoLevelApplySeconds',applySeconds, ...
    'constructionModelApplications',coarse.constructionModelApplications,'otherValidationModelApplications',modelApplications, ...
    'matlabBytesBeforeCache',before.MemUsedMATLAB,'matlabBytesAfterCoarse',afterCoarse.MemUsedMATLAB,'peakMemoryMeasured',false, ...
    'actualOriginalNSEvaluations',0,'newAcceptedIterations',0,'physicalMeshChanged',false,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
mkdir(folder);save(fullfile(folder,'least_squares_coarse_evidence.mat'),'report','sourceFile','coarse','projected','priorMatrix','fdMatrix', ...
    'rhs','initial','remaining','update','after','linear','second','other','combined','-v7');
fid=fopen(fullfile(folder,'least_squares_coarse_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
