function report=studyRocketCurrentReservoirParent(folder)
%STUDYROCKETCURRENTRESERVOIRPARENT 在真实Q20重建近似父矩阵和局部因子，求新的修正方向。
% 模块一：只使用当前已验解析缓存，全部图系数/WLS同线块/五量局部逆重新构造。
% 模块二：保留旧方向作数值初猜，但用当前完整右端/原范数求校正，旧方向资格不继承。
% 模块三：C是预条件近似，1e-3目标保持；最终完整模型及原NS资格各自单列。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
cacheFile=fullfile(root,'outputs','v83_validation','reservoirjacreal01','riemann_selected_direction_real_evidence.mat');
s=load(cacheFile,'cache','rate','direction','value','report','sourceIdentity');cache=s.cache;
assert(isequaln(s.sourceIdentity,rocketSelectedDirectionSourceIdentity)&&s.report.currentBaseQExactlyReplayed);
sourceBinding=struct('sourceFile',s.report.sourceFile,'sourceFileSHA256',s.report.sourceSHA256,'stateLocation','result');
assert(strcmp(rocketStreamingFileSHA256(sourceBinding.sourceFile),sourceBinding.sourceFileSHA256));
Q=cache.Q;gas=cache.gas;boundary=cache.boundary;scale=cache.scale;rhs=s.rate./scale;initial=s.direction;initialModel=s.value;sourceIdentity=s.sourceIdentity;clear s;
timer=tic;graph=buildRocketShellGraphPreconditioner(Q,cache.geometry,gas,boundary,scale,1e-7,1);
operator=buildRocketWlsNormalPreconditioner(graph,cache,512*1024^2);clear graph;constructionSeconds=toc(timer);
operatorFile=fullfile(folder,'current_operator_evidence.mat');save(operatorFile,'operator','Q','gas','boundary','scale','rhs','initial','initialModel','cacheFile','sourceBinding','sourceIdentity','-v7');
packed=packRocketSteadyGraph(operator,struct('ordering','symrcm'));clear operator cache;
packedFile=fullfile(folder,'current_packed_evidence.mat');save(packedFile,'packed','operatorFile','sourceBinding','sourceIdentity','-v7');
seedImage=applyRocketPackedGraph(packed,initial);
seedGraphDifference=norm(seedImage+initialModel,'fro')/max(1,norm(initialModel,'fro'));
partition=buildRocketStrongPatches(packed,8);
factor=factorRocketCoupledPatches(packed,partition,struct('failureEvidenceFile',fullfile(folder,'current_factor_failure.mat')));clear partition;
factorFile=fullfile(folder,'current_factor_evidence.mat');save(factorFile,'factor','packedFile','sourceBinding','sourceIdentity','-v7');
[columnScale,columnNorm,columnConstruction]=buildRocketParentImageScale(packed,factor);
scaleFile=fullfile(folder,'current_image_scale_evidence.mat');save(scaleFile,'columnScale','columnNorm','columnConstruction','packedFile','factorFile','sourceBinding','sourceIdentity','-v7');
storage=whos('packed','factor','columnScale','columnNorm');cacheBytes=sum([storage.bytes]);assert(cacheBytes<=1024*1024^2);
[value,detail,coordinates]=solveRocketParentLeastSquares(packed,factor,rhs,initial, ...
    struct('maximumIterations',256,'relativeTolerance',1e-3),columnScale);
image=applyRocketPackedGraph(packed,value);remaining=rhs-image;
save(fullfile(folder,'current_parent_candidate_evidence.mat'),'value','detail','coordinates','image','remaining','rhs','initial','initialModel', ...
    'seedGraphDifference','seedImage','sourceBinding','sourceIdentity','cacheFile','operatorFile','packedFile','factorFile','scaleFile','cacheBytes','storage','-v7');
factorReport=factor.factorizationReport;clear packed factor columnScale columnNorm;
c=load(cacheFile,'cache');[fullModel,fullDetail]=applyCutCellSelectedDirection(c.cache,value);clear c;
fullRemaining=rhs+fullModel;fullResidual=norm(fullRemaining,'fro')/norm(rhs,'fro');
report=struct('status','tested-fresh-current-Q20-graph-local-factor-image-scale-and-new-parent-candidate', ...
    'sourceBinding',sourceBinding,'currentStateCoefficientsRebuilt',true,'oldNumericMatrixOrFactorsReused',false, ...
    'graphConstructionSeconds',constructionSeconds,'seedGraphToCompleteModelDifference',seedGraphDifference, ...
    'factorReport',factorReport,'columnConstruction',columnConstruction,'innerReport',detail,'actualInputCacheArrayBytes',cacheBytes, ...
    'fullSelectedModelOriginalRHSRelativeResidual',fullResidual,'modelPointOneMet',fullResidual<=.1, ...
    'CInternalQualificationNotPromoted',true,'originalNSDirectionCheckedHere',false,'newSelectedModelApplications',1, ...
    'newActualNSEvaluations',0,'newFlowAcceptances',0,'physicalStateChanged',false,'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
save(fullfile(folder,'current_parent_model_probe.mat'),'value','fullModel','fullDetail','fullRemaining','report','sourceBinding','sourceIdentity','-v7');
fid=fopen(fullfile(folder,'current_parent_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
