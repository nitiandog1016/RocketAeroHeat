function report=studyRocketCoarseLayerEnrichment(folder,kind)
%STUDYROCKETCOARSELAYERENRICHMENT 将近壁层分别聚合，先测实际失败剩余的新增投影收益。
% 模块一：保留旧坐标分箱和真实正面连通关系，只细分法向层；核对旧粗空间严格嵌套。
% 模块二：构造当前完整细层作用的有界稀疏像；不增加物理网格或更新流场。
% 模块三：小域用原NS两个粗方向中心差分校核；真实域只用已保存端点及失败方向。
% 模块四：比较旧/新粗基投影、内存及条件数，收益不够时不应直接增加外层方向预算。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'LayerCoarse:Overwrite','不得覆盖逐层粗空间研究。');
assert(ismember(kind,{'unit','real'}),'LayerCoarse:Kind','研究类型须为unit或real。');
if strcmp(kind,'unit'),parentName='coarselsunit01';else,parentName='coarselsrocket01';end
parentFile=fullfile(root,'outputs','v83_validation',parentName,'least_squares_coarse_evidence.mat');
s=load(parentFile,'coarse');old=s.coarse;clear s;cut=old.space.geometry;Q=old.baseQ;gas=old.gas;boundary=old.boundary;scale=old.scale;
calls=0;sourceFile='';validationDirections={};validationPlus={};validationMinus={};validationModel={};h=1e-6;
if strcmp(kind,'unit')
    [rate,info]=cutCellNavierStokesResidual(Q,cut,gas,boundary);calls=calls+1;rhs=rate./scale;
    failureRemaining=reshape(sin((1:numel(Q))*.37),size(Q));
else
    sourceFile=fullfile(root,'outputs','v83_validation','krylov09','wall_velocity_krylov_evidence.mat');
    d=load(sourceFile,'result');assert(isequaln(Q,d.result.Q)&&isequaln(d.result.initialQ,d.result.Q));rhs=d.result.initialRate./scale;
    evidence=d.result.linearLog{end}.steadyDirectionEvidence;failureRemaining=rhs+evidence.value;
    info=struct;for key={'primitive','gradient','viscousFaceFlux'},info.(key{1})=d.result.initialInfo.(key{1});end;clear d evidence;
end
cache=buildCutCellSelectedDirectionCache(Q,info,cut,gas,boundary,scale);clear info;
[old,reuseReport]=loadRocketLeastSquaresCoarseOperator(parentFile,cache,old.space,256*1024^2);
modelApplications=reuseReport.reuseVerificationModelApplications;
options=old.space.options;nl=size(Q,1)/size(cut.inputGeometry.F,1);
options.layerBreaks=unique([options.layerBreaks 1:min(3,nl-2)]);options.maximumCoarseUnknowns=500;
space=buildRocketVolumeCoarseSpace(cut,options);parentTransfer=space.restriction*old.space.prolongation;
nesting=norm(space.prolongation*parentTransfer-old.space.prolongation,Inf);
assert(nesting<1e-12,'LayerCoarse:Nested','逐层空间没有保留原连通粗空间。');
before=memory;coarse=buildRocketLeastSquaresCoarseOperator(cache,space,256*1024^2);after=memory;
% 按列核对旧像是新像的线性组合，避免同时构造两份大稀疏父像。
parentImageDifference=0;parentImageNorm=0;ng=space.groupCount;oldNg=old.space.groupCount;
for b=1:5,for g=1:oldNg
    column=(b-1)*oldNg+g;weights=zeros(5*ng,1);ids=(b-1)*ng+(1:ng);
    weights(ids)=parentTransfer(:,g)*old.columnLiftScale(column)./coarse.columnLiftScale(ids);
    difference=coarse.image*weights-old.image(:,column);
    parentImageDifference=parentImageDifference+sum(difference.^2);parentImageNorm=parentImageNorm+sum(old.image(:,column).^2);
end,end
imageNesting=sqrt(parentImageDifference/parentImageNorm);assert(imageNesting<1e-9,'LayerCoarse:ImageNested','新粗像没有保留旧完整方向空间。');
[oldCorrection,oldReport]=solveRocketLeastSquaresCoarseSystem(old,failureRemaining);
[newCorrection,newReport]=solveRocketLeastSquaresCoarseSystem(coarse,failureRemaining);
oldPredicted=failureRemaining+applyCutCellSelectedDirection(cache,oldCorrection);newPredicted=failureRemaining+applyCutCellSelectedDirection(cache,newCorrection);modelApplications=modelApplications+2;
assert(norm(newPredicted,'fro')<=norm(oldPredicted,'fro')+1e-8*max(1,norm(oldPredicted,'fro')), ...
    'LayerCoarse:Projection','严格嵌套后的细层模型最小残差反而放大。');
if strcmp(kind,'unit')
    for j=1:2
        coefficients=reshape(sin((1:5*ng)*(.17+.11*j)),ng,5);direction=space.prolongation*coefficients;
        plus=cutCellNavierStokesResidual(Q+h*direction.*scale,cut,gas,boundary);
        minus=cutCellNavierStokesResidual(Q-h*direction.*scale,cut,gas,boundary);calls=calls+2;
        predicted=applyCutCellSelectedDirection(cache,direction);modelApplications=modelApplications+1;
        measured=(plus-minus)./scale/(2*h);error=norm(measured-predicted,'fro')/max(1,norm(measured,'fro'));
        assert(error<2e-6,'LayerCoarse:OriginalNS','逐层粗方向与原NS中心差分不一致。');
        validationDirections{j}=direction;validationPlus{j}=plus;validationMinus{j}=minus;validationModel{j}=predicted;
    end
end
graph=buildRocketShellGraphPreconditioner(Q,cut,gas,boundary,scale,1e-7,1);
[update,linear]=solveRocketMinResidualCoarsePreconditioner(graph,coarse,cache,rhs);modelApplications=modelApplications+1;
afterGraphCoarse=rhs+applyCutCellSelectedDirection(cache,update);modelApplications=modelApplications+1;
report=struct('status','tested-layer-enriched-nested-coarse-image-and-model-projection-not-physical-solution','kind',kind, ...
    'fineCells',size(Q,1),'layerBreaks',options.layerBreaks,'parentCoarseGroups',oldNg,'coarseGroups',ng,'coarseUnknowns',5*ng, ...
    'parentTransferInfinityError',nesting,'parentImageRelativeDifference',imageNesting,'imageNonzeros',coarse.imageNonzeros, ...
    'imageBytes',coarse.imageBytes,'maximumImageBytes',coarse.maximumImageBytes,'reciprocalConditionEstimate',coarse.reciprocalConditionEstimate, ...
    'failedRemainingNormRelativeToOriginalRHS',norm(failureRemaining,'fro')/norm(rhs,'fro'), ...
    'parentPredictedRemainingRelativeToOriginalRHS',norm(oldPredicted,'fro')/norm(rhs,'fro'), ...
    'enrichedPredictedRemainingRelativeToOriginalRHS',norm(newPredicted,'fro')/norm(rhs,'fro'), ...
    'additionalRemainingSquaredNormFractionRemoved',1-norm(newPredicted,'fro')^2/norm(oldPredicted,'fro')^2, ...
    'graphAndEnrichedCoarseModelRelativeResidual',norm(afterGraphCoarse,'fro')/norm(rhs,'fro'), ...
    'constructionModelApplications',coarse.constructionModelApplications,'otherModelApplications',modelApplications, ...
    'actualOriginalNSEvaluations',calls,'buildSeconds',coarse.buildSeconds,'matlabBytesBeforeBuild',before.MemUsedMATLAB, ...
    'matlabBytesAfterBuild',after.MemUsedMATLAB,'peakMemoryMeasured',false,'originalPhysicalMeshChanged',false, ...
    'actualNonlinearIterations',0,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
mkdir(folder);save(fullfile(folder,'coarse_layer_evidence.mat'),'report','sourceFile','parentFile','coarse','parentTransfer','rhs','failureRemaining', ...
    'oldCorrection','newCorrection','oldPredicted','newPredicted','oldReport','newReport','update','linear','afterGraphCoarse', ...
    'validationDirections','validationPlus','validationMinus','validationModel','h','-v7');
fid=fopen(fullfile(folder,'coarse_layer_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
