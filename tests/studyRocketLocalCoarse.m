function report=studyRocketLocalCoarse(folder,kind,numberOfCells)
%STUDYROCKETLOCALCOARSE 对高剩余原CV构造固定局部方向，核对内存、投影、完整作用及小域原NS。
% 模块一：真实场直接读取krylov12保存剩余，不新跑或冒称接受非线性步。
% 模块二：所有原CV/方程保留，局部像只用于求解加速；按原面/层保存明确身份。
% 模块三：小域原NS中心差分与真实区域的模型投影分列，均不计新增合格面积。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'LocalStudy:Overwrite','不得覆盖局部研究。');mkdir(folder);assert(ismember(kind,{'unit','real'}));
if strcmp(kind,'unit'),name='mixedunit01';else,name='mixedrocket02';end
parentFile=fullfile(root,'outputs','v83_validation',name,'mixed_coarse_evidence.mat');s=load(parentFile,'coarse');parent=s.coarse;clear s;
Q=parent.baseQ;cut=parent.space.geometry;gas=parent.gas;boundary=parent.boundary;scale=parent.scale;calls=0;sourceFile='';
if strcmp(kind,'unit')
    [rate,info]=cutCellNavierStokesResidual(Q,cut,gas,boundary);calls=calls+1;rhs=rate./scale;remaining=reshape(sin((1:numel(Q))*.37),size(Q));
    [~,order]=sort(sum(remaining.^2,2),'descend');
else
    sourceFile=fullfile(root,'outputs','v83_validation','krylov12','wall_velocity_krylov_evidence.mat');s=load(sourceFile,'result');
    assert(isequaln(s.result.Q,Q)&&isequaln(s.result.Q,s.result.initialQ));rhs=s.result.initialRate./scale;
    remaining=rhs+s.result.linearLog{end}.steadyDirectionEvidence.value;info=struct;
    for key={'primitive','gradient','viscousFaceFlux'},info.(key{1})=s.result.initialInfo.(key{1});end;clear s;
    selectionFile=fullfile(root,'outputs','v83_validation','localcells01','local_cells_selection.mat');s=load(selectionFile,'cellOrder','remaining');
    assert(isequaln(remaining,s.remaining));order=s.cellOrder;clear s;
end
validateattributes(numberOfCells,{'numeric'},{'real','finite','scalar','integer','positive','<=',numel(order)});cellIDs=int64(order(1:numberOfCells));
cache=buildCutCellSelectedDirectionCache(Q,info,cut,gas,boundary,scale);clear info;
[oldCorrection,oldLinear]=solveRocketMixedCoarseSystem(parent,remaining);oldPredicted=remaining+applyCutCellSelectedDirection(cache,oldCorrection);otherModel=1;
coarse=buildRocketLocalCoarseOperator(parent,cache,remaining,cellIDs,struct('maximumAdditionalColumns',5*numberOfCells));
% 先保存完整构造端点；后续投影/差分检查失败时仍可审计已有像，不能只留下运行日志。
save(fullfile(folder,'local_coarse_construction_evidence.mat'),'coarse','parentFile','sourceFile','rhs','remaining','cellIDs','oldCorrection','oldPredicted','oldLinear','-v7');
clear parent;[correction,linear]=solveRocketLocalCoarseSystem(coarse,remaining);predicted=remaining+applyCutCellSelectedDirection(cache,correction);otherModel=otherModel+1;
weights=coarse.upper\(coarse.upper'\(coarse.image'*remaining(:)));expected=remaining(:)-coarse.image*weights;
projectionError=norm(predicted(:)-expected)/max(1,norm(remaining(:)));assert(projectionError<1e-8,'LocalStudy:Projection','局部组合提升没有重放完整当前模型。');
probeCoefficients=sin((1:size(coarse.image,2))'*.27);probeDirection=liftRocketLocalCoarse(coarse,probeCoefficients);
probeActual=-applyCutCellSelectedDirection(cache,probeDirection);otherModel=otherModel+1;probeProduct=coarse.image*probeCoefficients;
probeError=norm(probeActual(:)-probeProduct)/max(1,norm(probeProduct));assert(probeError<1e-9,'LocalStudy:Probe','局部逐列提升和完整作用不符。');
secondRemaining=remaining.*reshape(cos((1:numel(remaining))*.43),size(remaining));other=solveRocketLocalCoarseSystem(coarse,secondRemaining);
combined=solveRocketLocalCoarseSystem(coarse,.7*remaining-1.3*secondRemaining);
linearity=norm(combined-(.7*correction-1.3*other),'fro')/max(1,norm(combined,'fro'));assert(linearity<1e-8,'LocalStudy:Linearity','局部冻结预条件不满足右端叠加。');
validationDirections={};validationPlus={};validationMinus={};validationModel={};fdErrors=[];h=1e-6;
directions=1;if strcmp(kind,'unit'),directions=2;end
    for j=1:directions
        coefficients=zeros(size(coarse.image,2),1);
        if strcmp(kind,'real'),coefficients=weights;
        elseif j==1,coefficients(numel(coarse.retainedMixedColumns)+1:end)=cos((1:numel(coarse.localCellIDs))'*.23);
        else,coefficients=sin((1:size(coarse.image,2))'*.31);end
        direction=liftRocketLocalCoarse(coarse,coefficients);assert(any(direction~=0,'all'),'LocalStudy:Direction','没有有效的局部方向。');
        direction=direction/max(abs(direction),[],'all');
        plus=cutCellNavierStokesResidual(Q+h*direction.*scale,cut,gas,boundary);minus=cutCellNavierStokesResidual(Q-h*direction.*scale,cut,gas,boundary);calls=calls+2;
        measured=(plus-minus)./scale/(2*h);model=applyCutCellSelectedDirection(cache,direction);otherModel=otherModel+1;
        error=norm(measured-model,'fro')/max(1,norm(measured,'fro'));assert(error<2e-6,'LocalStudy:NS','局部方向与原NS中心差分不符。');
        validationDirections{j}=direction;validationPlus{j}=plus;validationMinus{j}=minus;validationModel{j}=model;fdErrors(end+1)=error; %#ok<AGROW>
end
report=struct('status','tested-original-CV-local-conservative-enrichment-with-fixed-lift-and-original-model-replay','kind',kind, ...
    'candidateCells',numberOfCells,'candidateColumns',numel(coarse.candidateCellIDs),'selectedLocalColumns',numel(coarse.localCellIDs), ...
    'retainedMixedColumns',numel(coarse.retainedMixedColumns),'removedMixedColumns',coarse.removedMixedColumns,'coarseUnknowns',size(coarse.image,2), ...
    'imageBytes',coarse.imageBytes,'maximumImageBytes',coarse.maximumImageBytes,'imageNonzeros',coarse.imageNonzeros, ...
    'reciprocalConditionEstimate',coarse.reciprocalConditionEstimate,'selectionStopReason',coarse.selectionStopReason, ...
    'originalRemainingRelativeDefect',norm(remaining,'fro')/norm(rhs,'fro'),'oldMixedProjectedRelativeDefect',norm(oldPredicted,'fro')/norm(rhs,'fro'), ...
    'localProjectedRelativeDefect',norm(predicted,'fro')/norm(rhs,'fro'),'extraSquaredNormFractionRemoved',1-norm(predicted,'fro')^2/norm(oldPredicted,'fro')^2, ...
    'fullModelProjectionRelativeDifference',projectionError,'probeImageRelativeDifference',probeError,'fixedLinearityRelativeDifference',linearity, ...
    'originalNSDirectionFDRelativeDifferences',fdErrors,'actualOriginalNSEvaluations',calls, ...
    'constructionModelApplications',coarse.constructionModelApplications,'otherModelApplications',otherModel,'buildSeconds',coarse.buildSeconds, ...
    'physicalMeshChanged',false,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'local_coarse_evidence.mat'),'report','coarse','parentFile','sourceFile','rhs','remaining','cellIDs','oldCorrection','oldPredicted','oldLinear', ...
    'correction','predicted','linear','probeCoefficients','probeDirection','probeActual','probeProduct','secondRemaining','other','combined', ...
    'validationDirections','validationPlus','validationMinus','validationModel','h','-v7');
fid=fopen(fullfile(folder,'local_coarse_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
