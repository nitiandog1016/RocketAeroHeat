function report=inspectRocketTemperatureLookaheadRank(folder,stepFolder)
%INSPECTROCKETTEMPERATURELOOKAHEADRANK 检查前瞻拒绝的当前约束，不重算NS或自由方向。
% 模块一：绑定真实Q37、当前准备文件与已保存自由前缀，逐值核对原10K行。
% 模块二：分别观察原10K、全补充带和原行加全部预测失败单元所需面的秩。
% 模块三：物理坐标QR/SVD仅作诊断；是否可用仍由原QR右坐标投影器原门槛决定。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
contextFile=fullfile(stepFolder,'direction_context.mat');prefixFile=fullfile(stepFolder,'linear_checkpoints/lookahead_free_prefix.mat');
c=load(contextFile,'factor','initialRows','p','preparationFile','preparationSHA256');
f=load(prefixFile,'free','freeDetail','freeCoordinates','originalRows');
assert(isequaln(f.originalRows,c.initialRows));
assert(strcmp(rocketStreamingFileSHA256(c.preparationFile),c.preparationSHA256));
p=c.p;assert(strcmp(rocketStreamingFileSHA256(p.sourceBinding.sourceFile),p.sourceBinding.sourceFileSHA256));
assert(strcmp(rocketStreamingFileSHA256(p.geometryBinding.file),p.geometryBinding.sha256));
s=load(p.sourceBinding.sourceFile,'result','gas');g=load(p.geometryBinding.file,p.geometryBinding.variable);cut=g.(p.geometryBinding.variable);
[allRows,selection]=selectRocketTemperatureLookaheadRows(s.result.Q,s.result.info,cut,s.gas,s.result.residualScale,f.free,1/8);
sets={c.initialRows.faceIDs,allRows.faceIDs,unique([c.initialRows.faceIDs;selection.requiredOriginalStateFaceIDs])};
names={'original-current-10K','all-extended-band','original-plus-all-predicted-required-faces'};cases=cell(3,1);
for j=1:3
    [found,slots]=ismember(sets{j},allRows.faceIDs);assert(all(found));rows=allRows;
    rows.H=allRows.H(slots,:);rows.faceIDs=allRows.faceIDs(slots);rows.temperatureBoundSigns=allRows.temperatureBoundSigns(slots);
    tangent=rows.H*f.free(:);roundoff=128*eps*full(abs(rows.H)*abs(f.free(:)));
    active=find(rows.temperatureBoundSigns.*tangent>roundoff);selected=rows;selected.H=rows.H(active,:);selected.faceIDs=rows.faceIDs(active);
    Ht=selected.H';representatives=zeros(0,1);mapping=zeros(numel(active),1);
    for k=1:numel(active)
        match=0;for a=1:numel(representatives),if isequal(Ht(:,k),Ht(:,representatives(a))),match=a;break;end,end
        if match==0,representatives(end+1,1)=k;match=numel(representatives);end %#ok<AGROW>
        mapping(k)=match;
    end
    support=find(any(Ht,2));B=full(Ht(support,representatives));B=B./vecnorm(B,2,1);
    [~,R]=qr(B,0);physicalPivots=abs(diag(R));singularValues=svd(B);
    projector=[];identifier='';message='';
    try,projector=buildRocketTemperatureProjector(c.factor,selected);
    catch e,identifier=e.identifier;message=e.message;if ~strcmp(identifier,'TemperatureProjection:Rank'),rethrow(e);end,end
    result=struct('name',names{j},'originalFaceRows',numel(slots),'outwardActiveRows',numel(active), ...
        'exactDistinctActiveRows',numel(representatives),'physicalSupportUnknowns',numel(support), ...
        'physicalNormalizedQRDiagonal',physicalPivots,'physicalNormalizedSingularValues',singularValues, ...
        'originalRightProjectorPassed',isempty(identifier),'rejectionIdentifier',identifier,'rejectionMessage',message, ...
        'allOriginal10KFaceIDsRetained',all(ismember(c.initialRows.faceIDs,rows.faceIDs)), ...
        'allPredictedRequiredFaceIDsRetained',all(ismember(selection.requiredOriginalStateFaceIDs,rows.faceIDs)), ...
        'newActualNSEvaluations',0,'newFlowAcceptances',0,'newProductionSolvedAreaM2',0);
    cases{j}=struct('report',result,'rows',rows,'activeRowIndices',active,'representatives',representatives,'originalToUnique',mapping, ...
        'projector',projector,'physicalCompactColumns',B,'physicalSupport',support);
end
report=struct('status','observed-current-lookahead-rank-refusal-and-three-explicit-numerical-selection-sets', ...
    'sourceStateFile',p.sourceBinding.sourceFile,'sourceStateSHA256',p.sourceBinding.sourceFileSHA256, ...
    'contextSHA256',rocketStreamingFileSHA256(contextFile),'freePrefixSHA256',rocketStreamingFileSHA256(prefixFile), ...
    'newActualNSEvaluations',0,'newFlowAcceptances',0,'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'temperature_lookahead_rank.mat'),'cases','selection','report','contextFile','prefixFile','-v7');
summary=struct('report',report,'selection',selection,'cases',{cellfun(@(x)x.report,cases,'UniformOutput',false)});
fid=fopen(fullfile(folder,'temperature_lookahead_rank.json'),'w','n','UTF-8');assert(fid>=0);closer=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(summary));disp(report);
end
