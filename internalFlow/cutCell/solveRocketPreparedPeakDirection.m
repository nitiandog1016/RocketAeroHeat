function report=solveRocketPreparedPeakDirection(folder,preparationFile,options,previousCandidateFile)
%SOLVEROCKETPREPAREDPEAKDIRECTION 独立求解已绑定当前状态的方向，首次始终零初猜。
% 模块一：准备系统不含物理接受资格；所有产物/源状态/算法SHA在求解前重查。
% 模块二：同一系统同一惩罚可从保存解重启，显式记录各段，禁止悄悄改变目标。
% 模块三：完整模型0.1、原C1e-3与实际NS资格分列，不自动更新物理状态。
if nargin<3,options=struct;end
if nargin<4,previousCandidateFile='';end
defaults=struct('maximumIterations',128,'regularizationStrength',.1);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'PeakSystem:Options','未知方向求解选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
folder=char(java.io.File(folder).getCanonicalPath());preparationFile=char(java.io.File(preparationFile).getCanonicalPath());
assert(~exist(folder,'dir'),'PeakSystem:Overwrite','不得覆盖方向检查点。');p=loadRocketPeakSystem(preparationFile);preparationFileSHA256=rocketStreamingFileSHA256(preparationFile);
m=load(p.artifacts.metric.file);initial=zeros(size(m.rhs));segments=struct([]);previousDetail=struct([]);
if ~isempty(previousCandidateFile)
    previousCandidateFile=char(java.io.File(previousCandidateFile).getCanonicalPath());old=load(previousCandidateFile);
    assert(strcmp(old.preparationFileSHA256,preparationFileSHA256)&& ...
        old.detail.options.regularizationStrength==options.regularizationStrength&&old.detail.beatsZeroDirectionObjective&& ...
        old.detail.normalizedAugmentedObjective<1,'PeakSystem:Restart','重启检查点不属于同一系统/同一目标。');
    initial=old.value;segments=old.segments;previousDetail=old.detail;clear old;
end
packedFile=p.artifacts.packed.file;factorFile=p.artifacts.factor.file;scaleFile=p.artifacts.scale.file;cacheFile=p.artifacts.cache.file;
a=load(packedFile,'packed');b=load(factorFile,'factor');c=load(scaleFile,'columnScale');mkdir(folder);
[value,detail,coordinates]=solveRocketRegularizedParentDirection(a.packed,b.factor,m.rhs,initial,c.columnScale,m.metric, ...
    struct('maximumIterations',options.maximumIterations,'regularizationStrength',options.regularizationStrength,'metricNormalization',p.metricNormalization));
segment=struct('nativeIterations',detail.nativeIterations,'initialObjective',detail.initialNormalizedAugmentedObjective, ...
    'finalObjective',detail.normalizedAugmentedObjective,'relativeCResidual',detail.relativeCResidual);
if isempty(segments),segments=segment;else,segments(end+1)=segment;end
rhs=m.rhs;metric=m.metric;baseMetric=m.baseMetric;unweightedFinal=applyRocketPrimitiveUpdateMetric(baseMetric,value);
sourceBinding=p.sourceBinding;sourceIdentity=p.sourceIdentity;geometryBinding=p.geometryBinding;
save(fullfile(folder,'peak_zero_candidate_evidence.mat'),'value','detail','coordinates','initial','rhs','metric','baseMetric','unweightedFinal', ...
    'segments','previousDetail','previousCandidateFile','preparationFile','preparationFileSHA256','sourceBinding','sourceIdentity','geometryBinding', ...
    'cacheFile','packedFile','factorFile','scaleFile','-v7');clear a b c m;
c=load(cacheFile,'cache');[fullModel,fullDetail]=applyCutCellSelectedDirection(c.cache,value);clear c;
restartCount=numel(segments)-1;firstIterations=segments(1).nativeIterations;restartedIterations=sum([segments(2:end).nativeIterations]);
directionPreparation=struct('kind','current-state-peak-soft-zero-LSQR-v1','preparationFile',preparationFile, ...
    'preparationFileSHA256',preparationFileSHA256,'referenceReusedAsInitialGuess',false,'restartCount',restartCount, ...
    'firstSegmentIterations',firstIterations,'restartedSegmentIterations',restartedIterations, ...
    'totalNativeIterations',sum([segments.nativeIterations]),'hardUpdateLimitClaimed',false);
report=struct('status','tested-prepared-current-state-peak-soft-zero-start-direction','sourceBinding',sourceBinding, ...
    'sourceIdentity',sourceIdentity,'geometryBinding',geometryBinding,'linearReport',detail,'segments',segments, ...
    'directionPreparation',directionPreparation,'unrestartedLSQRClaimed',restartCount==0, ...
    'completeSelectedModelRelativeDefect',norm(rhs+fullModel,'fro')/norm(rhs,'fro'), ...
    'maximumRawPrimitiveComponents',max(abs(unweightedFinal),[],1),'rawPrimitiveUpdateNorm',norm(unweightedFinal,'fro'), ...
    'newSelectedModelApplications',1,'newActualNSEvaluations',0,'newFlowAcceptances',0,'physicalStateChanged',false, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
report.modelPointOneMet=report.completeSelectedModelRelativeDefect<=.1;
save(fullfile(folder,'peak_zero_model_evidence.mat'),'value','fullModel','fullDetail','report','-v7');
fid=fopen(fullfile(folder,'peak_zero_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
fprintf('PEAK_PREPARED_DIRECTION model=%.12g C=%.12g iterations=%d restartCount=%d\n',report.completeSelectedModelRelativeDefect,detail.relativeCResidual,detail.nativeIterations,restartCount);
end
