function report=studyRocketPeakZeroStart(folder)
%STUDYROCKETPEAKZEROSTART 保持当前Q22完整C/因子/权重/惩罚不变，从零方向求解。
% 模块一：旧256步结果的目标大于已知零方向，不能当有效候选；零初猜直接避免抵消巨量旧像。
% 模块二：显式固定同一metricNormalization，使penalty逐值相同；并不变更原物理/模型门槛。
% 模块三：完整C、原始状态幅度与完整模型分别记录，不自动执行NS或提交流场。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);base=fullfile(root,'outputs','v83_validation','peakweightreal01');
s=load(fullfile(base,'current_peak_candidate_evidence.mat'));assert(strcmp(rocketStreamingFileSHA256(s.sourceBinding.sourceFile),s.sourceBinding.sourceFileSHA256));
assert(isequaln(s.sourceIdentity,rocketSelectedDirectionSourceIdentity));
p=load(s.packedFile,'packed');f=load(s.factorFile,'factor');c=load(s.scaleFile,'columnScale');initial=zeros(size(s.initial));
[value,detail,coordinates]=solveRocketRegularizedParentDirection(p.packed,f.factor,s.rhs,initial,c.columnScale,s.metric, ...
    struct('maximumIterations',128,'regularizationStrength',.1,'metricNormalization',s.detail.metricNormalization));
assert(detail.penaltyInOriginalResidualUnits==s.detail.penaltyInOriginalResidualUnits);
rhs=s.rhs;metric=s.metric;baseMetric=s.baseMetric;sourceBinding=s.sourceBinding;sourceIdentity=s.sourceIdentity;cacheFile=s.cacheFile;
unweightedFinal=applyRocketPrimitiveUpdateMetric(baseMetric,value);previousDetail=s.detail;packedFile=s.packedFile;factorFile=s.factorFile;scaleFile=s.scaleFile;
save(fullfile(folder,'peak_zero_candidate_evidence.mat'),'value','detail','coordinates','initial','rhs','metric','baseMetric','unweightedFinal', ...
    'previousDetail','sourceBinding','sourceIdentity','cacheFile','packedFile','factorFile','scaleFile','-v7');clear p f c s;
c=load(cacheFile,'cache');[fullModel,fullDetail]=applyCutCellSelectedDirection(c.cache,value);clear c;
report=struct('status','tested-current-Q22-same-augmented-objective-zero-start-direction', ...
    'sourceBinding',sourceBinding,'sourceIdentity',sourceIdentity,'linearReport',detail, ...
    'samePenaltyAsWarmRun',true,'warmRunNormalizedObjective',previousDetail.normalizedAugmentedObjective, ...
    'knownZeroDirectionNormalizedObjective',1,'completeSelectedModelRelativeDefect',norm(rhs+fullModel,'fro')/norm(rhs,'fro'), ...
    'maximumRawPrimitiveComponents',max(abs(unweightedFinal),[],1),'rawPrimitiveUpdateNorm',norm(unweightedFinal,'fro'), ...
    'newSelectedModelApplications',1,'newActualNSEvaluations',0,'newFlowAcceptances',0, ...
    'physicalStateChanged',false,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
report.modelPointOneMet=report.completeSelectedModelRelativeDefect<=.1;
save(fullfile(folder,'peak_zero_model_evidence.mat'),'value','fullModel','fullDetail','report','-v7');
fid=fopen(fullfile(folder,'peak_zero_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
