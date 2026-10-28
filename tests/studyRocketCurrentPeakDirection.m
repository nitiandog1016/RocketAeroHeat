function report=studyRocketCurrentPeakDirection(folder,sourceCase)
%STUDYROCKETCURRENTPEAKDIRECTION 在当前真实端点重建全算子并研究局部峰值加权方向。
% 模块一：完整NS精确回放当前Q/rate/info，重建当前缓存、图/WLS块、局部逆及像列。
% 模块二：旧方向仅作数值初猜及冻结权重参考，旧矩阵/导数/探针资格不继承。
% 模块三：峰值软惩罚不保证硬更新上限；最终原C、完整模型和实际NS资格分开。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(ischar(sourceCase)&&~isempty(regexp(sourceCase,'^krylov[0-9]+$','once'))&&~exist(folder,'dir'));mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation',sourceCase,'wall_velocity_krylov_evidence.mat');s=load(sourceFile,'result','gas','boundary');
Q=s.result.Q;gas=s.gas;boundary=s.boundary;scale=s.result.residualScale;
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');cut=d.domain;
[rate,info]=cutCellNavierStokesResidual(Q,cut,gas,boundary);assert(isequaln(rate,s.result.rate)&&isequaln(info,s.result.info));clear s;
rhs=rate./scale;
seedFile=fullfile(root,'outputs','v83_validation','updateradiusreal01','regularized_parent_candidate_evidence.mat');seed=load(seedFile,'value');initial=seed.value;clear seed;
sourceBinding=struct('sourceFile',sourceFile,'sourceFileSHA256',rocketStreamingFileSHA256(sourceFile),'stateLocation','result');
sourceIdentity=rocketSelectedDirectionSourceIdentity;timer=tic;cache=buildCutCellSelectedDirectionCache(Q,info,cut,gas,boundary,scale);
cacheFile=fullfile(folder,'current_peak_cache_evidence.mat');save(cacheFile,'cache','Q','rate','scale','gas','boundary','sourceBinding','sourceIdentity','seedFile','-v7');
graph=buildRocketShellGraphPreconditioner(Q,cache.geometry,gas,boundary,scale,1e-7,1);
operator=buildRocketWlsNormalPreconditioner(graph,cache,512*1024^2);clear graph info cut d cache;
operatorFile=fullfile(folder,'current_peak_operator_evidence.mat');save(operatorFile,'operator','Q','rhs','sourceBinding','sourceIdentity','-v7');
packed=packRocketSteadyGraph(operator,struct('ordering','symrcm'));clear operator;
packedFile=fullfile(folder,'current_peak_packed_evidence.mat');save(packedFile,'packed','operatorFile','sourceBinding','sourceIdentity','-v7');
partition=buildRocketStrongPatches(packed,8);factor=factorRocketCoupledPatches(packed,partition,struct('failureEvidenceFile',fullfile(folder,'current_peak_factor_failure.mat')));clear partition;
factorFile=fullfile(folder,'current_peak_factor_evidence.mat');save(factorFile,'factor','sourceBinding','sourceIdentity','-v7');
[columnScale,columnNorm,columnConstruction]=buildRocketParentImageScale(packed,factor);
scaleFile=fullfile(folder,'current_peak_scale_evidence.mat');save(scaleFile,'columnScale','columnNorm','columnConstruction','sourceBinding','sourceIdentity','-v7');
baseMetric=buildRocketPrimitiveUpdateMetric(Q,scale,gas);[metric,weighting]=reweightRocketPrimitiveUpdateMetric(baseMetric,initial);
rhs=rate./scale;constructionSeconds=toc(timer);
[value,detail,coordinates]=solveRocketRegularizedParentDirection(packed,factor,rhs,initial,columnScale,metric, ...
    struct('maximumIterations',256,'regularizationStrength',.1));
unweightedInitial=applyRocketPrimitiveUpdateMetric(baseMetric,initial);unweightedFinal=applyRocketPrimitiveUpdateMetric(baseMetric,value);
factorReport=factor.factorizationReport;clear packed factor columnScale columnNorm;
save(fullfile(folder,'current_peak_candidate_evidence.mat'),'value','detail','coordinates','metric','baseMetric','weighting','rhs','initial', ...
    'unweightedInitial','unweightedFinal','sourceBinding','sourceIdentity','cacheFile','packedFile','factorFile','scaleFile','seedFile','-v7');
s=load(cacheFile,'cache');[fullModel,fullDetail]=applyCutCellSelectedDirection(s.cache,value);[seedModel,seedDetail]=applyCutCellSelectedDirection(s.cache,initial);clear s;
report=struct('status','tested-fresh-current-state-peak-weighted-total-direction-research', ...
    'sourceCase',sourceCase,'sourceBinding',sourceBinding,'sourceIdentity',sourceIdentity,'seedFile',seedFile, ...
    'seedFileSHA256',rocketStreamingFileSHA256(seedFile),'currentStateCompleteNSExactlyReplayed',true, ...
    'oldNumericMatricesFactorsOrDerivativesReused',false,'constructionSeconds',constructionSeconds,'factorReport',factorReport, ...
    'columnConstruction',columnConstruction,'weighting',weighting,'linearReport',detail, ...
    'currentSeedCompleteModelRelativeDefect',norm(rhs+seedModel,'fro')/norm(rhs,'fro'), ...
    'completeSelectedModelRelativeDefect',norm(rhs+fullModel,'fro')/norm(rhs,'fro'), ...
    'maximumRawPrimitiveComponentsBefore',max(abs(unweightedInitial),[],1), ...
    'maximumRawPrimitiveComponentsAfter',max(abs(unweightedFinal),[],1), ...
    'rawPrimitiveUpdateNormBefore',norm(unweightedInitial,'fro'),'rawPrimitiveUpdateNormAfter',norm(unweightedFinal,'fro'), ...
    'newSelectedModelApplications',2,'newActualNSEvaluations',1,'newFlowAcceptances',0, ...
    'physicalStateChanged',false,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
report.modelPointOneMet=report.completeSelectedModelRelativeDefect<=.1;
save(fullfile(folder,'current_peak_model_evidence.mat'),'value','fullModel','fullDetail','seedModel','seedDetail','report','-v7');
fid=fopen(fullfile(folder,'current_peak_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
