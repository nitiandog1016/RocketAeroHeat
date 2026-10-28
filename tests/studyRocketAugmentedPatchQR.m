function report=studyRocketAugmentedPatchQR(folder,preparationFile)
%STUDYROCKETAUGMENTEDPATCHQR 从真实当前检查点构造QR方向，再按原NS门槛有界试算一步。
% 模块一：准备系统/源状态/几何全部重查，固定同一峰值权重及penalty，仅改变右预处理。
% 模块二：128步LSQR后单列原C和完整模型；通过完整模型才进行新鲜原NS探针。
% 模块三：探针通过后由原主求解器试算；模型检查、探针和主NS调用分别记账。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
folder=char(java.io.File(folder).getCanonicalPath());assert(~exist(folder,'dir'));mkdir(folder);
p=loadRocketPeakSystem(preparationFile);preparationSHA256=rocketStreamingFileSHA256(preparationFile);
m=load(p.artifacts.metric.file);a=load(p.artifacts.packed.file,'packed');b=load(p.artifacts.factor.file,'factor');
penalty=.1*norm(m.rhs,'fro')/p.metricNormalization;
factor=buildRocketAugmentedPatchQR(a.packed,b.factor,m.metric,penalty);clear b;
save(fullfile(folder,'augmented_qr_factor.mat'),'factor','preparationFile','preparationSHA256','-v7');
[value,detail,coordinates]=solveRocketAugmentedQRDirection(a.packed,factor,m.rhs,m.metric);clear a;
construction=factor.constructionReport;clear factor;
rhs=m.rhs;metric=m.metric;baseMetric=m.baseMetric;raw=applyRocketPrimitiveUpdateMetric(baseMetric,value);clear m;
sourceBinding=p.sourceBinding;sourceIdentity=p.sourceIdentity;geometryBinding=p.geometryBinding;
save(fullfile(folder,'augmented_qr_direction.mat'),'value','detail','coordinates','rhs','metric','raw', ...
    'preparationFile','preparationSHA256','sourceBinding','sourceIdentity','geometryBinding','construction','-v7');
c=load(p.artifacts.cache.file,'cache');[fullModel,fullDetail]=applyCutCellSelectedDirection(c.cache,value);clear c;
modelDefect=norm(rhs+fullModel,'fro')/norm(rhs,'fro');
report=struct('status','tested-current-state-augmented-patch-QR-direction-and-gated-one-step-NS-study', ...
    'sourceBinding',sourceBinding,'preparationFile',preparationFile,'preparationSHA256',preparationSHA256, ...
    'construction',construction,'linearReport',detail,'completeSelectedModelRelativeDefect',modelDefect, ...
    'modelPointOneMet',modelDefect<=.1,'probe',[],'main',[],'newProbeActualNS',0,'newMainActualNS',0, ...
    'newSelectedModelApplications',1,'newProductionSolvedAreaM2',0,'physicalSolution',false);
save(fullfile(folder,'augmented_qr_model.mat'),'value','fullModel','fullDetail','report','-v7');
if modelDefect<=.1
    probeFolder=fullfile(folder,'probe');
    probe=createRocketNSDirectionProbe(probeFolder,p.sourceBinding.sourceFile,p.geometryBinding.file,p.geometryBinding.variable,value,'result');
    report.probe=probe;report.newProbeActualNS=probe.newActualNSEvaluations;
    if probe.originalNSDirectionTargetMet
        options=struct('linearEquation','steady-newton','normalEulerLinearization','exact-branch','jacobianAction','original-difference', ...
            'maximumIterations',1,'initialPseudoDt',1e-7,'maximumKrylovIterations',48,'linearRelativeTolerance',.1, ...
            'maximumStepRetries',0,'maximumLineSearch',16,'steadyInitialGuess', ...
            struct('file',fullfile(probeFolder,'direct_ns_direction_probe.mat'),'variant','ns-probe','scale',1));
        report.main=continueRocketWallVelocityKrylov(fullfile(folder,'krylov01'),fileparts(p.sourceBinding.sourceFile),options);
        report.newMainActualNS=report.main.residualEvaluations;
    end
end
fid=fopen(fullfile(folder,'augmented_qr_study.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));disp(report);
end
