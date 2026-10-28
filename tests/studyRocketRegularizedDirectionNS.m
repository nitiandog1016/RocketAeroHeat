function report=studyRocketRegularizedDirectionNS(folder)
%STUDYROCKETREGULARIZEDDIRECTIONNS 用原NS及真实线搜索检验同Q20的总方向正则化分支。
% 模块一：新候选与Q20/几何/源码绑定，先另作原NS探针，不继承解析模型资格。
% 模块二：通过后从Q20实际执行原求解器；与已有Q21是并列方向对照，不拼两条历史。
% 模块三：保留全部原0.1、正性、温域、RMS下降、守恒及稳态门槛。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','updateradiusreal01','regularized_parent_model_evidence.mat');s=load(sourceFile,'value','report');
assert(s.report.modelPointOneMet&&~s.report.physicalStateChanged&&s.report.newActualNSEvaluations==0);
baseFolder=fullfile(root,'outputs','v83_validation','krylov20');referenceFolder=fullfile(root,'outputs','v83_validation','krylov21');
probeFolder=fullfile(root,'outputs','v83_validation','directnsprobereal04');
probe=createRocketNSDirectionProbe(probeFolder,fullfile(baseFolder,'wall_velocity_krylov_evidence.mat'), ...
    fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain',s.value,'result');clear s;
candidateFolder=fullfile(root,'outputs','v83_validation','krylov22');main=struct([]);
if probe.originalNSDirectionTargetMet
    options=struct('linearEquation','steady-newton','normalEulerLinearization','exact-branch','jacobianAction','original-difference', ...
        'maximumIterations',1,'initialPseudoDt',1e-7,'maximumKrylovIterations',48,'linearRelativeTolerance',.1, ...
        'maximumStepRetries',0,'maximumLineSearch',16,'steadyInitialGuess', ...
        struct('file',fullfile(probeFolder,'direct_ns_direction_probe.mat'),'variant','ns-probe','scale',1));
    main=continueRocketWallVelocityKrylov(candidateFolder,baseFolder,options);
end
reference=jsondecode(fileread(fullfile(referenceFolder,'wall_velocity_krylov_summary.json')));
report=struct('status','tested-alternate-Q20-regularized-direction-with-original-NS-and-unchanged-main-line-search', ...
    'candidateDirectionSource',sourceFile,'candidateDirectionSourceSHA256',rocketStreamingFileSHA256(sourceFile), ...
    'alternatePredecessorCase','krylov20','referenceAcceptedCase','krylov21','candidateCase','krylov22', ...
    'sameBaseDirectionComparison',true,'alternateHistoriesMerged',false,'probe',probe,'main',main, ...
    'referenceRmsMerit',reference.rmsMerit,'referenceMaximumResidual',reference.maximumResidual, ...
    'referenceAcceptedAlpha',1/4096,'newProbeActualNSEvaluations',probe.newActualNSEvaluations, ...
    'newMainActualNSEvaluations',0,'candidateRMSBetterThanReference',false,'candidateMaximumBetterThanReference',false, ...
    'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
if ~isempty(main)
    report.newMainActualNSEvaluations=main.residualEvaluations;
    report.candidateRMSBetterThanReference=main.rmsMerit<reference.rmsMerit;
    report.candidateMaximumBetterThanReference=main.maximumResidual<reference.maximumResidual;
end
save(fullfile(folder,'regularized_direction_ns_evidence.mat'),'report','probe','main','-v7');
fid=fopen(fullfile(folder,'regularized_direction_ns_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
