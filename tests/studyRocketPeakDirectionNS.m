function report=studyRocketPeakDirectionNS(folder)
%STUDYROCKETPEAKDIRECTIONNS 将当前Q22已过模型门槛的零初猜续算方向交给原NS验收。
% 模块一：只读取明确的当前状态/源文件；模型通过不能替代新鲜原NS方向检查。
% 模块二：通过原0.1后从真实Q22接续主求解，状态/温域/下降/守恒/稳态门槛保持。
% 模块三：探针、模型研究和主求解调用分列，不能拼接Q20的并列Q23历史。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);candidateFile=fullfile(root,'outputs','v83_validation','peakzeroreal02','peak_zero_model_evidence.mat');
s=load(candidateFile,'value','report');assert(s.report.modelPointOneMet&&~s.report.physicalStateChanged&&s.report.newActualNSEvaluations==0);
assert(isequaln(s.report.sourceIdentity,rocketSelectedDirectionSourceIdentity)&& ...
    strcmp(rocketStreamingFileSHA256(s.report.sourceBinding.sourceFile),s.report.sourceBinding.sourceFileSHA256));
probeFolder=fullfile(root,'outputs','v83_validation','directnsprobereal05');
probe=createRocketNSDirectionProbe(probeFolder,s.report.sourceBinding.sourceFile, ...
    fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain',s.value,'result');clear s;
main=struct([]);
if probe.originalNSDirectionTargetMet
    options=struct('linearEquation','steady-newton','normalEulerLinearization','exact-branch','jacobianAction','original-difference', ...
        'maximumIterations',1,'initialPseudoDt',1e-7,'maximumKrylovIterations',48,'linearRelativeTolerance',.1, ...
        'maximumStepRetries',0,'maximumLineSearch',16,'steadyInitialGuess',struct('file',fullfile(probeFolder,'direct_ns_direction_probe.mat'),'variant','ns-probe','scale',1));
    main=continueRocketWallVelocityKrylov(fullfile(root,'outputs','v83_validation','krylov24'), ...
        fullfile(root,'outputs','v83_validation','krylov22'),options);
end
report=struct('status','tested-current-Q22-peak-weighted-zero-continuation-original-NS-direction-and-main-gates', ...
    'candidateFile',candidateFile,'candidateFileSHA256',rocketStreamingFileSHA256(candidateFile), ...
    'predecessorCase','krylov22','candidateCase','krylov24','probe',probe,'main',main, ...
    'parallelQ23HistoryAppended',false,'newProbeActualNS',probe.newActualNSEvaluations,'newMainActualNS',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
if ~isempty(main),report.newMainActualNS=main.residualEvaluations;end
save(fullfile(folder,'peak_direction_ns_evidence.mat'),'report','probe','main','-v7');
fid=fopen(fullfile(folder,'peak_direction_ns_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
