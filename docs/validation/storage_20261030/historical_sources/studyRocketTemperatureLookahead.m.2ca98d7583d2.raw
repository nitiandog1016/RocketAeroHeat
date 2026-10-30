function report=studyRocketTemperatureLookahead(folder,previousFolder,targetAlpha,selectionPolicy)
%STUDYROCKETTEMPERATURELOOKAHEAD 从显式真实端点重建一次前瞻选面方向并由原NS验收。
% 模块一：完整回放当前NS、重建C/W/QR；旧Q35夹具和Q36导数不用于新状态。
% 模块二：自由128步只算一次；显式逐轮策略最多三次受限重算，检查投影后新暴露的温域条件。
% 模块三：方向仍过原完整模型/新鲜NS0.1，候选仍从全步按原温域/下降规则搜索。
% 模块四：仅一次有界真实推进，保留自由前缀/全部原行/拒绝；不保证目标步长可接受。
if nargin<2,previousFolder='outputs/v83_validation/krylov37';end
if nargin<3,targetAlpha=1/8;end
if nargin<4,selectionPolicy='all-extended-band';end
assert(ischar(selectionPolicy)&&ismember(selectionPolicy,{'all-extended-band','original-and-predicted-required-faces', ...
    'original-and-predicted-cell-and-interface-faces','original-and-iteratively-predicted-cell-and-interface-faces'}), ...
    'TemperatureLookahead:Options','需要明确的前瞻选面策略。');
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
folder=char(java.io.File(folder).getCanonicalPath());previousFolder=char(java.io.File(previousFolder).getCanonicalPath());
stateFile=fullfile(previousFolder,'wall_velocity_krylov_evidence.mat');assert(isfile(stateFile)&&~exist(folder,'dir'));mkdir(folder);
geometry=fullfile(root,'outputs/v80_validation/domain02/shell_domain_evidence.mat');
stepFolder=fullfile(folder,'step01');mkdir(stepFolder);systemFolder=fullfile(stepFolder,'current');
prepareRocketCurrentPeakSystem(systemFolder,stateFile,geometry,'domain');
preparationFile=fullfile(systemFolder,'peak_system.json');p=loadRocketPeakSystem(preparationFile,stateFile);
preparationSHA256=rocketStreamingFileSHA256(preparationFile);
a=load(p.artifacts.packed.file,'packed');b=load(p.artifacts.factor.file,'factor');m=load(p.artifacts.metric.file);
factor=buildRocketAugmentedPatchQR(a.packed,b.factor,m.metric,.1*norm(m.rhs,'fro')/p.metricNormalization);clear b;
c=load(p.artifacts.cache.file,'cache');initialRows=selectRocketTemperatureBandRows(c.cache.Q,c.cache.baseInfo,c.cache.geometry,c.cache.gas,c.cache.scale);
preview=struct('Q',c.cache.Q,'info',c.cache.baseInfo,'cut',c.cache.geometry,'gas',c.cache.gas,'scale',c.cache.scale,'targetAlpha',targetAlpha, ...
    'selectionPolicy',selectionPolicy);clear c;
save(fullfile(stepFolder,'direction_context.mat'),'factor','initialRows','preparationFile','preparationSHA256','p','-v7');
checkpoints=fullfile(stepFolder,'linear_checkpoints');mkdir(checkpoints);
[value,directionReport,evidence]=solveRocketTemperatureActiveDirection(a.packed,factor,m.rhs,m.metric,initialRows, ...
    struct('checkpointFolder',checkpoints,'lookaheadPreview',preview));
rows=evidence.rows;rhs=m.rhs;metric=m.metric;clear preview a m factor;
sourceBinding=p.sourceBinding;sourceIdentity=p.sourceIdentity;geometryBinding=p.geometryBinding;
save(fullfile(stepFolder,'active_direction.mat'),'value','directionReport','evidence','rows','initialRows','rhs','metric', ...
    'sourceBinding','sourceIdentity','geometryBinding','preparationFile','preparationSHA256','-v7');clear evidence;
c=load(p.artifacts.cache.file,'cache');[model,modelDetail]=applyCutCellSelectedDirection(c.cache,value);clear c;
defect=norm(rhs+model,'fro')/norm(rhs,'fro');
row=struct('step',1,'sourceBinding',sourceBinding,'preparationFile',preparationFile,'preparationSHA256',preparationSHA256, ...
    'directionFile',fullfile(stepFolder,'active_direction.mat'),'directionFileSHA256',rocketStreamingFileSHA256(fullfile(stepFolder,'active_direction.mat')), ...
    'directionReport',directionReport,'completeSelectedModelRelativeDefect',defect,'newPreparationActualNS',1, ...
    'probe',[],'main',[],'resultFolder','','newProbeActualNS',0,'newMainActualNS',0, ...
    'stopReason','direction-or-model-gate','physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(stepFolder,'active_model.mat'),'value','model','modelDetail','row','-v7');
lastAccepted=previousFolder;
if directionReport.bandDirectionFeasible&&defect<=.1
    probeFolder=fullfile(stepFolder,'probe');probe=createRocketNSDirectionProbe(probeFolder,stateFile,geometry,'domain',value,'result');
    row.probe=probe;row.newProbeActualNS=probe.newActualNSEvaluations;
    if probe.originalNSDirectionTargetMet
        options=struct('linearEquation','steady-newton','normalEulerLinearization','exact-branch','jacobianAction','original-difference', ...
            'maximumIterations',1,'initialPseudoDt',1e-7,'maximumKrylovIterations',48,'linearRelativeTolerance',.1, ...
            'maximumStepRetries',0,'maximumLineSearch',16,'steadyInitialGuess', ...
            struct('file',fullfile(probeFolder,'direct_ns_direction_probe.mat'),'variant','ns-probe','scale',1));
        row.resultFolder=fullfile(stepFolder,'krylov01');row.main=continueRocketWallVelocityKrylov(row.resultFolder,previousFolder,options);
        row.newMainActualNS=row.main.residualEvaluations;row.stopReason=row.main.stopReason;
        if row.main.acceptedNonlinearIterations>0,lastAccepted=row.resultFolder;end
    else,row.stopReason='fresh-original-NS-direction-gate';end
end
write(fullfile(stepFolder,'step.json'),row);write(fullfile(folder,'progress.json'),{row});
report=struct('status','tested-bounded-current-state-predictive-temperature-selection-NS-study','initialStateFile',stateFile, ...
    'maximumSteps',1,'targetPreviewAlpha',targetAlpha,'steps',{{row}},'lastAcceptedFolder',lastAccepted, ...
    'sameEquationHistoryContinued',true,'oldActiveSetsInherited',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
write(fullfile(folder,'active_continuation.json'),report);disp(report);
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);closer=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(data));
end
