function report=continueRocketTemperatureActive(folder,previousFolder,maximumSteps)
%CONTINUEROCKETTEMPERATUREACTIVE 从真实端点逐态重建，最多三次连续活跃温度方向/原NS推进。
% 模块一：每一步都完整回放当前NS并重建C/W/QR/温度行，旧活跃集和方向资格不继承。
% 模块二：自由128步、最多三次各128步投影有界求方向；所有原物理/线性/稳态门槛保持。
% 模块三：逐步保存上下文、自由/活跃循环及真正接受结果；失败立即停，不把并列历史合并。
if nargin<3,maximumSteps=3;end
validateattributes(maximumSteps,{'numeric'},{'scalar','integer','positive','<=',3});
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
folder=char(java.io.File(folder).getCanonicalPath());previousFolder=char(java.io.File(previousFolder).getCanonicalPath());
assert(~exist(folder,'dir'));mkdir(folder);steps={};initialStateFile=fullfile(previousFolder,'wall_velocity_krylov_evidence.mat');
geometry=fullfile(root,'outputs/v80_validation/domain02/shell_domain_evidence.mat');
for k=1:maximumSteps
    stepFolder=fullfile(folder,sprintf('step%02d',k));mkdir(stepFolder);stateFile=fullfile(previousFolder,'wall_velocity_krylov_evidence.mat');
    systemFolder=fullfile(stepFolder,'current');prepareRocketCurrentPeakSystem(systemFolder,stateFile,geometry,'domain');
    preparationFile=fullfile(systemFolder,'peak_system.json');p=loadRocketPeakSystem(preparationFile,stateFile);preparationSHA256=rocketStreamingFileSHA256(preparationFile);
    a=load(p.artifacts.packed.file,'packed');b=load(p.artifacts.factor.file,'factor');m=load(p.artifacts.metric.file);
    factor=buildRocketAugmentedPatchQR(a.packed,b.factor,m.metric,.1*norm(m.rhs,'fro')/p.metricNormalization);clear b;
    c=load(p.artifacts.cache.file,'cache');rows=selectRocketTemperatureBandRows(c.cache.Q,c.cache.baseInfo,c.cache.geometry,c.cache.gas,c.cache.scale);clear c;
    save(fullfile(stepFolder,'direction_context.mat'),'factor','rows','preparationFile','preparationSHA256','p','-v7');
    checkpoints=fullfile(stepFolder,'linear_checkpoints');mkdir(checkpoints);
    [value,directionReport,evidence]=solveRocketTemperatureActiveDirection(a.packed,factor,m.rhs,m.metric,rows,struct('checkpointFolder',checkpoints));
    rhs=m.rhs;metric=m.metric;clear a factor m;
    sourceBinding=p.sourceBinding;sourceIdentity=p.sourceIdentity;geometryBinding=p.geometryBinding;
    save(fullfile(stepFolder,'active_direction.mat'),'value','directionReport','evidence','rhs','metric','rows', ...
        'sourceBinding','sourceIdentity','geometryBinding','preparationFile','preparationSHA256','-v7');clear evidence;
    c=load(p.artifacts.cache.file,'cache');[model,modelDetail]=applyCutCellSelectedDirection(c.cache,value);clear c;
    defect=norm(rhs+model,'fro')/norm(rhs,'fro');
    row=struct('step',k,'sourceBinding',sourceBinding,'preparationFile',preparationFile,'preparationSHA256',preparationSHA256, ...
        'directionFile',fullfile(stepFolder,'active_direction.mat'),'directionFileSHA256',rocketStreamingFileSHA256(fullfile(stepFolder,'active_direction.mat')), ...
        'directionReport',directionReport,'completeSelectedModelRelativeDefect',defect,'newPreparationActualNS',1, ...
        'probe',[],'main',[],'resultFolder','','newProbeActualNS',0,'newMainActualNS',0,'stopReason','direction-or-model-gate', ...
        'physicalSolution',false,'newProductionSolvedAreaM2',0);
    save(fullfile(stepFolder,'active_model.mat'),'value','model','modelDetail','row','-v7');
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
        else,row.stopReason='fresh-original-NS-direction-gate';end
    end
    fid=fopen(fullfile(stepFolder,'step.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(row));fclose(fid);steps{end+1}=row; %#ok<AGROW>
    fid=fopen(fullfile(folder,'progress.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(steps));fclose(fid);
    if isempty(row.main)||row.main.acceptedNonlinearIterations==0,break;end
    previousFolder=row.resultFolder;
end
report=struct('status','tested-bounded-current-state-active-temperature-continuous-NS-study','initialStateFile',initialStateFile, ...
    'maximumSteps',maximumSteps,'steps',{steps},'lastAcceptedFolder',previousFolder,'sameEquationHistoryContinued',true, ...
    'oldActiveSetsInherited',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
fid=fopen(fullfile(folder,'active_continuation.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(report));fclose(fid);disp(report);
end
