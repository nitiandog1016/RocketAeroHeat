function report=recoverRocketTemperatureRepreviewStep(stepFolder,bindingFile)
%RECOVERROCKETTEMPERATUREREPREVIEWSTEP 从同态失败前缀接续全部条件的逐轮投影，不重算自由方向。
% 模块一：绑定原失败、准备产物、上下文、自由方向及首轮预览；当前C/W/坐标逐项核对。
% 模块二：只变换约束表示，复用当前原行和自由128步；受限预算与原NS护栏保持。
% 模块三：成功或拒绝均写入真实研究记录，旧失败和来源原文不覆盖。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
stepFolder=char(java.io.File(stepFolder).getCanonicalPath());
binding=jsondecode(fileread(bindingFile));assert(strcmp(binding.kind,'same-Q40-free-prefix-and-first-preview-recovery-binding'));
for key={'failure','context','prefix','preview','independentProjectorProof'}
    entry=binding.(key{1});path=fullfile(root,entry.file);
    assert(strcmp(rocketStreamingFileSHA256(path),entry.sha256),'TemperatureRecovery:Binding','原失败或检查点SHA改变。');
end
failureFile=fullfile(root,binding.failure.file);failure=jsondecode(fileread(failureFile));
assert(strcmp(failure.status,'terminal')&&failure.exitCode==1&&all(structfun(@(x)logical(x),failure.sourceStableAfter)));
% 自由求解、增广映射和原C/W代码不能因恢复而改变；允许本轮已验的投影表示升级。
for name={'solveRocketAugmentedQRDirection','applyRocketAugmentedPatchQR','applyRocketPackedGraph', ...
        'applyRocketPackedGraphTranspose','applyRocketPrimitiveUpdateMetric','buildRocketAugmentedPatchQR'}
    rel=['internalFlow/cutCell/' name{1} '.m'];key=matlab.lang.makeValidName(rel);
    assert(isfield(failure.sourceSHA256,key)&&strcmp(rocketStreamingFileSHA256(fullfile(root,rel)),failure.sourceSHA256.(key)), ...
        'TemperatureRecovery:Source','当前自由目标或右映射与失败运行源码不同。');
end
contextFile=fullfile(root,binding.context.file);prefixFile=fullfile(root,binding.prefix.file);previewFile=fullfile(root,binding.preview.file);
c=load(contextFile,'factor','initialRows','p','preparationFile','preparationSHA256');
assert(strcmp(rocketStreamingFileSHA256(c.preparationFile),c.preparationSHA256));
p=loadRocketPeakSystem(c.preparationFile,c.p.sourceBinding.sourceFile);assert(isequaln(p,c.p));
assert(~isfile(fullfile(stepFolder,'active_direction.mat'))&&~isfile(fullfile(stepFolder,'step.json')));
f=load(prefixFile,'free','freeCoordinates','freeDetail','originalRows');s=load(previewFile,'rows','selection','summary','value');
assert(isequaln(f.originalRows,c.initialRows)&&isequaln(f.free,s.value)&&f.freeDetail.nativeIterations==128);
lift=applyRocketAugmentedPatchQR(c.factor,f.freeCoordinates);assert(isequaln(lift,f.free));clear lift;
a=load(p.artifacts.packed.file,'packed');m=load(p.artifacts.metric.file);rhs=m.rhs;metric=m.metric;clear m;
image=applyRocketPackedGraph(a.packed,f.free);weighted=applyRocketPrimitiveUpdateMetric(metric,f.free);
relative=norm(rhs-image,'fro')/norm(rhs,'fro');objective=relative^2+ ...
    (c.factor.penaltyInOriginalResidualUnits*norm(weighted,'fro')/norm(rhs,'fro'))^2;
assert(abs(relative-f.freeDetail.relativeCResidual)<=1e-10&&abs(objective-f.freeDetail.normalizedAugmentedObjective)<=1e-10);
clear image weighted;
cache=load(p.artifacts.cache.file,'cache');cc=cache.cache;
currentRows=selectRocketTemperatureBandRows(cc.Q,cc.baseInfo,cc.geometry,cc.gas,cc.scale);
assert(isequaln(currentRows,c.initialRows));
assert(strcmp(s.selection.trialStateCoordinates,binding.coordinates)&&s.selection.targetAlpha==binding.targetAlpha);
preview=struct('Q',cc.Q,'info',cc.baseInfo,'cut',cc.geometry,'gas',cc.gas,'scale',cc.scale, ...
    'targetAlpha',binding.targetAlpha,'selectionPolicy','original-and-iteratively-predicted-cell-and-interface-faces', ...
    'trialStateCoordinates',binding.coordinates);clear cache cc currentRows;
freeReplay=struct('file',prefixFile,'sha256',binding.prefix.sha256,'contextFile',contextFile, ...
    'contextSHA256',binding.context.sha256,'currentPreparationSHA256',c.preparationSHA256, ...
    'coordinatesLiftExactlyVerified',true,'currentOriginalCResidual',relative,'currentAugmentedObjective',objective, ...
    'allCurrentRowsAndSignsExactlyVerified',true,'reusedNativeIterations',128,'newFreeNativeIterations',0,'newActualNSEvaluations',0);
checkpoints=fullfile(stepFolder,'recovered_linear_checkpoints');assert(~exist(checkpoints,'dir'));mkdir(checkpoints);
options=struct('maximumIterations',128,'maximumActiveSolves',3,'checkpointFolder',checkpoints);
[value,directionReport,evidence]=finishRocketTemperatureRepreviewDirection(a.packed,c.factor,rhs,metric,s.rows,options, ...
    f.free,f.freeDetail,f.freeCoordinates,preview,s.selection,c.initialRows);
directionReport.freeCheckpointReplay=freeReplay;directionReport.reusedFreeNativeIterations=128;
directionReport.newNativeIterations=directionReport.totalNativeIterations-128;evidence.freeCheckpointReplay=freeReplay;
rows=evidence.rows;initialRows=c.initialRows;preparationFile=c.preparationFile;preparationSHA256=c.preparationSHA256;
sourceBinding=p.sourceBinding;sourceIdentity=p.sourceIdentity;geometryBinding=p.geometryBinding;clear a c f s preview;
save(fullfile(stepFolder,'active_direction.mat'),'value','directionReport','evidence','rows','initialRows','rhs','metric', ...
    'sourceBinding','sourceIdentity','geometryBinding','preparationFile','preparationSHA256','-v7');clear evidence;
cache=load(p.artifacts.cache.file,'cache');[model,modelDetail]=applyCutCellSelectedDirection(cache.cache,value);clear cache;
defect=norm(rhs+model,'fro')/norm(rhs,'fro');stateFile=sourceBinding.sourceFile;
row=struct('step',1,'sourceBinding',sourceBinding,'preparationFile',preparationFile,'preparationSHA256',preparationSHA256, ...
    'directionFile',fullfile(stepFolder,'active_direction.mat'),'directionFileSHA256',rocketStreamingFileSHA256(fullfile(stepFolder,'active_direction.mat')), ...
    'directionReport',directionReport,'completeSelectedModelRelativeDefect',defect,'newPreparationActualNS',1, ...
    'newPreparationActualNSInRecovery',0,'preparationFromInitialFailedAttempt',true, ...
    'probe',[],'main',[],'resultFolder','','newProbeActualNS',0,'newMainActualNS',0,'stopReason','direction-or-model-gate', ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(stepFolder,'active_model.mat'),'value','model','modelDetail','row','-v7');lastAccepted=fileparts(stateFile);
if directionReport.bandDirectionFeasible&&defect<=.1
    row.probe=createRocketNSDirectionProbe(fullfile(stepFolder,'probe'),stateFile,p.geometryBinding.file,p.geometryBinding.variable,value,'result');
    row.newProbeActualNS=row.probe.newActualNSEvaluations;
    if row.probe.originalNSDirectionTargetMet
        opts=struct('linearEquation','steady-newton','normalEulerLinearization','exact-branch','jacobianAction','original-difference', ...
            'maximumIterations',1,'initialPseudoDt',1e-7,'maximumKrylovIterations',48,'linearRelativeTolerance',.1, ...
            'maximumStepRetries',0,'maximumLineSearch',16,'trialStateCoordinates',binding.coordinates,'steadyInitialGuess', ...
            struct('file',fullfile(stepFolder,'probe','direct_ns_direction_probe.mat'),'variant','ns-probe','scale',1));
        row.resultFolder=fullfile(stepFolder,'krylov01');row.main=continueRocketWallVelocityKrylov(row.resultFolder,fileparts(stateFile),opts);
        row.newMainActualNS=row.main.residualEvaluations;row.stopReason=row.main.stopReason;
        if row.main.acceptedNonlinearIterations>0,lastAccepted=row.resultFolder;end
    else,row.stopReason='fresh-original-NS-direction-gate';end
end
write(fullfile(stepFolder,'step.json'),row);phase=fileparts(stepFolder);write(fullfile(phase,'progress.json'),{row});
report=struct('status','tested-bound-same-state-repreview-with-exact-condition-row-representation-recovery', ...
    'initialStateFile',stateFile,'maximumSteps',1,'targetPreviewAlpha',binding.targetAlpha,'trialStateCoordinates',binding.coordinates, ...
    'steps',{{row}},'lastAcceptedFolder',lastAccepted,'sameEquationHistoryContinued',true,'oldActiveSetsInherited',false, ...
    'recovery',struct('failureFile',failureFile,'failureSHA256',binding.failure.sha256,'freeCheckpointReplay',freeReplay), ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
write(fullfile(phase,'active_continuation.json'),report);disp(report);
end
function write(file,data)
assert(~isfile(file));fid=fopen(file,'w','n','UTF-8');assert(fid>=0);closer=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(data));
end
