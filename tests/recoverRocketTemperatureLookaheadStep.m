function report=recoverRocketTemperatureLookaheadStep(folder,failedStepFolder)
%RECOVERROCKETTEMPERATURELOOKAHEADSTEP 从已绑定同Q37自由前缀恢复目标驱动受限方向。
% 模块一：核对所有产物SHA、物理源码、几何、原行、QR坐标提升与完整原C/冻结W目标。
% 模块二：保留原10K全部行，补全部预测越域CV的原态面；选面不使用秩或NS结果。
% 模块三：与正常入口共用有界加约束循环；朝内行可以自由，不假定全部行朝外。
% 模块四：新方向依次通过完整当前模型、新鲜原NS与原线搜索；原失败和准备1NS单列保留。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
folder=char(java.io.File(folder).getCanonicalPath());failedStepFolder=char(java.io.File(failedStepFolder).getCanonicalPath());
assert(~exist(folder,'dir'));mkdir(folder);
contextFile=fullfile(failedStepFolder,'direction_context.mat');prefixFile=fullfile(failedStepFolder,'linear_checkpoints/lookahead_free_prefix.mat');
rankFile=fullfile(root,'outputs/v90_validation/rank01/temperature_lookahead_rank.mat');bound=load(rankFile,'report');
assert(strcmp(rocketStreamingFileSHA256(contextFile),bound.report.contextSHA256)&& ...
    strcmp(rocketStreamingFileSHA256(prefixFile),bound.report.freePrefixSHA256), ...
    'TemperatureRecovery:Prefix','失败上下文或自由前缀SHA与已完成来源核对不符。');
c=load(contextFile,'factor','initialRows','p','preparationFile','preparationSHA256');
assert(strcmp(rocketStreamingFileSHA256(c.preparationFile),c.preparationSHA256));
p=loadRocketPeakSystem(c.preparationFile,c.p.sourceBinding.sourceFile);assert(isequaln(p,c.p));
factor=c.factor;initialRows=c.initialRows;preparationFile=c.preparationFile;preparationSHA256=c.preparationSHA256;clear c;
f=load(prefixFile,'free','freeDetail','freeCoordinates','originalRows');assert(isequaln(f.originalRows,initialRows));
assert(strcmp(f.freeDetail.method,'complete-parent-C-fixed-primitive-penalty-augmented-patch-QR-LSQR')&&f.freeDetail.nativeIterations==128);
lift=applyRocketAugmentedPatchQR(factor,f.freeCoordinates);
assert(isequaln(lift,f.free),'TemperatureRecovery:Lift','当前QR坐标无法逐项还原保存的自由方向。');clear lift;
a=load(p.artifacts.packed.file,'packed');m=load(p.artifacts.metric.file);rhs=m.rhs;metric=m.metric;clear m;
image=applyRocketPackedGraph(a.packed,f.free);weighted=applyRocketPrimitiveUpdateMetric(metric,f.free);
relative=norm(rhs-image,'fro')/norm(rhs,'fro');objective=relative^2+(factor.penaltyInOriginalResidualUnits*norm(weighted,'fro')/norm(rhs,'fro'))^2;
assert(abs(relative-f.freeDetail.relativeCResidual)<=1e-10&&abs(objective-f.freeDetail.normalizedAugmentedObjective)<=1e-10, ...
    'TemperatureRecovery:Objective','自由前缀不属于当前完整C/冻结W目标。');clear image weighted;
c=load(p.artifacts.cache.file,'cache');currentRows=selectRocketTemperatureBandRows(c.cache.Q,c.cache.baseInfo,c.cache.geometry,c.cache.gas,c.cache.scale);
assert(isequaln(currentRows,initialRows),'TemperatureRecovery:Rows','当前原10K全部行、符号或元数据改变。');clear currentRows;
[rows,selection]=selectRocketTemperatureLookaheadRows(c.cache.Q,c.cache.baseInfo,c.cache.geometry,c.cache.gas,c.cache.scale,f.free,1/8, ...
    struct('selectionPolicy','original-and-predicted-required-faces'));clear c;
freeReplay=struct('file',prefixFile,'sha256',rocketStreamingFileSHA256(prefixFile), ...
    'contextFile',contextFile,'contextSHA256',rocketStreamingFileSHA256(contextFile), ...
    'currentPreparationSHA256',preparationSHA256,'coordinatesLiftExactlyVerified',true, ...
    'currentOriginalCResidual',relative,'currentAugmentedObjective',objective,'allCurrentRowsAndSignsExactlyVerified',true, ...
    'reusedNativeIterations',128,'newFreeNativeIterations',0,'newActualNSEvaluations',0);
checkpoints=fullfile(folder,'linear_checkpoints');mkdir(checkpoints);
options=struct('maximumIterations',128,'maximumActiveSolves',3,'checkpointFolder',checkpoints);
save(fullfile(folder,'direction_context.mat'),'factor','rows','initialRows','p','preparationFile','preparationSHA256','freeReplay','-v7');
[value,directionReport,evidence]=finishRocketTemperatureActiveDirection(a.packed,factor,rhs,metric,rows,options,f.free,f.freeDetail,f.freeCoordinates);
directionReport.lookaheadSelection=selection;directionReport.freeCheckpointReplay=freeReplay;
directionReport.reusedFreeNativeIterations=128;directionReport.newNativeIterations=directionReport.totalNativeIterations-128;
evidence.freeCheckpointReplay=freeReplay;clear a factor f;
sourceBinding=p.sourceBinding;sourceIdentity=p.sourceIdentity;geometryBinding=p.geometryBinding;
save(fullfile(folder,'active_direction.mat'),'value','directionReport','evidence','rhs','metric','rows','initialRows', ...
    'sourceBinding','sourceIdentity','geometryBinding','preparationFile','preparationSHA256','-v7');clear evidence;
c=load(p.artifacts.cache.file,'cache');[model,modelDetail]=applyCutCellSelectedDirection(c.cache,value);clear c;
defect=norm(rhs+model,'fro')/norm(rhs,'fro');stateFile=sourceBinding.sourceFile;
row=struct('step',1,'sourceBinding',sourceBinding,'preparationFile',preparationFile,'preparationSHA256',preparationSHA256, ...
    'directionFile',fullfile(folder,'active_direction.mat'),'directionFileSHA256',rocketStreamingFileSHA256(fullfile(folder,'active_direction.mat')), ...
    'directionReport',directionReport,'completeSelectedModelRelativeDefect',defect, ...
    'newPreparationActualNS',1,'newPreparationActualNSInRecovery',0,'preparationFromInitialFailedAttempt',true, ...
    'probe',[],'main',[],'resultFolder','','newProbeActualNS',0,'newMainActualNS',0,'stopReason','direction-or-model-gate', ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'active_model.mat'),'value','model','modelDetail','row','-v7');lastAccepted=fileparts(stateFile);
if directionReport.bandDirectionFeasible&&defect<=.1
    row.probe=createRocketNSDirectionProbe(fullfile(folder,'probe'),stateFile,p.geometryBinding.file,p.geometryBinding.variable,value,'result');
    row.newProbeActualNS=row.probe.newActualNSEvaluations;
    if row.probe.originalNSDirectionTargetMet
        options=struct('linearEquation','steady-newton','normalEulerLinearization','exact-branch','jacobianAction','original-difference', ...
            'maximumIterations',1,'initialPseudoDt',1e-7,'maximumKrylovIterations',48,'linearRelativeTolerance',.1, ...
            'maximumStepRetries',0,'maximumLineSearch',16,'steadyInitialGuess', ...
            struct('file',fullfile(folder,'probe','direct_ns_direction_probe.mat'),'variant','ns-probe','scale',1));
        row.resultFolder=fullfile(folder,'krylov01');row.main=continueRocketWallVelocityKrylov(row.resultFolder,fileparts(stateFile),options);
        row.newMainActualNS=row.main.residualEvaluations;row.stopReason=row.main.stopReason;
        if row.main.acceptedNonlinearIterations>0,lastAccepted=row.resultFolder;end
    else,row.stopReason='fresh-original-NS-direction-gate';end
end
write(fullfile(folder,'step.json'),row);
failureFile=fullfile(root,'outputs/v83_validation/lookaheadfull90_process.json');failure=jsondecode(fileread(failureFile));
assert(strcmp(failure.status,'terminal')&&failure.exitCode==1);
report=struct('status','tested-current-targeted-lookahead-with-bound-same-state-free-prefix-recovery', ...
    'initialStateFile',stateFile,'maximumSteps',1,'targetPreviewAlpha',1/8,'steps',{{row}}, ...
    'lastAcceptedFolder',lastAccepted,'sameEquationHistoryContinued',true,'oldActiveSetsInherited',false, ...
    'recovery',struct('originalFailureFile',failureFile,'originalFailureSHA256',rocketStreamingFileSHA256(failureFile), ...
        'originalFailureExitCode',1,'recoveryFolder',folder,'freeCheckpointReplay',freeReplay), ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
phase=fileparts(failedStepFolder);write(fullfile(phase,'progress.json'),{row});write(fullfile(phase,'active_continuation.json'),report);disp(report);
end
function write(file,data)
assert(~isfile(file));fid=fopen(file,'w','n','UTF-8');assert(fid>=0);closer=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(data));
end
