function report=testShockFittedNodalDynamics(folder,legacyResumeFile)
%TESTSHOCKFITTEDNODALDYNAMICS 检查原推进逐位兼容和新节点特征推进的实际续算。
% 模块一：同一v102六步原始输入，原默认数值不变，不复跑旧大算例。
% 模块二：新节点运动六步整算与三加三续算逐位核对，切换方法的检查点必须拒绝。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setupShockFittedFlow;
assert(~exist(folder,'dir'),'ShockFit:Overwrite','节点推进回归使用新目录。');mkdir(folder);
old=load(legacyResumeFile,'complete');reference=old.complete;
up=reference.upstream;gas=reference.gas;cap=reference.cap;
options=struct('maximumSteps',6,'endTimeS',1e-4);legacy=solveShockFittedCap3D(up,gas,cap,options);
assert(isequal(legacy.Q,reference.Q)&&isequal(legacy.nodes,reference.nodes)&& ...
    isequal(legacy.history,reference.history),'ShockFit:NodalLegacyReplay','原六步推进没有逐位复现。');
options.motionMethod='nodal-characteristic';complete=solveShockFittedCap3D(up,gas,cap,options);
partOptions=options;partOptions.maximumSteps=3;partOptions.checkpointFile=fullfile(folder,'first_checkpoint.mat');
part=solveShockFittedCap3D(up,gas,cap,partOptions);
continueOptions=partOptions;continueOptions.resumeFile=partOptions.checkpointFile;
continueOptions.checkpointFile=fullfile(folder,'resumed_checkpoint.mat');
resumed=solveShockFittedCap3D(up,gas,cap,continueOptions);
assert(isequal(complete.Q,resumed.Q)&&isequal(complete.nodes,resumed.nodes)&& ...
    isequal(complete.history,resumed.history)&& ...
    isequal(complete.shockTracePressureMismatchHistory,resumed.shockTracePressureMismatchHistory)&& ...
    isequal(complete.conservedIntegrals,resumed.conservedIntegrals)&& ...
    isequal(complete.boundaryFluxRates,resumed.boundaryFluxRates), ...
    'ShockFit:NodalResumeReplay','新节点方法分段续算没有逐位复现。');
different=continueOptions;different.motionMethod='face-riemann-L2';
id='';try,solveShockFittedCap3D(up,gas,cap,different);catch ex,id=ex.identifier;end
assert(strcmp(id,'ShockFit:Resume'),'ShockFit:NodalResumeMethod','不同运动方法的检查点未拒绝。');
assert(all(isfinite(complete.shockTracePressureMismatchHistory))&&~complete.productionQualified);
report=struct('status','tested-bitwise-legacy-compatibility-and-nodal-dynamics-resume', ...
    'legacyQNodesHistoryBitwise',true,'nodalQNodesHistoryBitwise',true, ...
    'nodalReconstructedPressureHistoryBitwise',true,'checkpointMethodMismatchRejected',true, ...
    'steps',6,'physicalQualification',false,'productionQualified',false,'newProductionSolvedAreaM2',0, ...
    'engineProcessID',feature('getpid'));
save(fullfile(folder,'nodal_dynamics_evidence.mat'),'report','legacy','complete','part','resumed','legacyResumeFile','-v7');
fid=fopen(fullfile(folder,'nodal_dynamics_summary.json'),'w','n','UTF-8');assert(fid>=0);
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));disp(report);
end
