function report=testSU2InitializationWorkbench(out)
%TESTSU2INITIALIZATIONWORKBENCH 复用原混合网格，仅检查主界面的初猜/保存上限及拒绝入口。
% 不重建网格、不启动求解器；完整初猜本身由本轮实际准备和独立守恒审核验证。
assert(~isfolder(out));mkdir(out);b=load('outputs/v118_validation/actual01/surface.mat');
c=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=c.cfg;
cfg.output.verbose=false;cfg.su2=struct('volumeMode','hybridTriangles','iterations',300,'outputEvery',150, ...
    'namedCheckpoints',true,'maximumSeconds',720,'pythonExecutable','D:/anacondapython/python.exe');
fig=rocketAeroHeatWorkbench(cfg);g=onCleanup(@()close(fig));
fig.UserData.actions.useSurfaceForSU2(b.s);fig.UserData.actions.readHybridMesh('outputs/v118_validation/actual01/mesh');
fig.UserData.actions.selectStep(3);controls=fig.UserData.controls;
assert(strcmp(controls.initialMode.Value,'freestream') && controls.checkpointLimit.Value==100);
controls.checkpointLimit.Value=256;controls.checkpointLimit.ValueChangedFcn(controls.checkpointLimit,[]);
controls=fig.UserData.controls;controls.initialThickness.Value=11.829243601117;
controls.initialThickness.ValueChangedFcn(controls.initialThickness,[]);
controls=fig.UserData.controls;controls.initialMode.Value='wallEnthalpy';controls.initialMode.ValueChangedFcn(controls.initialMode,[]);
fig.UserData.actions.confirmBoundary();state=fig.UserData.workflowState;
assert(strcmp(state.cfg.su2.initialMode,'wallEnthalpy') && state.cfg.su2.checkpointLimitMiB==256);
assert(abs(state.cfg.su2.initialThicknessM-.011829243601117)<1e-14 && isempty(fig.UserData.engineJob));
fig.UserData.actions.selectStep(3);controls=fig.UserData.controls;
controls.initialThickness.Value=100;controls.initialThickness.ValueChangedFcn(controls.initialThickness,[]);
assert(contains(fig.UserData.status.Text,'超出完整棱柱区') && ~fig.UserData.workflowState.busy);
assert(abs(fig.UserData.workflowState.cfg.su2.initialThicknessM-.011829243601117)<1e-14);
controls.initialThickness.Value=11.829243601117;controls.thermalMode.Value='adiabatic';controls.thermalMode.ValueChangedFcn(controls.thermalMode,[]);
assert(contains(fig.UserData.status.Text,'只支持给定统一壁温'));
% 默认限制行为保留，显式扩展只改资源上限，不改变物理状态或残差门槛。
v=struct('vertices',zeros(210914,3),'wall',{{1:9878}});
old=estimateSU2CheckpointStorage(v,300,150,true);new=estimateSU2CheckpointStorage(v,300,150,true,256);
assert(old.limitBytes==100*1024^2 && ~old.withinLimit && new.withinLimit && old.estimatedBytes==new.estimatedBytes);
report=struct('mainBoundaryControlsTested',true,'initialModeConfirmed',true,'tooThickRejected',true, ...
    'unsupportedThermalModeRejected',true,'default100MiBPreserved',true,'explicit256MiBAccepted',true, ...
    'mainSolveLaunchRetested',false,'solverLaunched',false,'mesherLaunched',false,'physicsQualified',false);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);f=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Main initialization controls and explicit checkpoint budget checked; no solver or mesher repeated.\n');
end
