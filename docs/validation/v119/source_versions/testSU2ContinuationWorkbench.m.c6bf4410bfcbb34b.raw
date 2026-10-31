function report=testSU2ContinuationWorkbench(out)
%TESTSU2CONTINUATIONWORKBENCH 检查源场入口及MATLAB适配器，不重复启动引擎。
% 模块一：默认自由流保持；缺少源目录拒绝；有效源目录进入集中配置。
% 模块二：实际调用Python准备器，确认精确原生场及累计身份。
assert(~isfolder(out));mkdir(out);d=load('outputs/v118_validation/actual01/surface.mat');
c=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=c.cfg;
cfg.output.verbose=false;cfg.su2=struct('volumeMode','hybridTriangles','iterations',300,'outputEvery',150, ...
    'namedCheckpoints',true,'maximumSeconds',720,'checkpointLimitMiB',256,'pythonExecutable','D:/anacondapython/python.exe');
fig=rocketAeroHeatWorkbench(cfg);g=onCleanup(@()close(fig));
fig.UserData.actions.useSurfaceForSU2(d.s);fig.UserData.actions.readHybridMesh('outputs/v118_validation/actual01/mesh');
fig.UserData.actions.selectStep(3);controls=fig.UserData.controls;
assert(strcmp(controls.initialMode.Value,'freestream'));
controls.initialMode.Value='checkpoint';controls.initialMode.ValueChangedFcn(controls.initialMode,[]);
assert(contains(fig.UserData.status.Text,'请选择'));
source=fullfile(pwd,'outputs','v118_validation','actual01','budget_run');
controls.checkpointSource.Value=source;controls.checkpointSource.ValueChangedFcn(controls.checkpointSource,[]);
fig.UserData.actions.confirmBoundary();state=fig.UserData.workflowState;
assert(strcmp(state.cfg.su2.initialMode,'checkpoint') && strcmp(state.cfg.su2.checkpointSource,source));
assert(isempty(fig.UserData.engineJob));
fig.UserData.actions.selectStep(3);controls=fig.UserData.controls;
controls.checkpointMode.Value='overwrite';controls.checkpointMode.ValueChangedFcn(controls.checkpointMode,[]);
assert(contains(fig.UserData.status.Text,'具名检查点'));
target=fullfile(pwd,'outputs','v119_validation','matlab_adapter01');
lineage=writeSU2CheckpointInitialState(target,source,'D:/anacondapython/python.exe');
assert(lineage.completedUpdatesBeforeRun==300 && lineage.sourceSavedLocalIteration==299 && isfile(fullfile(target,'initial_flow.csv')));
report=struct('defaultFreestreamPreserved',true,'missingSourceRejected',true,'indexedOutputRequired',true, ...
    'mainBoundarySourceConfirmed',true,'matlabAdapterExecuted',true,'lineage',lineage, ...
    'solverLaunched',false,'mainEndToEndLaunchClaimed',false,'physicsQualified',false);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);g2=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Continuation controls, source refusal and actual MATLAB adapter passed; no engine run.\n');
end
