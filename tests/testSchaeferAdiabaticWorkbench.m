function report=testSchaeferAdiabaticWorkbench(out,capture)
%TESTSCHAEFERADIABATICWORKBENCH 从主五步骤实际生成经典体网格并执行有界绝热研究。
% 模块一：显式实验来流/绝热边界随输入冻结，未指定Tw不进入给定温度场。
% 模块二：实际进程终态后保留成功/超时身份；未收敛场仅作研究。
% 模块三：检查绝热零值在云图/截面/导出标为边界输入，不是热流预测。
if nargin<2,capture=false;end
assert(~isfolder(out));mkdir(out);evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.su2=struct('pythonExecutable','D:/anacondapython/python.exe','iterations',3000, ...
    'maximumSeconds',720,'outputEvery',1000,'namedCheckpoints',true);
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()close(fig));fig.UserData.actions.loadClassic();
state=fig.UserData.workflowState;assert(strcmp(state.modelSource.kind,'schaefer-ferguson-1961'));
fig.UserData.actions.generateMesh();state=fig.UserData.workflowState;
assert(~isempty(state.mesh),fig.UserData.status.Text);fig.UserData.actions.selectStep(3);
assert(strcmp(fig.UserData.controls.thermalMode.Value,'adiabatic') && strcmp(fig.UserData.controls.wall.Enable,'off'));
assert(strcmp(fig.UserData.controls.altitude.Enable,'off'));
fig.UserData.actions.confirmBoundary();assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
assert(abs(fig.UserData.workflowState.freestream.temperature-state.modelSource.condition.atmosphere.temperature)<1e-13);
fig.UserData.actions.run();assert(fig.UserData.workflowState.busy,fig.UserData.status.Text);
write(fullfile(out,'active.json'),struct('engineFolder',fig.UserData.engineJob.folder));
t=tic;while fig.UserData.workflowState.busy && toc(t)<780,drawnow;pause(.1);end
state=fig.UserData.workflowState;assert(~state.busy,'实际任务仍运行，不能结束审查。');
folder=fig.UserData.engineJob.folder;receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
assert(strcmp(receipt.status,'terminal') && ~receipt.cancelled);
% 超时场只在身份可证明时通过既有显式入口重开；无证据时如实保留终态，禁止伪造迭代。
if receipt.timedOut,fig.UserData.actions.readTimeoutCheckpoint(folder);state=fig.UserData.workflowState;end
assert(~isempty(state.researchResult),fig.UserData.status.Text);r=state.researchResult;
assert(r.heatFluxIsBoundaryInput && all(r.heatFluxRaw==0) && all(isnan(r.wallTemperatureInput)));
ui=fig.UserData.researchUI;ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);
ui.setSectionX(state.modelSource.truth.flareLength);waitRender(ui);
assert(ui.geometry.html.UserData.payload.sectionView.valueStatus==2);
ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'section.csv'));ui.inspectFace(10);
report=struct('engineFolder',folder,'receipt',receipt,'geometry',state.modelSource.truth, ...
    'condition',state.modelSource.condition,'volumeInfo',fig.UserData.volumeInfo, ...
    'convergence',r.convergenceAudit,'rendered',ui.geometry.html.UserData.rendered, ...
    'allHeatFluxIsBoundaryInput',true,'physicsQualified',false,'newQualifiedArea',0);
save(fullfile(out,'classic_result.mat'),'r','report','-v7');write(fullfile(out,'report.json'),report);
if capture
    write(fullfile(out,'ready_for_capture.json'),struct('name',fig.Name));fprintf('V114_CLASSIC_UI_READY\n');
    t=tic;while ~isfile(fullfile(out,'capture_done.json')) && toc(t)<60,drawnow;pause(.1);end
end
ui.select.Value='wallInput';ui.select.ValueChangedFcn(ui.select,[]);waitRender(ui);
assert(all(isnan(ui.geometry.html.UserData.payload.faceField.values)) && ui.geometry.html.UserData.payload.sectionView.valueStatus==0);
fig.UserData.actions.selectStep(3);ctl=fig.UserData.controls.thermalMode;ctl.Value='isothermal';ctl.ValueChangedFcn(ctl,[]);
assert(isempty(fig.UserData.workflowState.researchResult) && ~fig.UserData.workflowState.boundaryConfirmed);
write(fullfile(out,'ui_checks.json'),struct('adiabaticInputZero',true,'unknownGivenTw',true,'thermalEditInvalidatesResult',true));
fprintf('Classic main workflow and adiabatic field identity passed; no experimental/physical qualification.\n');
end
function waitRender(ui)
t=tic;while toc(t)<40
    drawnow;s=ui.geometry.html.UserData;
    if ~isempty(s.rendered) && s.rendered.revision==s.revision && ~isempty(ui.getSection()),break;end
    pause(.1);
end
assert(~isempty(s.rendered) && s.rendered.revision==s.revision && s.rendered.webglError==0 && ~isempty(ui.getSection()));
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
