function report=testSU2SSTOutputWorkbench(out,mode,capture)
%TESTSU2SSTOUTPUTWORKBENCH 实际主工作台接收非紧凑SST多步场，并验证同工况重开。
% 模块一：只运行事前计划的完整域50步，不自动重试；review模式只读原任务。
% 模块二：模型选择置顶/参数按需显示，阶段p/T与更新k/ω分别审核。
% 模块三：源局部0/上游99与本段49分列，原面/截面/探针/未知壁温及导出。
if nargin<2,mode='run';end
if nargin<3,capture=false;end
assert(~isfolder(out));mkdir(out);
d=load('outputs/v118_validation/actual01/surface.mat');c=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=c.cfg;
cfg.output.verbose=false;cfg.su2=struct('volumeMode','hybridTriangles','iterations',50,'outputEvery',25, ...
    'namedCheckpoints',true,'maximumSeconds',240,'checkpointLimitMiB',950,'pythonExecutable','D:/anacondapython/python.exe', ...
    'flowModel','SST-2003m','initialMode','sstTransfer','checkpointSource',fullfile(pwd,'outputs','v125_validation','phase01'));
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));
fig.UserData.actions.useSurfaceForSU2(d.s);fig.UserData.actions.readHybridMesh('outputs/v120_validation/mesh32');
fig.UserData.actions.selectStep(3);controls=fig.UserData.controls;
assert(controls.flowModel.Layout.Row==2 && strcmp(controls.heosSource.Visible,'off'));
controls.flowModel.Value='laminar';controls.flowModel.ValueChangedFcn(controls.flowModel,[]);
assert(strcmp(controls.turbulenceIntensity.Visible,'off'));
controls.flowModel.Value='SST-2003m';controls.flowModel.ValueChangedFcn(controls.flowModel,[]);
assert(strcmp(controls.turbulenceIntensity.Visible,'on') && numel(controls.boundaryLayout.grid.RowHeight)<40);
fig.UserData.actions.confirmBoundary();assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
if capture
    fig.UserData.actions.selectStep(3);waitCapture(fig,out,'boundary');fig.UserData.actions.selectStep(4);
end
if strcmp(mode,'run')
    fig.UserData.actions.run();assert(fig.UserData.workflowState.busy,fig.UserData.status.Text);
    folder=fig.UserData.engineJob.folder;write(fullfile(out,'engine_folder.json'),struct('folder',folder));
    t=tic;while fig.UserData.workflowState.busy && toc(t)<650,drawnow;pause(.1);end
    if fig.UserData.workflowState.busy,fig.UserData.actions.cancel();end
    assert(~fig.UserData.workflowState.busy,'后台尚未终态，已请求取消。');
    assert(strcmp(fig.UserData.workflowState.status,'research-ready'),fig.UserData.status.Text);
else
    q=jsondecode(fileread('outputs/v126_validation/main01/engine_folder.json'));folder=q.folder;
    fig.UserData.actions.readResearch(fullfile(folder,'mapped_research.mat'));
end
r=fig.UserData.workflowState.researchResult;a=r.convergenceAudit;
assert(a.checkpointIdentityConfirmed,a.reason);assert(a.checkpointIteration==49 && a.sourceModelLocalIteration==0);
assert(isfield(a,'turbulenceRestart') && a.turbulenceRestart.sourceWasDiagnostic);
assert(a.turbulenceRestart.upstreamInitialization.sourceLocalIteration==99 && ~isfield(a,'cumulativeCheckpointIteration'));
assert(strcmp(r.thermalStateAudit.thermalStatePhase,'native-noncompact-flow-primitive-before-turbulence-v1'));
assert(~r.productionQualified && numel(r.heatFluxRaw)==19752);
ui=fig.UserData.researchUI;ui.setSectionX(-1.5);waitGL(ui);ui.inspectFace(12000,[]);
assert(contains(strjoin(ui.probe.Value,newline),'分别审核'));
ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui);
ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'section.csv'));section=ui.getSection();assert(size(section.segments,1)==92);
ui.select.Value='wallSolved';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui);
assert(all(isnan(ui.geometry.html.UserData.payload.faceField.values)));
ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui);
if capture,waitCapture(fig,out,'result');end
save(fullfile(out,'result.mat'),'r','-v7');
fig.UserData.actions.selectStep(4);fig.UserData.actions.readResearch(fullfile(folder,'mapped_research.mat'));
assert(fig.UserData.workflowState.researchResult.convergenceAudit.checkpointIteration==49);
% 重开会替换旧研究视图；报告必须读取新视图，不能再访问已删除的HTML句柄。
ui=fig.UserData.researchUI;waitGL(ui);
report=struct('mode',mode,'engineFolder',folder,'mainPipelineAccepted',true,'sameConditionReopened',true, ...
    'sourceDiagnosticLocalIteration',0,'upstreamSSTLocalIteration',99,'checkpointLocalIteration',49, ...
    'modelControlsAtTop',true,'irrelevantControlsHidden',true,'originalFaces',19752,'sectionSegments',92, ...
    'webgl',ui.geometry.html.UserData.rendered,'physicsQualified',false,'newQualifiedArea',0);
write(fullfile(out,'report.json'),report);
end
function waitGL(ui)
t=tic;while toc(t)<40
    drawnow;s=ui.geometry.html.UserData;
    if ~isempty(s.rendered)&&s.rendered.revision==s.revision&&~isempty(ui.getSection()),break;end
    pause(.1);
end
assert(~isempty(s.rendered)&&s.rendered.revision==s.revision&&s.rendered.faces==19752&&s.rendered.webglError==0&&~isempty(ui.getSection()));
end
function waitCapture(fig,out,stage)
write(fullfile(out,[stage '_ready.json']),struct('figureName',fig.Name));
t=tic;while ~isfile(fullfile(out,[stage '_captured.json'])) && toc(t)<180,drawnow;pause(.1);end
end
function closeIfPresent(fig)
if isgraphics(fig),close(fig);end
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
