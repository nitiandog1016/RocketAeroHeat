function report=testSU2SSTWorkbench(out,mode,capture)
%TESTSU2SSTWORKBENCH 主五步骤SST真实完整域研究与原面结果，不自动重试数值。
% 模块一：复用已存32层混合域，显式模型/入口/Prt/来源，HEOS组合拒绝。
% 模块二：run只做一次100步实际任务；review只读原成功目录，不补跑流场。
% 模块三：检查k/ω残差、模型步号、原19752面/截面/CSV与未知材料壁温。
if nargin<2,mode='run';end
if nargin<3,capture=false;end
assert(~isfolder(out));mkdir(out);
d=load('outputs/v118_validation/actual01/surface.mat');c=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=c.cfg;
cfg.output.verbose=false;cfg.su2=struct('volumeMode','hybridTriangles','iterations',100,'outputEvery',50, ...
    'namedCheckpoints',true,'maximumSeconds',600,'checkpointLimitMiB',600,'pythonExecutable','D:/anacondapython/python.exe');
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));
fig.UserData.actions.useSurfaceForSU2(d.s);fig.UserData.actions.readHybridMesh('outputs/v120_validation/mesh32');
fig.UserData.actions.selectStep(3);controls=fig.UserData.controls;
fig.UserData.actions.confirmBoundary();assert(~isfield(fig.UserData.workflowState.cfg.su2,'flowModel'));
controls=fig.UserData.controls;
controls.flowModel.Value='SST-2003m';controls.flowModel.ValueChangedFcn(controls.flowModel,[]);
assert(~fig.UserData.workflowState.boundaryConfirmed,'缺源场时必须撤销边界确认。');
controls.checkpointSource.Value=fullfile(pwd,'outputs','v120_validation','study01','candidate');
controls.checkpointSource.ValueChangedFcn(controls.checkpointSource,[]);fig.UserData.actions.confirmBoundary();
assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
assert(strcmp(fig.UserData.workflowState.cfg.su2.initialMode,'sstTransfer'));
controls=fig.UserData.controls;
% 当前物性/流动组合应拒绝，不能静默降为另一种计算。
controls.thermodynamics.Value='HEOS-Air';controls.thermodynamics.ValueChangedFcn(controls.thermodynamics,[]);
assert(~fig.UserData.workflowState.boundaryConfirmed && contains(fig.UserData.status.Text,'HEOS RANS'));
controls.thermodynamics.Value='idealGas';controls.thermodynamics.ValueChangedFcn(controls.thermodynamics,[]);
controls.flowModel.ValueChangedFcn(controls.flowModel,[]);fig.UserData.actions.confirmBoundary();
assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
if strcmp(mode,'run')
    fig.UserData.actions.run();assert(fig.UserData.workflowState.busy,fig.UserData.status.Text);
    folder=fig.UserData.engineJob.folder;write(fullfile(out,'engine_folder.json'),struct('folder',folder));
    t=tic;while fig.UserData.workflowState.busy && toc(t)<950,drawnow;pause(.1);end
    if fig.UserData.workflowState.busy,fig.UserData.actions.cancel();end
    assert(~fig.UserData.workflowState.busy,'实际流程尚未终态，已请求取消。');
    assert(strcmp(fig.UserData.workflowState.status,'research-ready'),fig.UserData.status.Text);
else
    q=jsondecode(fileread('outputs/v124_validation/workbench02/engine_folder.json'));folder=q.folder;
    fig.UserData.actions.selectStep(4);fig.UserData.actions.readResearch(fullfile(folder,'mapped_research.mat'));
end
r=fig.UserData.workflowState.researchResult;a=r.convergenceAudit;
assert(a.checkpointIdentityConfirmed,a.reason);assert(a.checkpointIteration==99 && numel(a.checkpointResiduals)==4);
assert(a.sourceModelLocalIteration==199 && ~isfield(a,'cumulativeCheckpointIteration'));
assert(a.turbulenceTransfer.upstreamBeforeMeshTransferIteration==599 && isempty(a.turbulenceTransfer.sourceCumulativeIteration));
assert(~r.productionQualified && numel(r.heatFluxRaw)==19752 && abs(sum(r.faceArea)-1.9646292565715553)<1e-12);
ui=fig.UserData.researchUI;ui.setSectionX(-1.5);waitGL(ui);
ui.inspectFace(12000,[]);assert(contains(strjoin(ui.probe.Value,newline),'SST-2003m'));
ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'section.csv'));section=ui.getSection();
assert(size(section.segments,1)==92);
ui.select.Value='wallSolved';ui.select.ValueChangedFcn(ui.select,[]);assert(all(isnan(ui.geometry.html.UserData.payload.faceField.values)));
ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui);
save(fullfile(out,'result.mat'),'r','-v7');
rendered=ui.geometry.html.UserData.rendered;
if capture
    write(fullfile(out,'ready_for_capture.json'),struct('figureName',fig.Name,'timeLimitSeconds',180));
    t=tic;while ~isfile(fullfile(out,'capture_done.json')) && toc(t)<180,drawnow;pause(.1);end
end
% 旧层流/HEOS审核按实际原场只读回归，防止新列名规则误伤已有路径。
old=load('outputs/v123_validation/ui02/result.mat','r');oldAudit=auditSU2ResultConvergence('outputs/v122_validation/actual04',old.r);
assert(oldAudit.checkpointIdentityConfirmed && oldAudit.checkpointIteration==19,oldAudit.reason);
fig.UserData.actions.selectStep(4);fig.UserData.actions.readResearch(fullfile(folder,'mapped_research.mat'));
assert(fig.UserData.workflowState.researchResult.convergenceAudit.checkpointIteration==99);
report=struct('mode',mode,'engineFolder',folder,'mainFiveStepSST',true,'HEOSRANSRejected',true, ...
    'defaultFlowFieldAbsent',true,'existingHEOSReadOnlyPassed',true,'originalFaces',19752,'sectionRows',92, ...
    'webgl',rendered,'physicsQualified',false,'newQualifiedArea',0);
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
function closeIfPresent(fig)
if isgraphics(fig),close(fig);end
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
