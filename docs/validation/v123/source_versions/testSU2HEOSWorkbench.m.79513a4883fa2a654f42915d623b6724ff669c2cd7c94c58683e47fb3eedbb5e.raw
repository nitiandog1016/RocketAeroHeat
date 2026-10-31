function report=testSU2HEOSWorkbench(out,mode,capture)
%TESTSU2HEOSWORKBENCH 主五步骤真实HEOS流程及已存场复用，保持研究资格。
% 模块一：检查显式物性/来源/预算入口、失败和后台物性查询取消。
% 模块二：仅run模式启动一次完整域两步任务；review模式从原目录只读恢复。
% 模块三：原19752面、探针、截面、导出和未知壁温；可等待原生窗口检查。
if nargin<2,mode='run';end
if nargin<3,capture=false;end
assert(~isfolder(out));mkdir(out);
d=load('outputs/v118_validation/actual01/surface.mat');c=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=c.cfg;
cfg.output.verbose=false;cfg.su2=struct('volumeMode','hybridTriangles','iterations',2,'outputEvery',1, ...
    'namedCheckpoints',true,'maximumSeconds',180,'checkpointLimitMiB',256,'pythonExecutable','D:/anacondapython/python.exe');
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));
fig.UserData.actions.useSurfaceForSU2(d.s);fig.UserData.actions.readHybridMesh('outputs/v118_validation/actual01/mesh');
fig.UserData.actions.selectStep(3);controls=fig.UserData.controls;
assert(strcmp(controls.thermodynamics.Value,'idealGas') && strcmp(controls.initialMode.Value,'freestream'));
fig.UserData.actions.confirmBoundary();s=fig.UserData.workflowState;
assert(s.boundaryConfirmed && ~isfield(s.cfg.su2,'thermodynamicModel') && ~isfield(s.cfg.su2,'heosSource'));
fig.UserData.actions.selectStep(3);controls=fig.UserData.controls;
controls.thermodynamics.Value='HEOS-Air';controls.thermodynamics.ValueChangedFcn(controls.thermodynamics,[]);
assert(contains(fig.UserData.status.Text,'请选择已完成的定比热源场') && isempty(fig.UserData.engineJob));
assert(~fig.UserData.workflowState.boundaryConfirmed && ~fig.UserData.workflowState.steps(4).available);
source=fullfile(pwd,'outputs','v120_validation','study01','control');
controls.heosSource.Value=source;controls.heosSource.ValueChangedFcn(controls.heosSource,[]);
fig.UserData.actions.confirmBoundary();assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
assert(strcmp(fig.UserData.workflowState.cfg.su2.initialMode,'propertyTransfer'));
engineFolder='';failureFolder='';cancelFolder='';
if strcmp(mode,'run')
    % 故意选择另一个HEOS场，真实后台须在准备阶段拒绝，不能伪造SU2已运行。
    fig.UserData.actions.selectStep(3);controls=fig.UserData.controls;
    controls.heosSource.Value=fullfile(pwd,'outputs','v122_validation','actual04');
    controls.heosSource.ValueChangedFcn(controls.heosSource,[]);fig.UserData.actions.confirmBoundary();fig.UserData.actions.run();
    failureFolder=fig.UserData.engineJob.folder;waitTerminal(fig,40);
    assert(strcmp(fig.UserData.workflowState.status,'failed') && isempty(fig.UserData.workflowState.researchResult));
    p=jsondecode(fileread(fullfile(failureFolder,'workflow','pipeline_process.json')));
    assert(strcmp(p.stage,'preparation') && p.exitCode~=0 && ~p.engineStarted && ~isfile(fullfile(failureFolder,'engine_process.json')));
    % 等真实物性子进程运行后取消；主窗口仍能响应，原生引擎尚未启动。
    fig.UserData.actions.selectStep(3);controls=fig.UserData.controls;
    controls.heosSource.Value=source;controls.heosSource.ValueChangedFcn(controls.heosSource,[]);
    fig.UserData.actions.confirmBoundary();fig.UserData.actions.run();cancelFolder=fig.UserData.engineJob.folder;
    t=tic;found=false;
    while toc(t)<45 && fig.UserData.workflowState.busy
        drawnow;path=fullfile(cancelFolder,'initial_thermodynamics','PT_process.json');
        if isfile(path),p=jsondecode(fileread(path));if strcmp(p.status,'running'),found=true;break;end,end
        pause(.05);
    end
    assert(found && strcmp(fig.UserData.controls.cancel.Enable,'on'));
    fig.UserData.actions.cancel();waitTerminal(fig,20);p=jsondecode(fileread(fullfile(cancelFolder,'workflow','pipeline_process.json')));
    assert(p.cancelled && p.exitCode==125 && ~p.engineStarted && ~isfile(fullfile(cancelFolder,'engine_process.json')));
    q=jsondecode(fileread(fullfile(cancelFolder,'initial_thermodynamics','PT_process.json')));
    assert(strcmp(q.status,'terminal') && q.cancelled && q.exitCode~=0 && fig.UserData.engineJob.process.HasExited);
    % 唯一实际完整域两步任务；失败时不自动再次启动。
    fig.UserData.actions.run();assert(fig.UserData.workflowState.busy,fig.UserData.status.Text);
    engineFolder=fig.UserData.engineJob.folder;write(fullfile(out,'engine_folder.json'),struct('folder',engineFolder));
    waitTerminal(fig,370);s=fig.UserData.workflowState;assert(strcmp(s.status,'research-ready'),fig.UserData.status.Text);
    short=s.researchResult;assert(short.convergenceAudit.checkpointIdentityConfirmed && short.convergenceAudit.checkpointIteration==1);
    assert(short.convergenceAudit.sourceCumulativeIteration==799 && ~isfield(short.convergenceAudit,'cumulativeCheckpointIteration'));
    p=short.workflowReceipt;assert(p.exitCode==0 && numel(p.phases)==3 && p.engineStarted && ~p.cancelled);
    assert(short.receipt.elapsedSeconds<=180 && short.flight.referenceMach==5 && abs(short.flight.Mach-4.998034974966865)<1e-12);
    ui=fig.UserData.researchUI;ui.setSectionX(-1.5);waitGL(ui);ui.exportTo(fullfile(out,'short_fields.csv'));
    save(fullfile(out,'short_result.mat'),'short','-v7');
    % 同一当前快照重开沿用真实后台终态，不重新转换/启动。
    fig.UserData.actions.selectStep(4);fig.UserData.actions.readResearch(fullfile(engineFolder,'mapped_research.mat'));
    assert(fig.UserData.workflowState.researchResult.convergenceAudit.checkpointIteration==1);
end
fig.UserData.actions.selectStep(4);fig.UserData.actions.readVerifiedHEOS(fullfile(pwd,'outputs','v122_validation','actual04'));
s=fig.UserData.workflowState;r=s.researchResult;assert(r.convergenceAudit.checkpointIteration==19 && r.convergenceAudit.sourceCumulativeIteration==799);
assert(~r.productionQualified && ~r.convergenceAudit.checkpointResidualPassed && r.verifiedCaseReview.solverRerun==false);
ui=fig.UserData.researchUI;ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);ui.setSectionX(-1.5);waitGL(ui);
ui.inspectFace(12000,[]);assert(contains(strjoin(ui.probe.Value,newline),'CoolProp'));
assert(contains(ui.geometry.html.UserData.payload.caption,'EOS M'));
ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'section.csv'));section=ui.getSection();
assert(size(section.segments,1)==92 && numel(r.heatFluxRaw)==19752);
ui.select.Value='wallSolved';ui.select.ValueChangedFcn(ui.select,[]);assert(all(isnan(ui.geometry.html.UserData.payload.faceField.values)));
ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui);
save(fullfile(out,'result.mat'),'r','-v7');
rendered=ui.geometry.html.UserData.rendered;
if capture
    write(fullfile(out,'ready_for_capture.json'),struct('figureName',fig.Name,'timeLimitSeconds',180));fprintf('V123_MAIN_UI_READY\n');
    t=tic;while ~isfile(fullfile(out,'capture_done.json')) && toc(t)<180,drawnow;pause(.1);end
end
% 相同表面但不同物理边界拒绝，避免“查看旧场”绕过当前工况。
fig.UserData.actions.selectStep(3);controls=fig.UserData.controls;
controls.alpha.Value=21;controls.alpha.ValueChangedFcn(controls.alpha,[]);fig.UserData.actions.confirmBoundary();refused=false;
try,fig.UserData.actions.readVerifiedHEOS(fullfile(pwd,'outputs','v122_validation','actual04'));catch,refused=true;end
assert(refused && isempty(fig.UserData.workflowState.researchResult));
old=load('outputs/v119_validation/review01/result.mat','r');lines=describeSU2Convergence(old.r);assert(contains(lines{2},'累计 599'));
report=struct('mode',mode,'mainFiveStepWorkflow',true,'engineFolder',engineFolder,'preparationFailureFolder',failureFolder, ...
    'preparationCancellationFolder',cancelFolder,'existingLocal19Reused',true,'defaultSnapshotFieldsPreserved',true, ...
    'wrongPhysicalConditionRejected',true,'invalidEditRevokesOldConfirmation',true,'legacyCumulative599Preserved',true,'originalFaces',19752,'sectionRows',92, ...
    'webgl',rendered,'givenTwDistinctFromUnknownSolvedTw',true,'physicsQualified',false,'newQualifiedArea',0);
write(fullfile(out,'report.json'),report);fprintf('Main HEOS workflow and verified-case review passed; new qualified area 0.\n');
end
function waitTerminal(fig,limit)
t=tic;while fig.UserData.workflowState.busy && toc(t)<limit,drawnow;pause(.1);end
if fig.UserData.workflowState.busy,fig.UserData.actions.cancel();end
assert(~fig.UserData.workflowState.busy,'主工作台尚未收到终态；已提交取消，不自动重跑。');
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
