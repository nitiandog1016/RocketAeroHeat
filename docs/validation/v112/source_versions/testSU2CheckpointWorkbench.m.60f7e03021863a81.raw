function testSU2CheckpointWorkbench(out,source,capture)
%TESTSU2CHECKPOINTWORKBENCH 实测保存选项、预算限制和超时检查点显式重开。
% 模块一：一次小型真实超时任务；复核时传入原目录，不重复启动求解。
% 模块二：完整原面身份、拒绝不同工况、导出和未知壁温保持研究语义。
if nargin<2,source='';end
if nargin<3,capture=false;end
assert(~isfolder(out));mkdir(out);evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.flight.Mach=4;cfg.flight.alpha=20;cfg.flight.beta=5;cfg.flight.Twall=300;
cfg.su2=struct('volumeMode','starRay','pythonExecutable','D:/anacondapython/python.exe', ...
    'iterations',300,'maximumSeconds',5,'outputEvery',20,'namedCheckpoints',true);
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()close(fig));
fig.UserData.actions.loadNative();fig.UserData.actions.generateMesh();fig.UserData.actions.selectStep(3);
assert(strcmp(fig.UserData.controls.checkpointMode.Value,'indexed'));
assert(contains(fig.UserData.controls.storageEstimate.Text,'MiB'));
% 预算估计越界须在创建进程/输出目录之前拒绝；改回原参数再运行。
c=fig.UserData.controls.outputEvery;c.Value=1;c.ValueChangedFcn(c,[]);
fig.UserData.actions.confirmBoundary();fig.UserData.actions.run();
assert(isempty(fig.UserData.engineJob) && ~fig.UserData.workflowState.busy && contains(fig.UserData.status.Text,'100 MiB'));
fig.UserData.actions.selectStep(3);c=fig.UserData.controls.outputEvery;c.Value=20;c.ValueChangedFcn(c,[]);
fig.UserData.actions.confirmBoundary();
if isempty(source)
    fig.UserData.actions.run();assert(fig.UserData.workflowState.busy,fig.UserData.status.Text);
    start=tic;while fig.UserData.workflowState.busy && toc(start)<45,drawnow;pause(.1);end
    assert(~fig.UserData.workflowState.busy,'实际任务尚未终态。');source=fig.UserData.engineJob.folder;
end
receipt=jsondecode(fileread(fullfile(source,'engine_process.json')));
assert(strcmp(receipt.status,'terminal') && receipt.timedOut && ~receipt.cancelled);
assert(isempty(fig.UserData.workflowState.researchResult),'超时结果被自动接收。');
fig.UserData.actions.readTimeoutCheckpoint(source);state=fig.UserData.workflowState;r=state.researchResult;
assert(strcmp(state.status,'research-ready') && r.receipt.timedOut && ~state.steps(4).completed);
assert(r.convergenceAudit.checkpointIdentityConfirmed && ~r.productionQualified && size(r.F,1)==2448);
assert(strcmp(r.convergenceAudit.checkpointIdentityMethod,'native-indexed-pair-SHA-and-EOS'));
% 普通成功入口仍拒绝真实超时回执；新入口也拒绝未确认身份。
probe=aeroHeatWorkflowState('beginSolve',state);payload=struct('result',r,'solveId',probe.solveId,'sourceRevision',probe.revision);
rejected=false;try,aeroHeatWorkflowState('completeResearch',probe,payload);catch,rejected=true;end;assert(rejected);
payload.result.convergenceAudit.checkpointIdentityConfirmed=false;
rejected=false;try,aeroHeatWorkflowState('completeCheckpointResearch',probe,payload);catch,rejected=true;end;assert(rejected);
ui=fig.UserData.researchUI;assert(contains(ui.auditText.Text,'预算超时') && contains(ui.auditText.Text,'显示检查点'));
ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);ui.setSectionX(.7);
start=tic;while toc(start)<40
    drawnow;web=ui.geometry.html.UserData;
    if ~isempty(web.rendered) && web.rendered.revision==web.revision && ~isempty(ui.getSection()),break;end
    pause(.1);
end
assert(~isempty(web.rendered) && web.rendered.webglError==0 && web.rendered.visiblePixelSamples>0);
ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'section.csv'));
for file={'fields.csv.json','section.csv.json'}
    meta=jsondecode(fileread(fullfile(out,file{1})));assert(isequaln(meta.convergenceAudit,jsondecode(jsonencode(r.convergenceAudit))));
end
ui.inspectFace(1235);
if capture
    write(fullfile(out,'ready_for_capture.json'),struct('source',source));fprintf('V112_CHECKPOINT_UI_READY\n');
    start=tic;while ~isfile(fullfile(out,'capture_done.json')) && toc(start)<60,drawnow;pause(.1);end
end
fig.UserData.actions.selectStep(3);c=fig.UserData.controls.alpha;c.Value=10;c.ValueChangedFcn(c,[]);
fig.UserData.actions.confirmBoundary();refused=false;
try,fig.UserData.actions.readTimeoutCheckpoint(source);catch,refused=true;end
assert(refused && isempty(fig.UserData.workflowState.researchResult));
% 将新输出模式改成覆盖写，菜单状态必须同步并使旧输入失效。
fig.UserData.actions.selectStep(3);c=fig.UserData.controls.checkpointMode;c.Value='overwrite';c.ValueChangedFcn(c,[]);
assert(~fig.UserData.workflowState.cfg.su2.namedCheckpoints && ~fig.UserData.workflowState.boundaryConfirmed);
report=struct('actualTimeoutFolder',source,'noAutomaticAcceptance',true,'explicitReopenPassed',true, ...
    'budgetRejectedBeforeProcess',true,'differentConditionRejected',true,'unconfirmedStateRejected',true, ...
    'ordinarySuccessPathStillRejectsTimeout',true,'rendered',web.rendered,'audit',r.convergenceAudit, ...
    'newQualifiedArea',0,'solverInvoked',nargin<2 || isempty(source));
% solverInvoked明确按入口实参记录，不能因source在运行中被赋值而误报。
report.solverInvoked=~isempty(fig.UserData.engineJob);
write(fullfile(out,'report.json'),report);save(fullfile(out,'reopened.mat'),'r','report','-v7');
fprintf('Main checkpoint options, bounded actual timeout, explicit research reopen and export verified.\n');
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
