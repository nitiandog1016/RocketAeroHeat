function testEmbeddedSU2Workbench(out,source,capture)
%TESTEMBEDDEDSU2WORKBENCH 真实主窗口五步骤/进程/取消/原面研究查看验收。
% 模块一：使用短迭代真实SU2检查界面链路，不将其当收敛物理算例。
% 模块二：核对WebGL、探针、CSV、输入失效及关闭时子进程实际退出。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if nargin<2,source='';end
if nargin<3,capture=false;end
assert(~isfolder(out));mkdir(out);
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.flight.Mach=4;cfg.flight.alpha=20;cfg.flight.beta=5;cfg.flight.Twall=300;
cfg.su2=struct('pythonExecutable','D:/anacondapython/python.exe','iterations',2,'maximumSeconds',60);
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));
assert(numel(fig.UserData.workflowButtons)==5 && isempty(fig.UserData.workflowState.result));
fig.UserData.actions.loadNative();s=fig.UserData.workflowState;
assert(s.currentStep==2 && size(s.rawMesh.F,1)==2448);
fig.UserData.actions.generateMesh();s=fig.UserData.workflowState;
assert(s.steps(2).completed && isequal(s.mesh.V,s.rawMesh.V));
fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();
assert(fig.UserData.workflowState.currentStep==4);
if isempty(source)
    fig.UserData.actions.run();s=fig.UserData.workflowState;assert(s.busy);
    job=fig.UserData.engineJob;assert(~isempty(job));
    assert(strcmp(fig.UserData.controls.cancel.Enable,'on') && strcmp(fig.UserData.editorPanel.Enable,'on'));
    awaitTerminal(fig,90);
else
    fig.UserData.actions.readResearch(source);job=struct('folder',fileparts(source));
end
s=fig.UserData.workflowState;
assert(strcmp(s.status,'research-ready'),fig.UserData.status.Text);
assert(isempty(s.result) && ~isempty(s.researchResult) && ~s.steps(4).completed);
r=s.researchResult;assert(~r.productionQualified && ~r.eligibleForSurfaceCoverage);
ui=fig.UserData.researchUI;waitGL(ui.geometry,2448);initialRender=ui.geometry.html.UserData.rendered;
ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui.geometry,2448);
assert(isequal(ui.geometry.html.UserData.payload.faceField.values,r.heatFluxRaw));
ui.low.Value=-10000;ui.high.Value=10000;ui.apply.ButtonPushedFcn(ui.apply,[]);
assert(isequal(ui.geometry.html.UserData.payload.faceField.range,[-10000 10000]));
ui.inspectFace(123);assert(contains(strjoin(ui.probe.Value,newline),'原面元 123'));
ui.select.Value='wallSolved';ui.select.ValueChangedFcn(ui.select,[]);
assert(all(isnan(ui.geometry.html.UserData.payload.faceField.values)));
ui.select.Value='fluidNodeTemperatureMeanRaw';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui.geometry,2448);
csv=fullfile(out,'research_snapshot.csv');ui.exportTo(csv);t=readtable(csv);
assert(max(abs(t.heatFluxRaw_W_m2-r.heatFluxRaw))<1e-7*max(1,max(abs(r.heatFluxRaw))));
assert(all(isnan(t.solvedWallTemperature_K)) && ~any(t.productionQualified));
% MATLAB exportapp对该uihtml宿主失败；保留失败记录，使用原窗口截图核对。
if capture
    write(fullfile(out,'ready_for_capture.json'),struct('figureName',fig.Name));
    fprintf('V107_UI_READY_FOR_CAPTURE\n');start=tic;
    while ~isfile(fullfile(out,'capture_done.json')) && toc(start)<60,drawnow;pause(.1);end
end
% 浏览不重算；边界修改后研究结果立即失效。
id=s.solveId;rev=s.revision;fig.UserData.actions.selectStep(3);fig.UserData.actions.selectStep(5);
s=fig.UserData.workflowState;assert(s.solveId==id && s.revision==rev);
fig.UserData.actions.selectStep(3);c=fig.UserData.controls.alpha;c.Value=10;c.ValueChangedFcn(c,[]);
s=fig.UserData.workflowState;assert(isempty(s.researchResult) && ~s.steps(5).available);
% 同几何不同边界不得重开为当前场。
fig.UserData.actions.confirmBoundary();refused=false;
try,fig.UserData.actions.readResearch(fullfile(job.folder,'mapped_research.mat'));catch,refused=true;end
assert(refused && isempty(fig.UserData.workflowState.researchResult));fig.UserData.actions.selectStep(3);
% 新任务明确改迭代预算，等待实际running之后取消；不依赖固定睡眠猜进程已启动。
c=fig.UserData.controls.iterations;c.Value=30000;c.ValueChangedFcn(c,[]);
fig.UserData.actions.confirmBoundary();fig.UserData.actions.run();awaitRunning(fig,30);
cancelFolder=fig.UserData.engineJob.folder;fig.UserData.actions.cancel();awaitTerminal(fig,30);
s=fig.UserData.workflowState;assert(strcmp(s.status,'failed') && isempty(s.researchResult));
cancelled=jsondecode(fileread(fullfile(cancelFolder,'engine_process.json')));assert(cancelled.cancelled && ~cancelled.timedOut);
% 窗口关闭同样请求取消并等待监督器退出；不能留下后台CPU任务。
fig.UserData.actions.run();awaitRunning(fig,30);closeJob=fig.UserData.engineJob;close(fig);
start=tic;while isgraphics(fig) && toc(start)<30,drawnow;pause(.1);end
assert(~isgraphics(fig) && closeJob.process.HasExited);
closed=jsondecode(fileread(fullfile(closeJob.folder,'engine_process.json')));assert(closed.cancelled);
report=struct('fiveStepMainWindow',true,'originalFaceCount',2448,'actualEngineFolder',job.folder, ...
    'cancelFolder',cancelFolder,'closeFolder',closeJob.folder,'rendered',initialRender, ...
    'csvBoundToSameSnapshot',true,'researchDistinctFromProduction',true,'newQualifiedArea',0);
write(fullfile(out,'report.json'),report);save(fullfile(out,'ui_result.mat'),'r','report','-v7');
fprintf('Main five-step embedded SU2: actual short run, same-face research WebGL/CSV, invalidation, cancellation and close terminal passed. No physical qualification.\n');
end
function awaitTerminal(fig,limit)
start=tic;while fig.UserData.workflowState.busy && toc(start)<limit,drawnow;pause(.1);end
assert(~fig.UserData.workflowState.busy,'主工作台未取得终态。');
end
function awaitRunning(fig,limit)
start=tic;
while toc(start)<limit
    drawnow;file=fullfile(fig.UserData.engineJob.folder,'engine_process.json');
    if isfile(file),s=jsondecode(fileread(file));if isfield(s,'status') && strcmp(s.status,'running'),return;end,end
    pause(.1);
end
error('引擎未进入实际running状态。');
end
function waitGL(view,n)
start=tic;
while toc(start)<40
    drawnow;s=view.html.UserData;
    if ~isempty(s.rendered) && s.rendered.revision==s.revision && s.rendered.faces==n
        assert(s.rendered.webglError==0 && s.rendered.visiblePixelSamples>0);return;
    end
    pause(.1);
end
error('当前版本未完成WebGL实际绘制。');
end
function closeIfPresent(fig)
if isgraphics(fig),close(fig);end
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
