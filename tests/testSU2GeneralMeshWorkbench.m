function report=testSU2GeneralMeshWorkbench(out,mode)
%TESTSU2GENERALMESHWORKBENCH 新完整网格从主五步骤启动、接收并定量导出。
% 模块一：原壁/一般列坐标/尾缘扇与已验收几何身份保持，独立初场不用旧场编号。
% 模块二：只启动事前约定的一项150更新，review仅复用其保存结果。
% 模块三：原19752面、周向截面、方法/局部步号及未知求解壁温保持可核对。
if nargin<2,mode='run';end
assert(~isfolder(out));mkdir(out);d=load('outputs/v118_validation/actual01/surface.mat');
c=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=c.cfg;cfg.output.verbose=false;
cfg.su2=struct('volumeMode','hybridTriangles','iterations',150,'outputEvery',75,'namedCheckpoints',true, ...
    'maximumSeconds',600,'checkpointLimitMiB',950,'pythonExecutable','D:/anacondapython/python.exe', ...
    'flowModel','SST-2003m','initialMode','sstWallDistance','initialThicknessM',.001,'turbulenceIntensity',.01, ...
    'turbulentViscosityRatio',10,'turbulentPrandtl',.9, ...
    'turbulenceProvenance','v145 explicit fully turbulent startup assumption; not measured transition');
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));
fig.UserData.actions.useSurfaceForSU2(d.s);fig.UserData.actions.readHybridMesh('outputs/v145_validation/mesh');
info=fig.UserData.volumeInfo;assert(info.nodes==335726 && info.cells==705490 && info.construction.tetsBelowSICN01==0);
assert(strcmp(info.construction.method,'original-triangle-general-prisms-and-tetgen-plc-exterior'));
assert(contains(fig.UserData.controls.meshMethod.Text,'切向偏移'));
fig.UserData.actions.selectStep(3);ctrl=fig.UserData.controls;
assert(strcmp(ctrl.initialMode.Value,'sstWallDistance') && strcmp(ctrl.checkpointSource.Visible,'off'));
fig.UserData.actions.confirmBoundary();assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
fig.UserData.actions.selectStep(4);
if strcmp(mode,'run')
    fig.UserData.actions.run();assert(fig.UserData.workflowState.busy,fig.UserData.status.Text);
    folder=fig.UserData.engineJob.folder;write(fullfile(out,'engine_folder.json'),struct('folder',folder));
    t=tic;while fig.UserData.workflowState.busy && toc(t)<960,drawnow;pause(.2);end
    if fig.UserData.workflowState.busy,fig.UserData.actions.cancel();end
    assert(~fig.UserData.workflowState.busy,'后台未终态，已请求取消。');
    assert(strcmp(fig.UserData.workflowState.status,'research-ready'),fig.UserData.status.Text);
else
    source=jsondecode(fileread('outputs/v145_validation/main01/engine_folder.json'));folder=source.folder;
    fig.UserData.actions.readResearch(fullfile(folder,'mapped_research.mat'));
end
r=fig.UserData.workflowState.researchResult;a=r.convergenceAudit;assert(a.checkpointIdentityConfirmed,a.reason);
assert(a.checkpointIteration==149 && isfield(a,'independentInitialization') && ~isfield(a,'sourceModelLocalIteration'));
assert(~r.productionQualified && numel(r.heatFluxRaw)==19752 && r.thermalStateAudit.interfacePassed);
ui=fig.UserData.researchUI;ui.setSectionX(-1.5);waitGL(ui);ui.inspectFace(12000,[]);
ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui);
ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'section.csv'));section=ui.getSection();assert(size(section.segments,1)==92);
ui.select.Value='wallSolved';ui.select.ValueChangedFcn(ui.select,[]);waitGL(ui);
assert(all(isnan(ui.geometry.html.UserData.payload.faceField.values)));save(fullfile(out,'result.mat'),'r','-v7');
fig.UserData.actions.selectStep(4);fig.UserData.actions.readResearch(fullfile(folder,'mapped_research.mat'));
assert(fig.UserData.workflowState.researchResult.convergenceAudit.checkpointIteration==149);
waitGL(fig.UserData.researchUI);
report=struct('passed',true,'mode',mode,'engineFolder',folder,'nodes',335726,'cells',705490,'originalFaces',19752, ...
    'sectionSegments',92,'independentInitialState',true,'localIteration',149,'sameConditionReopened',true, ...
    'webgl',fig.UserData.researchUI.geometry.html.UserData.rendered,'physicsQualified',false,'newQualifiedArea',0);
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
