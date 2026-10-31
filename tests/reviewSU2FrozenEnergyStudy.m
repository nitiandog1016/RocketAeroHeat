function report=reviewSU2FrozenEnergyStudy(out)
%REVIEWSU2FROZENENERGYSTUDY 只读真实终态，复用原结果组件核对运行时来源与原面导出。
% 模块一：具名检查点和原三角映射；给定绝热零不当作热流预测。
% 模块二：实验/官方来源贯穿三维云图、探针、截面和CSV元数据。
% 模块三：保留未知求解壁温；不重跑引擎，也不授予物理资格。
assert(~isfolder(out));mkdir(out);[surface,~]=buildSchaeferFergusonBody();
folder=fullfile('outputs','v117_validation','apply');
if ~isfile(fullfile(folder,'engine_process.json')),folder=fullfile('outputs','v117_validation','diagnostic');end
receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
assert(strcmp(receipt.status,'terminal') && receipt.exitCode==0 && ~receipt.timedOut);
saved=load(fullfile(folder,'projection.mat'));manifest=jsondecode(fileread(fullfile(folder,'case_manifest.json')));
r=mapSU2WallResult(readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk')),saved.projection);
r.V=surface.V;r.F=surface.F;r.receipt=receipt;r.flight=manifest.flight;
r.sourceRevision='v117-frozen-energy-research';r.solveId=receipt.request.jobToken;
r=attachSU2ThermalBoundary(r,resolveSU2ThermalBoundary(surface,struct('mode','adiabatic')));
r.convergenceAudit=auditSU2ResultConvergence(folder,r);
assert(r.convergenceAudit.checkpointIdentityConfirmed && ~r.productionQualified);
assert(r.convergenceAudit.checkpointIteration==manifest.flight.iterations-1);
info=describeSU2Runtime(receipt);assert(strcmp(info.flavor,'rah-frozen-energy-v1'));
fig=uifigure('Name','RocketAeroHeat · 更新机制研究复核','Position',[20 30 1460 850]);guard=onCleanup(@()close(fig));
host=uipanel(fig,'Units','normalized','Position',[0 0 1 1]);ui=createSU2ResearchView(host,r);
ui.setSectionX(.158496);waitRender(ui);ui.inspectFace(1454,[]);
assert(contains(strjoin(ui.probe.Value,newline),info.method));
assert(contains(ui.geometry.html.UserData.payload.sectionView.method,info.method));
assert(isequal(ui.geometry.html.UserData.payload.faceField.values,r.pressureRaw));
ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'section.csv'));
for name={'fields.csv.json','section.csv.json'}
    m=jsondecode(fileread(fullfile(out,name{1})));assert(isequal(m.algorithmRuntime,info));
end
save(fullfile(out,'result.mat'),'r','-v7');
rendered=ui.geometry.html.UserData.rendered;
ui.select.Value='wallSolved';ui.select.ValueChangedFcn(ui.select,[]);waitRender(ui);
assert(all(isnan(ui.geometry.html.UserData.payload.faceField.values)));
% 旧官方快照沿用相同组件，不依赖刚构建运行时替换其历史身份。
old=load(fullfile('outputs','v116_validation','review02','control_result.mat'));
delete(host.Children);legacy=createSU2ResearchView(host,old.r);waitRender(legacy);legacy.inspectFace(1,[]);
assert(contains(strjoin(legacy.probe.Value,newline),'官方原生更新'));
report=struct('originalFaces',size(r.F,1),'sourceFolder',folder,'algorithmRuntime',info, ...
    'convergence',r.convergenceAudit,'rendered',rendered,'oldOfficialSnapshotPassed',true, ...
    'unknownWallTemperaturePreserved',true,'mainWorkflowAcceptanceClaimed',false, ...
    'newQualifiedArea',0,'physicsQualified',false,'solverRerun',false);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Frozen-energy original-face display, runtime identity, CSV, unknown Tw and old official snapshot passed.\n');
end
function waitRender(ui)
t=tic;while toc(t)<40
    drawnow;s=ui.geometry.html.UserData;
    if ~isempty(s.rendered)&&s.rendered.revision==s.revision&&~isempty(ui.getSection()),break;end
    pause(.1);
end
assert(~isempty(s.rendered)&&s.rendered.revision==s.revision&&s.rendered.webglError==0&&~isempty(ui.getSection()));
end
