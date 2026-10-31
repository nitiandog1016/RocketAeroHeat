function report=reviewActualRocketFullDomain(out)
%REVIEWACTUALROCKETFULLDOMAIN 只读真实完整域终态，把全部原面研究值接回现有工作台。
% 模块一：原面映射、原生具名时刻、热边界与实际终态；超时不得冒充正常结束。
% 模块二：同一组件显示p/q/T/来源，导出三个X截面及全表面数据。
% 模块三：数值有限与全表面物理资格分开；给定Tw与未知求解Tw分开。
assert(~isfolder(out));mkdir(out);folder='outputs/v118_validation/actual01/budget_run';
receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
assert(strcmp(receipt.status,'terminal') && (receipt.exitCode==0 || receipt.timedOut));
surface=load('outputs/v118_validation/actual01/surface.mat');saved=load(fullfile(folder,'projection.mat'));
r=mapSU2WallResult(readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk')),saved.projection);
r.V=surface.s.V;r.F=surface.s.F;r.receipt=receipt;
manifest=jsondecode(fileread(fullfile(folder,'case_manifest.json')));r.flight=manifest.flight;
r.sourceRevision='v118-original19752-exact30km';r.solveId=receipt.request.jobToken;
r=attachSU2ThermalBoundary(r,resolveSU2ThermalBoundary(surface.s,struct('mode','isothermal','temperature',300)));
r.convergenceAudit=auditSU2ResultConvergence(folder,r);
assert(r.convergenceAudit.checkpointIdentityConfirmed && ~r.productionQualified);
assert(r.convergenceAudit.checkpointIteration<=r.convergenceAudit.terminalIteration);
assert(numel(r.pressureRaw)==19752 && abs(sum(r.faceArea)-1.9646292565715553)<1e-12);
fig=uifigure('Name','RocketAeroHeat · 真实20度完整域研究','Position',[20 30 1460 850]);guard=onCleanup(@()close(fig));
host=uipanel(fig,'Units','normalized','Position',[0 0 1 1]);ui=createSU2ResearchView(host,r);
records=cell(3,1);xs=[-.2,-1.5,-2.65];fields={'pressureRaw','heatFluxRaw','fluidNodeTemperatureMeanRaw'};
for k=1:3
    ui.select.Value=fields{k};ui.select.ValueChangedFcn(ui.select,[]);ui.setSectionX(xs(k));waitRender(ui);
    assert(isequal(ui.geometry.html.UserData.payload.faceField.values,r.(fields{k})));
    ui.exportSectionTo(fullfile(out,sprintf('section%d.csv',k)));
    records{k}=struct('X_m',xs(k),'field',fields{k},'rendered',ui.geometry.html.UserData.rendered);
end
ui.inspectFace(12000,[]);assert(contains(strjoin(ui.probe.Value,newline),'官方原生更新'));
ui.exportTo(fullfile(out,'fields.csv'));save(fullfile(out,'result.mat'),'r','-v7');
ui.select.Value='wallInput';ui.select.ValueChangedFcn(ui.select,[]);waitRender(ui);
assert(all(ui.geometry.html.UserData.payload.faceField.values==300));
ui.select.Value='wallSolved';ui.select.ValueChangedFcn(ui.select,[]);waitRender(ui);
assert(all(isnan(ui.geometry.html.UserData.payload.faceField.values)));
report=struct('originalFaces',19752,'totalAreaM2',sum(r.faceArea),'convergence',r.convergenceAudit, ...
    'flight',r.flight,'sections',{records},'processExitCode',receipt.exitCode,'timedOut',receipt.timedOut, ...
    'fullPhysicalAreaPowerRaw_W',r.targetHeatPowerRaw,'nativeNodalPowerRaw_W',r.conservativeHeatLoadPowerRaw, ...
    'productionResultComponentReused',true,'mainWorkflowAcceptanceClaimed',false, ...
    'newQualifiedArea',0,'physicsQualified',false,'solverRerun',false);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Actual original 19752 faces, three sections, runtime and given/unknown Tw checked; qualification remains separate.\n');
end
function waitRender(ui)
t=tic;while toc(t)<40
    drawnow;s=ui.geometry.html.UserData;
    if ~isempty(s.rendered)&&s.rendered.revision==s.revision&&~isempty(ui.getSection()),break;end
    pause(.1);
end
assert(~isempty(s.rendered)&&s.rendered.revision==s.revision&&s.rendered.webglError==0&&~isempty(ui.getSection()));
end
