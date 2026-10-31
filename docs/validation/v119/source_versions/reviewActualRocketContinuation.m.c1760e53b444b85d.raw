function report=reviewActualRocketContinuation(out)
%REVIEWACTUALROCKETCONTINUATION 将续算真实终态接回原面、累计身份和科研组件。
% 模块一：局部/累计迭代、热边界及源场身份全部沿实际回执核对。
% 模块二：旧快照文案保持；新云图/探针/截面和导出绑定同一未合格研究场。
assert(~isfolder(out));mkdir(out);folder='outputs/v119_validation/continuation01';
receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
assert(strcmp(receipt.status,'terminal') && (receipt.exitCode==0 || receipt.timedOut));
old=load('outputs/v118_validation/review01/result.mat','r');projection=load(fullfile(folder,'projection.mat'),'projection');
r=mapSU2WallResult(readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk')),projection.projection);
r.V=old.r.V;r.F=old.r.F;r.receipt=receipt;r.flight=old.r.flight;
r.sourceRevision='v119-original19752-checkpoint-continuation';r.solveId=receipt.request.jobToken;
r=attachSU2ThermalBoundary(r,old.r.thermalBoundary);
r.convergenceAudit=auditSU2ResultConvergence(folder,r);a=r.convergenceAudit;
assert(a.checkpointIdentityConfirmed,a.reason);
assert(a.completedUpdatesBeforeRun==300 && a.cumulativeCheckpointIteration==300+a.checkpointIteration);
assert(~r.productionQualified && numel(r.pressureRaw)==19752);
lines=describeSU2Convergence(old.r);assert(contains(lines{2},'显示检查点 299 步') && ~contains(lines{2},'累计'));
fig=uifigure('Name','RocketAeroHeat · 完整域续算研究','Position',[30 30 1460 850]);g=onCleanup(@()close(fig));
host=uipanel(fig,'Units','normalized','Position',[0 0 1 1]);ui=createSU2ResearchView(host,r);
assert(contains(ui.auditText.Text,'累计') && contains(ui.auditText.Text,num2str(a.cumulativeCheckpointIteration)));
ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);ui.setSectionX(-1.5);
t=tic;while toc(t)<40
    drawnow;s=ui.geometry.html.UserData;
    if ~isempty(s.rendered)&&s.rendered.revision==s.revision&&~isempty(ui.getSection()),break;end
    pause(.1);
end
assert(s.rendered.webglError==0 && s.rendered.revision==s.revision && ~isempty(ui.getSection()));
ui.inspectFace(12000,[]);ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'section.csv'));
for name={'fields.csv.json','section.csv.json'}
    m=jsondecode(fileread(fullfile(out,name{1})));assert(isequaln(m.convergenceAudit,jsondecode(jsonencode(a))));
end
save(fullfile(out,'result.mat'),'r','-v7');
section=ui.getSection();
report=struct('audit',a,'originalFaces',19752,'fullArea',sum(r.faceArea), ...
    'givenTwRange',[min(r.wallTemperatureInput),max(r.wallTemperatureInput)],'legacySnapshotIdentityPreserved',true, ...
    'componentWebGLPassed',true,'sectionRows',size(section.segments,1), ...
    'newQualifiedArea',0,'physicsQualified',false,'mainEndToEndLaunchClaimed',false,'solverRerun',false);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);g2=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Original faces and local/cumulative identities bound to WebGL and exports. No qualification implied.\n');
end
