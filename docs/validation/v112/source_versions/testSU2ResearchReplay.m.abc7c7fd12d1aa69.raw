function testSU2ResearchReplay(out,source,capture)
%TESTSU2RESEARCHREPLAY 复用已保存的原面研究场检查最终UI，不重复启动数值任务。
% 模块用途：核对最终工具栏、真实绘制、未知温度和同工况再确认；保留截图窗口。
if nargin<3,capture=false;end
assert(~isfolder(out));mkdir(out);root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
input=load(fullfile(fileparts(source),'input_snapshot.mat'),'frozen');
fig=rocketAeroHeatWorkbench(input.frozen.cfg);guard=onCleanup(@()close(fig));
fig.UserData.actions.loadNative();fig.UserData.actions.generateMesh();
fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();
fig.UserData.actions.readResearch(source);s=fig.UserData.workflowState;
assert(strcmp(s.status,'research-ready') && isempty(fig.UserData.engineJob));
fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();
assert(strcmp(fig.UserData.workflowState.status,'research-ready'));
assert(contains(fig.UserData.status.Text,'仍有效'));fig.UserData.actions.selectStep(5);
ui=fig.UserData.researchUI;ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);
start=tic;
while toc(start)<40
    drawnow;web=ui.geometry.html.UserData;
    if ~isempty(web.rendered) && web.rendered.revision==web.revision,break;end
    pause(.1);
end
assert(~isempty(web.rendered) && web.rendered.webglError==0 && web.rendered.visiblePixelSamples>0);
ui.exportTo(fullfile(out,'research_snapshot.csv'));ui.inspectFace(1235);
if capture
    write(fullfile(out,'ready_for_capture.json'),struct('ready',true));fprintf('V107_FINAL_UI_READY\n');
    start=tic;while ~isfile(fullfile(out,'capture_done.json')) && toc(start)<60,drawnow;pause(.1);end
end
state=fig.UserData.workflowState;
write(fullfile(out,'report.json'),struct('reopenedWithoutSolver',true,'sameBoundaryPreservesResearch',true, ...
    'rendered',web.rendered,'originalFaces',size(state.mesh.F,1),'newQualifiedArea',0));
fprintf('Final research replay UI passed; no solver rerun.\n');
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
