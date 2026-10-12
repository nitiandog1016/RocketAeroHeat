function testAerothermalWebGUI(out)
%TESTAEROTHERMALWEBGUI 用导入箭体检查离线导出、默认页及 MATLAB/网页握手。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if nargin<1,out=fullfile(tempdir,'RocketAeroHeat_web_test');end
if ~exist(out,'dir'),mkdir(out);end
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.flight.alpha=20;
raw=calcFaceGeometry(loadSTLMesh(cfg));[raw,~]=alignGeometryAxis(raw,'pca');
evalc('mesh=autoRemeshSurface(raw,cfg);mesh=fixNormalDirection(mesh,cfg);');
mesh=buildSurfaceCoordinates(mesh,cfg);mesh=analyzeLocalGeometry(mesh,cfg);
mesh=calcFlightStationState(mesh,cfg.flight,cfg);
data=buildAerothermalViewerData(mesh);
assert(size(data.faces,1)==size(mesh.F,1));
path=exportAerothermalWeb(mesh,fullfile(out,'actual_rocket'));
fid=fopen(fullfile(out,'actual_rocket','results.json'),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(data,'ConvertInfAndNaN',true));fclose(fid);
assert(isfile(path));
save(fullfile(out,'actual_rocket.mat'),'mesh','cfg');
fig=pressureGUI(mesh,cfg);clean=onCleanup(@()deleteIfValid(fig)); %#ok<NASGU>
drawnow;ui=fig.UserData.webUI;
assert(fig.UserData.tabs.SelectedTab==ui.tab);
assert(strcmp(fig.UserData.selector.Enable,'off'));
started=tic;
while (isempty(ui.html.UserData.rendered) || ~ui.html.UserData.ready) && toc(started)<40,drawnow;pause(.1);end
assert(ui.html.UserData.ready,'网页未完成真实 uihtml 握手。');
gpu=ui.html.UserData.rendered;
assert(~isempty(gpu) && gpu.webglError==0 && gpu.visiblePixelSamples>0, ...
    '内嵌 WebGL 未确认实际绘制有效像素。');
assert(gpu.faces==size(mesh.F,1));
fprintf('Embedded GPU evidence: %d faces, %d visible pixel samples, %dx%d, GL error %d.\n', ...
    gpu.faces,gpu.visiblePixelSamples,gpu.width,gpu.height,gpu.webglError);
% 更新新计算值，再切回专题页，检查数据同步及统一窗口保留。
updated=mesh;updated.pressure.staticPressure=mesh.pressure.staticPressure+123;
updateAerothermalWebTab(ui,updated);drawnow;
payload=ui.html.UserData.payload;
p=payload.fields(strcmp({payload.fields.key},'pressure'));
assert(max(abs(p.values-updated.pressure.staticPressure/1000))<1e-10);
started=tic;
while (~isfield(ui.html.UserData.rendered,'pressureSample') || ...
        abs(ui.html.UserData.rendered.pressureSample-p.values(1))>1e-10) && toc(started)<20
    drawnow;pause(.1);
end
assert(abs(ui.html.UserData.rendered.pressureSample-p.values(1))<1e-10, ...
    '网页实际渲染的数据没有同步新压力。');
updateAerothermalWebTab(ui,mesh);
fig.UserData.tabs.SelectedTab=fig.UserData.mapTab;
fig.UserData.tabs.SelectionChangedFcn([],[]);
assert(strcmp(fig.UserData.selector.Enable,'on'));
fig.UserData.tabs.SelectedTab=ui.tab;
fig.UserData.tabs.SelectionChangedFcn([],[]);drawnow;pause(2);
% exportapp 在当前 MATLAB 的 uihtml 整窗导出不受支持；用 GPU 回执验证内嵌绘制。
% 浏览器截图与交互另行验证，不能把握手通过当作绘制通过。
fprintf('Web GUI passed: %d actual faces, uihtml handshake, default tab, update, legacy controls.\n',size(mesh.F,1));
end
function deleteIfValid(fig)
if isgraphics(fig),delete(fig);end
end
