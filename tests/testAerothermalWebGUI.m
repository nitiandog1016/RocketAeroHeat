function testAerothermalWebGUI(out,sourceFile)
%TESTAEROTHERMALWEBGUI 用导入箭体检查离线导出、默认页及 MATLAB/网页握手。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if nargin<1,out=fullfile(tempdir,'RocketAeroHeat_web_test');end
if ~exist(out,'dir'),mkdir(out);end
% 可复用已保存的真实箭体计算状态，专门复测网页修复，避免重复外流求解。
if nargin>=2
    saved=load(sourceFile,'mesh','cfg');mesh=saved.mesh;cfg=saved.cfg;
else
    evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
    cfg.flight.alpha=20;
    raw=calcFaceGeometry(loadSTLMesh(cfg));[raw,~]=alignGeometryAxis(raw,'pca');
    evalc('mesh=autoRemeshSurface(raw,cfg);mesh=fixNormalDirection(mesh,cfg);');
    mesh=buildSurfaceCoordinates(mesh,cfg);mesh=analyzeLocalGeometry(mesh,cfg);
    mesh=calcFlightStationState(mesh,cfg.flight,cfg);
end
data=buildAerothermalViewerData(mesh);
assert(size(data.faces,1)==size(mesh.F,1));
path=exportAerothermalWeb(mesh,fullfile(out,'actual_rocket'));
fid=fopen(fullfile(out,'actual_rocket','results.json'),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(data,'ConvertInfAndNaN',true));fclose(fid);
assert(isfile(path));
assert(isfile(fullfile(out,'actual_rocket','export.js')), ...
    '离线包缺少当前场量导出资源。');
% 复测读取的基线文件保持原样；首次端到端计算才保存供后续使用。
if nargin<2,save(fullfile(out,'actual_rocket.mat'),'mesh','cfg');end
fig=pressureGUI(mesh,cfg);clean=onCleanup(@()deleteIfValid(fig)); %#ok<NASGU>
drawnow;ui=fig.UserData.webUI;
assert(fig.UserData.tabs.SelectedTab==ui.tab);
assert(strcmp(fig.UserData.selector.Enable,'off'));
assert(strcmp(fig.UserData.selector.Visible,'off'));
assert(strcmp(fig.UserData.workbenchLabel.Visible,'on'));
assert(strcmp(fig.UserData.loadToggle.Enable,'on'));
assert(strcmp(fig.UserData.changeFlight.Visible,'on'));
assert(all(cellfun(@(x)isequal(x,0),fig.UserData.resultToolbar.ColumnWidth(4:9))));
assert(strcmp(fig.UserData.basePressureButton.Visible,'on'));
assert(strcmp(fig.UserData.noseSelector.Visible,'on'));
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
assert(strcmp(fig.UserData.selector.Visible,'on'));
assert(strcmp(fig.UserData.workbenchLabel.Visible,'off'));
assert(strcmp(fig.UserData.noseSelector.Visible,'off'));
fig.UserData.tabs.SelectedTab=ui.tab;
fig.UserData.tabs.SelectionChangedFcn([],[]);drawnow;pause(2);
assert(strcmp(fig.UserData.selector.Visible,'off'));
% 气动力/参考量可在工作台直接展开；不要求先切换旧版压力页。
fig.UserData.loadToggle.ButtonPushedFcn([],[]);drawnow;
assert(strcmp(fig.UserData.loadsPanel.Visible,'on'));
assert(fig.UserData.tabs.SelectedTab==ui.tab);
fig.UserData.loadToggle.ButtonPushedFcn([],[]);drawnow;
assert(strcmp(fig.UserData.loadsPanel.Visible,'off'));
% exportapp 在当前 MATLAB 的 uihtml 整窗导出不受支持；用 GPU 回执验证内嵌绘制。
% 浏览器截图与交互另行验证，不能把握手通过当作绘制通过。
fprintf('Web GUI passed: %d actual faces, uihtml handshake, default tab, update, legacy controls.\n',size(mesh.F,1));
end
function deleteIfValid(fig)
if isgraphics(fig),delete(fig);end
end
