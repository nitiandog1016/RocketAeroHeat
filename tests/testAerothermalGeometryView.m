function testAerothermalGeometryView(out,sourceFile)
%TESTAEROTHERMALGEOMETRYVIEW 核验真实导入/更新/清空网格与 uihtml 绘制。
% 模块一：选择实际 STL 或已有真实网格快照，不借用任何压力/热流字段。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if nargin<1,out=fullfile(tempdir,'RocketAeroHeat_geometry_view_test');end
if ~exist(out,'dir'),mkdir(out);end
if nargin>=2,saved=load(sourceFile,'mesh');mesh=saved.mesh;
else,evalc('cfg=createCaseConfig();');mesh=loadSTLMesh(cfg);end
mesh=struct('V',mesh.V,'F',mesh.F);
fig=uifigure('Position',[40 40 1150 740],'Name','几何网页真实测试');
clean=onCleanup(@()deleteIfValid(fig)); %#ok<NASGU>
host=uigridlayout(fig,[1 1]);host.Padding=[0 0 0 0];
view=createAerothermalGeometryView(host,[]);waitFor(view,1,40);
assert(view.html.UserData.ready,'空工作台未完成握手。');
assert(view.html.UserData.rendered.faces==0);
assert(view.html.UserData.rendered.visiblePixelSamples==0,'清空场景仍残留几何。');
% 模块二：输入真网格，核对网页当前版本、非零像素和无 WebGL 错误。
updateAerothermalGeometryView(view,mesh,'实际模型输入',struct('showWire',true));
waitFor(view,2,30);state=view.html.UserData;checkGeometry(state,mesh);
assert(~isfield(state.payload,'fields')&&~isfield(state.payload,'pressure')&& ...
    ~isfield(state.payload,'heatFlux'),'几何预览不允许带入伪物理量。');
assert(isequal(state.payload.vertices,mesh.V)&&isequal(state.payload.faces,mesh.F-1));
fprintf('Input GPU: %d faces, %d visible samples, %dx%d, GL %d.\n', ...
    state.rendered.faces,state.rendered.visiblePixelSamples,state.rendered.width, ...
    state.rendered.height,state.rendered.webglError);
% 模块三：改变实际几何及比例，网页应回执新坐标而不是陈旧帧。
updated=mesh;updated.V(:,1)=1.3*updated.V(:,1)+.125;
updateAerothermalGeometryView(view,updated,'更新后的实际模型', ...
    struct('showWire',false,'resetCamera',true));
waitFor(view,3,30);state=view.html.UserData;checkGeometry(state,updated);
payload=state.payload;gpu=state.rendered;save(fullfile(out,'geometry_view_evidence.mat'), ...
    'payload','gpu','mesh','updated');
% 无效网格不得覆盖已确认的当前几何。
bad=updated;bad.F(1,1)=size(updated.V,1)+1;rejected=false;
try,updateAerothermalGeometryView(view,bad,'无效几何');catch err
    rejected=strcmp(err.identifier,'RocketAeroHeat:GeometryView:Mesh');
end
assert(rejected&&view.html.UserData.revision==3);
% 模块四：清空输入必须删除 GPU 中的旧几何，保留空工作台而不是旧结果。
updateAerothermalGeometryView(view,[],'等待输入');waitFor(view,4,30);
state=view.html.UserData;assert(state.rendered.faces==0&& ...
    state.rendered.vertices==0&&state.rendered.visiblePixelSamples==0&& ...
    state.rendered.webglError==0);
fprintf('Geometry view passed: input/update real V/F, no physics fields, empty reset, rejected invalid mesh.\n');
end

function waitFor(view,revision,seconds)
started=tic;
while toc(started)<seconds
    drawnow;state=view.html.UserData;
    if ~isempty(state.error),error('RocketAeroHeat:GeometryView:Browser','%s',state.error.message);end
    if state.ready&&~isempty(state.rendered)&&state.rendered.revision==revision,return;end
    pause(.1);
end
error('RocketAeroHeat:GeometryView:Timeout','未收到几何版本 %d 的真实绘制回执。',revision);
end
function checkGeometry(state,mesh)
gpu=state.rendered;
assert(gpu.faces==size(mesh.F,1)&&gpu.vertices==size(mesh.V,1));
assert(gpu.visiblePixelSamples>0&&gpu.webglError==0&&gpu.width>0&&gpu.height>0);
assert(max(abs(gpu.firstVertex(:)-mesh.V(1,:).'))<1e-12);
expectedBounds=[min(mesh.V,[],1);max(mesh.V,[],1)];
assert(max(abs(gpu.bounds(:)-expectedBounds(:)))<1e-12);
end
function deleteIfValid(fig)
if isgraphics(fig),delete(fig);end
end
