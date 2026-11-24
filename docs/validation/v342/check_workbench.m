function check_workbench(out)
% 模块一：从工作台的算例按钮读取几何，独立核对各区的文件面号以及真实 WebGL 回执。
fig=rocketAeroHeatWorkbench();cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.rfzModel);
s=fig.UserData.workflowState;assert(isfield(s.rawMesh,'rfzRegions'),'%s',fig.UserData.status.Text);
mesh=s.rawMesh;regions=mesh.rfzRegions;html=fig.UserData.geometryView.html;
gl=rendered(html);assert(gl.faces==size(mesh.F,1) && gl.vertices==size(mesh.V,1));
assert(isempty(s.result) && isempty(s.mesh) && isempty(s.freestream) && ~s.boundaryConfirmed);
assert(html.UserData.payload.enableFaceProbe && ~isfield(html.UserData.payload,'faceField'));
for j=[1 13 204]
    face=regions.lastFace(j);notifyFace(html,face,html.UserData.revision);
    assert(contains(fig.UserData.status.Text,sprintf('原面 %d · CAD 面 %d ·',face,j)));
end
last=fig.UserData.status.Text;notifyFace(html,1,html.UserData.revision-1);
assert(strcmp(fig.UserData.status.Text,last));
press(fig.UserData.controls.confirm);s=fig.UserData.workflowState;
assert(s.currentStep==2 && isequal(s.orientedMesh.F,mesh.F) && isequal(s.orientedMesh.V,mesh.V) && isempty(s.result));
glOriented=rendered(html);

% 模块二：用同一窗口的普通文件入口重开，并验证尺度变化使旧编号资格失效。
press(fig.UserData.workflowButtons(1));c=fig.UserData.controls;file=c.file.Value;
c.scale.Value=.001;press(c.load);s=fig.UserData.workflowState;
assert(~isfield(s.rawMesh,'rfzRegions') && isempty(s.orientedMesh) && isempty(s.result));
c=fig.UserData.controls;c.file.Value=file;c.scale.Value=1;press(c.load);
s=fig.UserData.workflowState;assert(isequal(s.rawMesh.rfzRegions,regions) && isequal(s.rawMesh.F,mesh.F));
glReopen=rendered(html);
record=struct('stlSHA256',s.modelSource.sha256,'faces',size(mesh.F,1),'vertices',size(mesh.V,1), ...
    'regions',height(regions),'areaM2',sum(mesh.area),'presetPixels',gl.visiblePixelSamples, ...
    'orientationPixels',glOriented.visiblePixelSamples,'fileReopenPixels',glReopen.visiblePixelSamples, ...
    'probeCallbackCADFaces',[1 13 204],'staleProbeIgnored',true,'changedScaleUnbound',true, ...
    'computedFaces',0,'jointExperimentalAreaM2',0);
file=fullfile(out,'workbench.json');assert(~isfile(file));
f=fopen(file,'w','n','UTF-8');assert(f>=0);closer=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',jsonencode(record,'PrettyPrint',true));
end

function msg=rendered(html)
start=tic;
while toc(start)<55
    drawnow;u=html.UserData;
    assert(isempty(u.error),'RAH:RFZ:WebGL','%s',jsonencode(u.error));
    msg=u.rendered;
    if ~isempty(msg) && msg.revision==u.revision && msg.webglError==0 && msg.visiblePixelSamples>0,return;end
    pause(.1);
end
error('RAH:RFZ:WebGL','未收到当前几何的有效绘制回执。');
end

function press(button)
button.ButtonPushedFcn(button,[]);drawnow;
end

function notifyFace(html,face,revision)
% 从网页消息边界测试面号转换，不以程序发送的消息代替鼠标实测。
html.Data=struct('kind','geometry-face-probe','revision',revision,'faceIndex',face-1);
html.DataChangedFcn(html,[]);
end
