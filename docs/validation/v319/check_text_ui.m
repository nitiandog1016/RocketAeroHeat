function check_text_ui(out)
% 通过普通文件页加载十进制米制坐标，用用户可操作的 +90° 绕 X 旋转确认姿态。
root=fileparts(which('rocketAeroHeatWorkbench'));path=fullfile(out,'translated-m.stl');
check_identity(out); % 同目录实际写出文本试件，同时检查七种拒绝路径。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.file.Value=path;c.scale.Value=1;press(c.load);
c=fig.UserData.controls;c.rx.Value=90;press(c.rotate);press(fig.UserData.controls.confirm);
c=fig.UserData.controls;c.meshBackend.Value='thermistorB';c.meshBackend.ValueChangedFcn(c.meshBackend,[]);
press(fig.UserData.controls.generate);
s=fig.UserData.workflowState;assert(~isempty(s.mesh),'%s',fig.UserData.status.Text);
ids=readmatrix(fullfile(root,'outputs','v319_validation','fixture','reference_ids.txt'));
b=s.mesh.thermistorImport;
assert(isequal(b.referenceFaceID,ids) && b.maxVertexErrorM<1e-16 && b.maxRelativeAreaError<1e-11);
shown=fig.UserData.geometryView.html.UserData.payload;
assert(isequal(shown.vertices,s.mesh.thermistorPrepared.preview.V+b.translationM));
rejected=false;
try
    aeroHeatWorkflowState('setMesh',s,struct('mesh',rmfield(s.mesh,'thermistorImport')));
catch ex
    rejected=strcmp(ex.identifier,'RAH:Bead:ImportMesh');
end
assert(rejected,'普通 STL 入口接收了缺少源文件绑定的体网。');
press(fig.UserData.controls.nextBoundary);press(fig.UserData.controls.confirm);press(fig.UserData.controls.readBead);
s=fig.UserData.workflowState;assert(s.currentStep==5,'%s',fig.UserData.status.Text);
assert(isequal(s.result.F,s.orientedMesh.F) && isequal(s.result.V,s.orientedMesh.V));
% 原面屏幕拾取及截线在用户坐标中完成，记录平移后的真实截面位置。
v=fig.UserData.webUI.evidenceView;v.selectField('pressure');id=s.result.thermistorMeasuredPoint.peakFaceID;
v.inspectFace(id,[]);u=v.getState();assert(u.selectedFaceID==id);
x=b.translationM(1)+.000172;v.setSectionX(x);timer=tic;
while toc(timer)<40
    drawnow;cut=v.getSection();
    if isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'pressure')&&abs(cut.x-x)<1e-13,break;end
    pause(.1);
end
assert(isstruct(cut)&&isfield(cut,'csv')&&strcmp(cut.fieldKey,'pressure')&&abs(cut.x-x)<1e-13);
v.exportSectionTo(fullfile(out,'text-section.csv'));
g=fig.UserData.webUI.html.UserData.rendered;assert(g.faces==numel(ids) && g.visiblePixelSamples>0 && g.webglError==0);
proof=struct('faces',numel(ids),'sha256',b.source.sha256,'scale',b.source.scale,'peak',id,'referencePeak',ids(id), ...
    'xM',x,'translationM',b.translationM,'rotationMatrix',b.rotationMatrix,'areaM2',sum(b.areaM2), ...
    'maxVertexErrorM',b.maxVertexErrorM,'maxAreaRelativeError',b.maxRelativeAreaError,'visiblePixels',g.visiblePixelSamples);
fid=fopen(fullfile(out,'text-ui.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
