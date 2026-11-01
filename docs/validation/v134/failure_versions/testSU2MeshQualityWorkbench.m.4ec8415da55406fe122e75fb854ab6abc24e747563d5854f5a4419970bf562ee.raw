function report=testSU2MeshQualityWorkbench(out)
%TESTSU2MESHQUALITYWORKBENCH 主工作台验证低质量体定位、数值来源与旧帧拒绝，不求流场。
% 模块一：原面和真实体网格保持，读取独立重算SICN的彩色单元及灰色原面。
% 模块二：验证WebGL当前帧、面到体的探针消息、切换及上游失效。
assert(~isfolder(out));mkdir(out);d=load('outputs/v118_validation/actual01/surface.mat');
c=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=c.cfg;cfg.output.verbose=false;
cfg.su2=struct('volumeMode','hybridTriangles','iterations',2,'outputEvery',2,'maximumSeconds',120,'pythonExecutable','D:/anacondapython/python.exe');
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));
fig.UserData.actions.useSurfaceForSU2(d.s);fig.UserData.actions.readHybridMesh('outputs/v133_validation/fan_workbench');
button=fig.UserData.controls.meshQuality;button.ButtonPushedFcn(button,[]);
waitGeometry(fig);info=fig.UserData.meshQuality;payload=fig.UserData.geometryView.html.UserData.payload;
assert(info.lowQualityCells==17 && abs(info.minimumSICN-.0139122574504664)<1e-12);
assert(strcmp(payload.faceField.context,'mesh-quality') && isequal(payload.faceField.range,[0 .1]));
assert(all(isnan(payload.faceField.values(info.visibleCellFaces+1:end))) && numel(info.contextOriginalFaceIDs)>0);
assert(isequal(fig.UserData.workflowState.mesh.V,d.s.V) && isequal(fig.UserData.workflowState.mesh.F,d.s.F));
assert(isempty(fig.UserData.workflowState.researchResult));
probe=struct('kind','geometry-face-probe','revision',payload.revision,'faceIndex',0);
% 使用与网页同一消息接收入口；此处为定向桥接测试，不声称原生鼠标已操作。
html=fig.UserData.geometryView.html;html.Data=probe;html.DataChangedFcn(html,[]);
picked=fig.UserData.meshQualityProbe;assert(picked.cellID==info.displaySourceCellIDs(1));
assert(abs(picked.SICN-info.displayValues(1))<1e-15 && picked.geometryOnly);
old=fig.UserData.controls.meshQualityProbe.Text;probe.revision=payload.revision-1;
fig.UserData.actions.probeMeshQuality(probe);assert(strcmp(old,fig.UserData.controls.meshQualityProbe.Text));
probe.revision=payload.revision;probe.faceIndex=info.visibleCellFaces;fig.UserData.actions.probeMeshQuality(probe);
assert(contains(fig.UserData.controls.meshQualityProbe.Text,'原壁面定位参考'));
save(fullfile(out,'quality_preview.mat'),'info','payload','picked','-v7');
button=fig.UserData.controls.creasePreview;button.ButtonPushedFcn(button,[]);waitGeometry(fig);
assert(isempty(fig.UserData.meshQuality) && fig.UserData.creaseSlice.faces==62);
button=fig.UserData.controls.creasePreview;button.ButtonPushedFcn(button,[]);waitGeometry(fig);
assert(size(fig.UserData.geometryView.html.UserData.payload.faces,1)==19752);
ctrl=fig.UserData.controls.volumeMode;ctrl.Value='normalLayers';ctrl.ValueChangedFcn(ctrl,[]);
assert(isempty(fig.UserData.workflowState.mesh) && isempty(fig.UserData.meshQuality));
report=struct('lowQualityCells',17,'visibleCellFaces',info.visibleCellFaces,'contextFaces',numel(info.contextOriginalFaceIDs), ...
    'minimumSICN',info.minimumSICN,'webglPassed',true,'bridgeProbePassed',true,'staleRevisionRejected',true, ...
    'contextNotMislabelled',true,'sourceSurfaceUnchanged',true,'upstreamInvalidationPassed',true, ...
    'nativeMouseTest',false,'nativeFlowRuns',0,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
end
function waitGeometry(fig)
t=tic;while toc(t)<40
    drawnow;s=fig.UserData.geometryView.html.UserData;
    if ~isempty(s.rendered) && s.rendered.revision==s.revision,break;end
    pause(.1);
end
assert(~isempty(s.rendered) && s.rendered.revision==s.revision && s.rendered.faces==size(s.payload.faces,1) && s.rendered.webglError==0);
end
function closeIfPresent(fig)
if isgraphics(fig),close(fig);end
end
