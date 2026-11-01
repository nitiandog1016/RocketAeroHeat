function report=testSU2CreaseFanWorkbench(out)
%TESTSU2CREASEFANWORKBENCH 真实主工作台核对原面混合网格与棱柱扇截面，不启动流场。
% 模块一：读取已验证几何，独立核对完整体积与原面；保留质量失败状态。
% 模块二：在同一结果区切换整箭与体网格截面，验证真实来源和上游失效。
assert(~isfolder(out));mkdir(out);d=load('outputs/v118_validation/actual01/surface.mat');
c=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=c.cfg;cfg.output.verbose=false;
cfg.su2=struct('volumeMode','hybridTriangles','iterations',2,'outputEvery',2,'maximumSeconds',120,'pythonExecutable','D:/anacondapython/python.exe');
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));
fig.UserData.actions.useSurfaceForSU2(d.s);
fig.UserData.actions.readHybridMesh('outputs/v133_validation/fan_workbench');
state=fig.UserData.workflowState;vi=fig.UserData.volumeInfo;
assert(isequal(state.mesh.V,d.s.V) && isequal(state.mesh.F,d.s.F));
assert(vi.nodes==350632 && vi.cells==799494 && vi.construction.tetsBelowSICN01==17);
assert(strcmp(vi.construction.method,'original-triangle-crease-fan-and-gmsh-exterior-tetrahedra'));
assert(contains(fig.UserData.controls.meshMethod.Text,'未获物理资格'));
assert(strcmp(fig.UserData.controls.generate.Text,'重新生成常规体网格'));
waitGeometry(fig,19752);
button=fig.UserData.controls.creasePreview;button.ButtonPushedFcn(button,[]);waitGeometry(fig,62);
info=fig.UserData.creaseSlice;assert(info.geometryOnly && info.sourceVertexID==18 && info.faces==62);
payload=fig.UserData.geometryView.html.UserData.payload;
assert(contains(payload.caption,'仅几何') && isempty(state.researchResult));
save(fullfile(out,'slice.mat'),'info','payload','-v7');
button=fig.UserData.controls.creasePreview;button.ButtonPushedFcn(button,[]);waitGeometry(fig,19752);
assert(isempty(fig.UserData.creaseSlice));
fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();
assert(fig.UserData.workflowState.boundaryConfirmed && ~fig.UserData.workflowState.busy);
fig.UserData.actions.selectStep(2);ctrl=fig.UserData.controls.volumeMode;ctrl.Value='normalLayers';ctrl.ValueChangedFcn(ctrl,[]);
assert(isempty(fig.UserData.workflowState.mesh) && isempty(fig.UserData.workflowState.researchResult));
report=struct('originalFaces',19752,'nodes',350632,'cells',799494,'sliceFaces',62,'sliceVertexID',18, ...
    'mainWorkbenchReadOnlyAccepted',true,'webglPassed',true,'originalSurfaceUnchanged',true, ...
    'upstreamInvalidationPassed',true,'qualityNotClaimed',true,'nativeFlowRuns',0,'physicsQualified',false,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
end
function waitGeometry(fig,faces)
t=tic;while toc(t)<40
    drawnow;s=fig.UserData.geometryView.html.UserData;
    if ~isempty(s.rendered) && s.rendered.revision==s.revision && s.rendered.faces==faces,break;end
    pause(.1);
end
assert(~isempty(s.rendered) && s.rendered.revision==s.revision && s.rendered.faces==faces && s.rendered.webglError==0);
end
function closeIfPresent(fig)
if isgraphics(fig),close(fig);end
end
