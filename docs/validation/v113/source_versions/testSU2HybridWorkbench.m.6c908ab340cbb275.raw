function report=testSU2HybridWorkbench(out,capture)
%TESTSU2HYBRIDWORKBENCH 实际主工作台创建突起、混合建网格、运行和导出完整原面研究场。
% 模块一：已确认表面转内嵌的原面身份及旧模式配置回归。
% 模块二：带整流罩完整模型走原五步骤，两步真实NS只作接口验证。
% 模块三：实际WebGL、原面拾取及周向曲线与给定/未知温度语义核对。
if nargin<2,capture=false;end
assert(~isfolder(out));mkdir(out);evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.flight.Mach=5;cfg.flight.alpha=20;cfg.flight.beta=0;cfg.flight.Twall=300;
cfg.su2=struct('pythonExecutable','D:/anacondapython/python.exe','iterations',2,'maximumSeconds',120,'outputEvery',2);
[~,old]=su2WorkbenchMeshOptions(cfg);assert(strcmp(old,'starRay'));
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()close(fig));
% 真实19752面先验证接入动作原样保留，不重复生成其已通过的体网格。
d=load('outputs/v113_validation/actual_input.mat','V','F');fig.UserData.actions.useSurfaceForSU2(struct('V',d.V,'F',d.F));
s=fig.UserData.workflowState;assert(isequal(s.orientedMesh.V,d.V) && isequal(s.orientedMesh.F,d.F));
assert(strcmp(fig.UserData.controls.volumeMode.Value,'hybridTriangles'));
fig.UserData.actions.loadNativeFeature();before=fig.UserData.workflowState;
assert(~isempty(before.modelSource.truth.features));
fig.UserData.actions.generateMesh();state=fig.UserData.workflowState;
assert(~isempty(state.mesh),fig.UserData.status.Text);
assert(isequal(state.mesh.V,before.orientedMesh.V) && isequal(state.mesh.F,before.orientedMesh.F));
volumeInfo=fig.UserData.volumeInfo;assert(volumeInfo.construction.interfaceExactlyPreserved);
fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();
fig.UserData.actions.run();assert(fig.UserData.workflowState.busy,fig.UserData.status.Text);
t=tic;while fig.UserData.workflowState.busy && toc(t)<160,drawnow;pause(.1);end
state=fig.UserData.workflowState;assert(strcmp(state.status,'research-ready'),fig.UserData.status.Text);
r=state.researchResult;assert(numel(r.pressureRaw)==size(before.orientedMesh.F,1) && ~r.productionQualified);
ui=fig.UserData.researchUI;ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);
ui.setSectionX(before.modelSource.truth.features(1).centerX);
t=tic;while toc(t)<40
    drawnow;web=ui.geometry.html.UserData;
    if ~isempty(web.rendered) && web.rendered.revision==web.revision && ~isempty(ui.getSection()),break;end
    pause(.1);
end
assert(~isempty(web.rendered) && web.rendered.webglError==0 && web.rendered.visiblePixelSamples>0);
ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'section.csv'));
featureFace=find(before.modelSource.truth.faceFeatureMask(:,1),1);ui.inspectFace(featureFace);
report=struct('actualOriginalSurfaceRetained',true,'featureFace',featureFace, ...
    'originalFaces',size(r.F,1),'maximumProtrusionHeight',max(before.modelSource.truth.nodeRadialOffset), ...
    'volumeInfo',volumeInfo,'engineFolder',r.resultDirectory,'rendered',web.rendered, ...
    'newQualifiedArea',0,'physicsQualified',false,'flowIterations',2);
save(fullfile(out,'feature_geometry.mat'),'before','report','-v7');
write(fullfile(out,'report.json'),report);
if capture
    write(fullfile(out,'ready_for_capture.json'),struct('engineFolder',r.resultDirectory));fprintf('V113_HYBRID_UI_READY\n');
    t=tic;while ~isfile(fullfile(out,'capture_done.json')) && toc(t)<60,drawnow;pause(.1);end
end
% 上游模式改变应撤销网格和下游研究结果。
fig.UserData.actions.selectStep(2);c=fig.UserData.controls.volumeMode;c.Value='normalLayers';c.ValueChangedFcn(c,[]);
assert(isempty(fig.UserData.workflowState.researchResult) && isempty(fig.UserData.workflowState.mesh));
fprintf('Actual protrusion main workbench, original faces, mesh/solve/result/section/export and invalidation passed.\n');
end
function write(path,data)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
