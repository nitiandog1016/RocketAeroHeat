function testWorkflowWorkbench(out)
%TESTWORKFLOWWORKBENCH 实际五步骤路径及上游修改使下游失效的集成校核。
% 模块用途：用真实STL/重构网格/当前计算和实际WebGL检验单窗口操作。
% 浏览不计算；无效参数、失败输入或求解不接受旧结果。不是实验验证。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if nargin<1,out=fullfile(root,'outputs','v56_validation','workflow_final');end
if ~exist(out,'dir'),mkdir(out);end
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;cfg.flight.alpha=20;
fig=rocketAeroHeatWorkbench(cfg);cleanup=onCleanup(@()closeFigure(fig)); %#ok<NASGU>
initial=fig.UserData.workflowState;
assert(numel(fig.UserData.workflowButtons)==5 && initial.currentStep==1);
assert(~any([initial.steps(2:5).available]) && isempty(initial.result));
geometry=fig.UserData.geometryView;
assert(~isfield(geometry.html.UserData.payload,'fields'));

%% 1. 原始STL输入、真实PCA朝向确认，生成同一算法的求解网格
press(fig.UserData.controls.load);s=fig.UserData.workflowState;
assert(size(s.rawMesh.F,1)==3176 && isempty(s.mesh) && isempty(s.orientedMesh));
press(fig.UserData.controls.pca);press(fig.UserData.controls.confirm);
s=fig.UserData.workflowState;assert(s.currentStep==2 && s.steps(1).completed);
assert(s.cfg.geometry.isAligned && isempty(s.result));
press(fig.UserData.controls.generate);s=fig.UserData.workflowState;
assert(~isempty(s.mesh) && s.steps(2).completed,fig.UserData.status.Text);
assert(~isfield(s.mesh,'heating') && isempty(s.result) && ~s.boundaryConfirmed);
meshFaces=size(s.mesh.F,1);assert(meshFaces>3000);
waitGeometry(fig.UserData.geometryView,size(s.mesh.F,1));

%% 2. 边界确认只求自由流；明确按钮才建立压力/热流和完整当前快照
press(fig.UserData.workflowButtons(3));s=fig.UserData.workflowState;
scopeRevision=s.revision;scopeSolveId=s.solveId;
press(fig.UserData.controls.methodToggle);s=fig.UserData.workflowState;
assert(s.revision==scopeRevision && s.solveId==scopeSolveId && isempty(s.result));
description=strjoin(fig.UserData.controls.methodScope.Value,newline);
assert(contains(description,'不是程序求得尾迹压力为零') && contains(description,'不求解壁体温度'));
assert(contains(description,'只有有效外缘热力状态才能用于Kn流态判断'));
% 真实选择回调给出与当前20度工况冲突的解释，且不启动昂贵鼻区BVP。
field=fig.UserData.controls.nose;field.Value='sphericalLocalSimilarity';field.ValueChangedFcn(field,[]);
description=strjoin(fig.UserData.controls.methodScope.Value,newline);
assert(contains(fig.UserData.controls.methodSummary.Text,'球鼻轴对称局部相似热边界不适用'));
assert(contains(description,'δ/R≤0.05') && contains(description,'鼻部/后体接口'));
s=fig.UserData.workflowState;assert(s.solveId==scopeSolveId && isempty(s.result) && ~isfield(s.mesh,'heating'));
field.Value='cosineExtension';field.ValueChangedFcn(field,[]);
press(fig.UserData.controls.confirm);
s=fig.UserData.workflowState;assert(s.currentStep==4 && s.boundaryConfirmed && isempty(s.result));
assert(~isfield(s.mesh,'pressure') && ~isfield(s.mesh,'heating'));
before=s.revision;press(fig.UserData.controls.run);s=fig.UserData.workflowState;
assert(s.currentStep==5 && ~isempty(s.result) && s.resultSourceRevision==before,fig.UserData.status.Text);
assert(s.steps(4).completed && s.steps(5).completed);
assert(size(s.result.F,1)==meshFaces && ~s.busy);
assert(strcmp(fig.UserData.editorPanel.Enable,'on'), ...
    '计算完成后父编辑面板必须恢复，不能只检查子控件Enable属性。');
assert(strcmp(fig.UserData.controls.export.Enable,'on'));
result=s.result;report=fig.UserData.caseReport;
assert(report.faceCount==meshFaces && report.heatComputedAreaFraction>0);
ui=fig.UserData.webUI;waitResult(ui,meshFaces);
assert(ancestor(ui.html,'figure')==fig && ancestor(geometry.html,'figure')==fig);
data=buildAerothermalViewerData(result);
wall=data.fields(strcmp({data.fields.key},'wallInput'));solved=data.fields(strcmp({data.fields.key},'wallSolved'));
assert(all(wall.status==2) && all(solved.status==0) && all(isnan(solved.values)));
% 检查真实20度工况显示：域外流向占位不可染成已求外缘场，未知上游不可输出热流。
placeholder=result.boundaryLayerEdge.diagnosticPlaceholderFace;
assert(any(placeholder),'本真实工况应覆盖域外热力占位的显示路径。');
for key={'edgePressure','edgeTemperature'}
    fieldData=data.fields(strcmp({data.fields.key},key{1}));
    assert(all(fieldData.status(placeholder)==0) && all(isnan(fieldData.values(placeholder))));
end
blocked=result.continuumRegime.stateUnresolvedFace;
assert(any(blocked) && all(~result.heating.eligibleFace(blocked)));
heat=data.fields(strcmp({data.fields.key},'heatFlux'));
assert(all(heat.status(blocked)==0) && all(isnan(heat.values(blocked))));
assert(abs(heat.computedAreaFraction-report.heatComputedAreaFraction)<1e-12);
assert(contains(fig.UserData.controls.export.Text,'导出'));
geometryRender=geometry.html.UserData.rendered;resultRender=ui.html.UserData.rendered;

%% 3. 查看已完成设置不触发求解或改变结果；编辑边界使旧场量立即失效
solveId=s.solveId;press(fig.UserData.workflowButtons(3));press(fig.UserData.workflowButtons(5));
s=fig.UserData.workflowState;assert(s.solveId==solveId && s.resultSourceRevision==before);
assert(isequaln(s.result.heating.heatFlux,result.heating.heatFlux));
press(fig.UserData.workflowButtons(3));press(fig.UserData.controls.confirm);s=fig.UserData.workflowState;
assert(s.solveId==solveId && s.currentStep==4 && ~isempty(s.result));
assert(contains(fig.UserData.status.Text,'结果仍有效'), ...
    '未改变输入的确认不得把有效快照错误描述为尚未计算。');
press(fig.UserData.workflowButtons(3));field=fig.UserData.controls.alpha;
field.Value=10;field.ValueChangedFcn(field,[]);s=fig.UserData.workflowState;
assert(isempty(s.result) && ~s.boundaryConfirmed && ~s.steps(5).available);
assert(~isfield(s.mesh,'heating') && strcmp(fig.UserData.workflowButtons(5).Enable,'off'));
press(fig.UserData.controls.confirm);press(fig.UserData.controls.run);s=fig.UserData.workflowState;
assert(s.result.pressure.freestream.alpha==10 && s.solveId>solveId && s.resultSourceRevision==s.revision);
assert(~isequaln(s.result.pressure.staticPressure,result.pressure.staticPressure));

%% 4. 真实失败计算留在同宿主，无过期结果；用户可以修改条件重试
press(fig.UserData.workflowButtons(3));field=fig.UserData.controls.base;
field.Value='sand93Laminar';field.ValueChangedFcn(field,[]);
description=strjoin(fig.UserData.controls.methodScope.Value,newline);
assert(contains(description,'半锥角9–10°') && contains(description,'球鼻、锥柱、尾锥不是这个几何族'));
assert(contains(description,'当前α 10°') && contains(description,'当前变比热物性不满足'));
assert(strcmp(fig.UserData.controls.baseCp.Enable,'off'), ...
    '选择关联时不能把未被读取的给定底Cp控件展示为有效求解参数。');
press(fig.UserData.controls.confirm);press(fig.UserData.controls.run);s=fig.UserData.workflowState;
assert(isempty(s.result) && strcmp(s.status,'failed') && ~s.busy);
assert(~s.steps(5).available && strcmp(fig.UserData.controls.run.Enable,'on'));
assert(strcmp(fig.UserData.editorPanel.Enable,'on'));
assert(~isempty(s.error.message));
failureText=strjoin(fig.UserData.controls.failure.Value,newline);
assert(contains(failureText,s.error.message) && contains(failureText,s.error.identifier));
failedRevision=s.revision;failedSolveId=s.solveId;failedReason=s.error.message;
press(fig.UserData.controls.methodToggle);s=fig.UserData.workflowState;
assert(s.revision==failedRevision && s.solveId==failedSolveId && strcmp(s.error.message,failedReason));
assert(contains(strjoin(fig.UserData.controls.failure.Value,newline),failedReason));
press(fig.UserData.controls.adjustBoundary);assert(fig.UserData.workflowState.currentStep==3);
press(fig.UserData.workflowButtons(4));s=fig.UserData.workflowState;
assert(s.revision==failedRevision && s.solveId==failedSolveId && isempty(s.result));
assert(contains(strjoin(fig.UserData.controls.failure.Value,newline),failedReason));
press(fig.UserData.controls.adjustBoundary);field=fig.UserData.controls.base;
field.Value='prescribed';field.ValueChangedFcn(field,[]);
assert(strcmp(fig.UserData.controls.baseCp.Enable,'on'));
press(fig.UserData.controls.confirm);press(fig.UserData.controls.run);s=fig.UserData.workflowState;
assert(~isempty(s.result) && strcmp(s.status,'solved'));

%% 5. 修改网格、实际朝向及模型尺度，分别证明下游清空并可继续操作
press(fig.UserData.workflowButtons(2));field=fig.UserData.controls.hMin;
field.Value=field.Value+1;field.ValueChangedFcn(field,[]);s=fig.UserData.workflowState;
assert(isempty(s.mesh) && isempty(s.result) && ~s.steps(3).available && ~isempty(s.orientedMesh));
press(fig.UserData.controls.generate);s=fig.UserData.workflowState;assert(~isempty(s.mesh));
press(fig.UserData.workflowButtons(1));press(fig.UserData.controls.flip);s=fig.UserData.workflowState;
assert(isempty(s.orientedMesh) && isempty(s.mesh) && isempty(s.result) && ~isempty(s.rawMesh));
press(fig.UserData.controls.confirm);s=fig.UserData.workflowState;assert(~isempty(s.orientedMesh));
press(fig.UserData.workflowButtons(1));field=fig.UserData.controls.scale;
field.Value=field.Value*2;field.ValueChangedFcn(field,[]);s=fig.UserData.workflowState;
assert(isempty(s.rawMesh) && isempty(s.orientedMesh) && isempty(s.mesh) && isempty(s.result));
assert(~s.steps(2).available && s.cfg.geometry.scale==2*cfg.geometry.scale);
% 输入错误不会保留先前模型的物理结果或禁用用户修复控件。
field=fig.UserData.controls.file;field.Value='missing-workflow-test.stl';field.ValueChangedFcn(field,[]);
press(fig.UserData.controls.load);s=fig.UserData.workflowState;
assert(isempty(s.rawMesh) && isempty(s.result) && strcmp(fig.UserData.controls.load.Enable,'on'));
exportAerothermalWeb(result,fullfile(out,'workflow_actual_rocket'));
save(fullfile(out,'workflow_workbench.mat'),'result','report','cfg','meshFaces','data','geometryRender','resultRender');
fprintf('Fresh actual 20-degree snapshot: heat computed area %.10f percent; diagnostic placeholder faces %d; unknown-state/downstream blocked faces %d; geometry GL %d and result GL %d.\n', ...
    100*report.heatComputedAreaFraction,nnz(placeholder),nnz(blocked),geometryRender.webglError,resultRender.webglError);
fprintf('Workflow workbench passed: 5 actual sidebar steps, %d true mesh faces, explicit solve, live WebGL, browsing/scope expansion without recomputation; real nose/base selection limitations, prescribed-field identity, original failure reasons preserved across pages, parent reenabled, retry and upstream invalidation in one window.\n',meshFaces);
end
function press(control)
control.ButtonPushedFcn(control,[]);drawnow;
end
function waitGeometry(view,n)
start=tic;
while toc(start)<40
    drawnow;r=view.html.UserData.rendered;
    if ~isempty(r) && r.faces==n && r.revision==view.html.UserData.revision,break;end
    pause(.1);
end
r=view.html.UserData.rendered;
assert(view.html.UserData.ready && ~isempty(r) && r.faces==n && r.webglError==0 && r.visiblePixelSamples>0);
end
function waitResult(ui,n)
start=tic;
while (~ui.html.UserData.ready || isempty(ui.html.UserData.rendered)) && toc(start)<40,drawnow;pause(.1);end
r=ui.html.UserData.rendered;
assert(ui.html.UserData.ready && ~isempty(r) && r.faces==n && r.webglError==0 && r.visiblePixelSamples>0);
end
function closeFigure(fig)
if isgraphics(fig),delete(fig);end
end
