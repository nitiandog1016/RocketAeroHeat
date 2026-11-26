function testWorkflowWorkbench(out)
%TESTWORKFLOWWORKBENCH 实际五步骤路径及上游修改使下游失效的集成校核。
% 模块用途：用真实STL/重构网格/当前计算和实际WebGL检验单窗口操作。
% 浏览不计算；无效参数、失败输入或求解不接受旧结果。不是实验验证。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
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
assert(s.currentStep==4 && strcmp(s.error.identifier,'RocketAeroHeat:WorkflowIncompleteSurface'),fig.UserData.status.Text);
assert(s.boundaryConfirmed && ~s.busy && ~s.steps(5).available && isempty(s.result));
assert(strcmp(fig.UserData.workflowButtons(5).Enable,'off') && isempty(fig.UserData.webUI));
assert(~isfield(fig.UserData.controls,'export') && isempty(fig.UserData.caseReport));
assert(strcmp(fig.UserData.editorPanel.Enable,'on') && strcmp(fig.UserData.controls.run.Enable,'on'));

%% 3. 查看设置不重算，编辑边界后需要用户再确认。
solveId=s.solveId;press(fig.UserData.workflowButtons(3));press(fig.UserData.workflowButtons(4));
s=fig.UserData.workflowState;assert(s.solveId==solveId && s.revision==before && isempty(s.result));
press(fig.UserData.workflowButtons(3));field=fig.UserData.controls.alpha;
field.Value=10;field.ValueChangedFcn(field,[]);s=fig.UserData.workflowState;
assert(isempty(s.result) && ~s.boundaryConfirmed && ~s.steps(5).available);
press(fig.UserData.controls.confirm);press(fig.UserData.controls.run);s=fig.UserData.workflowState;
assert(s.solveId>solveId && strcmp(s.error.identifier,'RocketAeroHeat:WorkflowIncompleteSurface'));

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
assert(isempty(s.result) && strcmp(s.error.identifier,'RocketAeroHeat:WorkflowIncompleteSurface'));
assert(~s.steps(5).available && isempty(fig.UserData.webUI));

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
if nargin>0
    if ~isfolder(out),mkdir(out);end
    save(fullfile(out,'workflow_workbench.mat'),'meshFaces','failedReason');
end
fprintf('工作台复核通过：%d 面，缺项工况停在计算页，结果和导出关闭，边界可重新调整。\n',meshFaces);
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
function closeFigure(fig)
if isgraphics(fig),delete(fig);end
end
