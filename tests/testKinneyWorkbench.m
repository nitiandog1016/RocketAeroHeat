function testKinneyWorkbench(out)
%TESTKINNEYWORKBENCH 验证统一五步骤中的方法选择、当前计算与失效流程。
% 模块一：真实STL、朝向、重网格在同一宿主完成，不创建新结果窗口。
% 模块二：切换侧壁公式只修改边界配置；显式计算才产生新场量。
% 模块三：两种方法返回的缺项结果均须留在计算页，切换方法后可重新操作。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
assert(~exist(out,'dir'));mkdir(out);
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;cfg.flight.alpha=20;
fig=rocketAeroHeatWorkbench(cfg);fig.Visible='off';cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
assert(numel(fig.UserData.workflowButtons)==5);
press(fig.UserData.controls.load);press(fig.UserData.controls.pca);
press(fig.UserData.controls.confirm);press(fig.UserData.controls.generate);
state=fig.UserData.workflowState;assert(state.steps(2).completed && ~isempty(state.mesh));
press(fig.UserData.workflowButtons(3));press(fig.UserData.controls.methodToggle);
controls=fig.UserData.controls;
assert(strcmp(controls.sidewall.Value,'boundaryLayerIntegral'));
before=fig.UserData.workflowState;control=controls.sidewall;
control.Value='kinneyAcreage';control.ValueChangedFcn(control,[]);drawnow;
state=fig.UserData.workflowState;
assert(strcmp(state.cfg.heating.sidewallModel,'kinneyAcreage'));
assert(state.solveId==before.solveId && isempty(state.result) && ~state.boundaryConfirmed);
text=strjoin(fig.UserData.controls.methodScope.Value,newline);
assert(contains(text,'Kinney') && contains(text,'首段局部常状态') && contains(text,'外缘未知'));
press(fig.UserData.controls.confirm);press(fig.UserData.controls.run);
state=fig.UserData.workflowState;
assert(strcmp(state.error.identifier,'RocketAeroHeat:WorkflowIncompleteSurface'),fig.UserData.status.Text);
assert(state.currentStep==4 && isempty(state.result) && isempty(fig.UserData.webUI));
assert(~state.steps(5).available && ~isfield(fig.UserData.controls,'export'));
firstSolve=state.solveId;
press(fig.UserData.workflowButtons(3));control=fig.UserData.controls.sidewall;
control.Value='boundaryLayerIntegral';control.ValueChangedFcn(control,[]);drawnow;
state=fig.UserData.workflowState;
assert(isempty(state.result) && ~state.boundaryConfirmed && ~state.steps(5).available);
press(fig.UserData.controls.confirm);press(fig.UserData.controls.run);
state=fig.UserData.workflowState;
assert(state.solveId>firstSolve && strcmp(state.error.identifier,'RocketAeroHeat:WorkflowIncompleteSurface'));
assert(state.currentStep==4 && isempty(state.result) && ~state.steps(5).available && isempty(fig.UserData.webUI));
summary=struct('status','passed','steps',5,'faceCount',size(state.mesh.F,1), ...
    'methodsRecomputed',2,'partialResultsPublished',0,'singleHost',true);
fid=fopen(fullfile(out,'summary.json'),'w','n','UTF-8');fprintf(fid,'%s\n',jsonencode(summary));fclose(fid);
fprintf('两种侧壁方法缺项时均留在计算页，未发布云图。\n');
end
function press(control)
control.ButtonPushedFcn(control,[]);drawnow;
end
