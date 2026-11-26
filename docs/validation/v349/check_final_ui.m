function check_final_ui(out)
% 模块一：执行现有状态回归，重放归档以检查普通结果页的逐面接收条件。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(root,fullfile(root,'tests'));
testAeroHeatWorkflowState([],fullfile(root,'outputs','v107_validation','state01','workflow_state_regression.mat'));
% 模块二：从五步侧栏读取低速试件已有完整归档，确认结果和导出可用。
fig=rocketAeroHeatWorkbench();fig.Visible='off';
closeFig=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.file.Value=fullfile(root,'outputs','v341_validation','ui_final','input.stl');
c.scale.Value=1;press(c.load);press(fig.UserData.controls.confirm);
c=fig.UserData.controls;choose(c.meshBackend,'shortSupport65');press(fig.UserData.controls.generate);
press(fig.UserData.controls.nextBoundary);choose(fig.UserData.controls.shortHeatBoundary,'film');
press(fig.UserData.controls.confirm);press(fig.UserData.controls.readShort);
s=fig.UserData.workflowState;r=fig.UserData.caseReport;
assert(s.currentStep==5 && r.fullSurfaceComputed && r.jointComputedFaceCount==32480,fig.UserData.status.Text);
assert(strcmp(fig.UserData.controls.export.Enable,'on') && ~isempty(fig.UserData.webUI));
% 模块三：在同一窗口开始普通工况。用隔离目录中的单面缺值样本替代求解调用。
% 返回的正负数仅用于软件测试，不代表模型物理量。临时测试函数只在本次进程生效。
press(fig.UserData.workflowButtons(1));c=fig.UserData.controls;
cfg=createCaseConfig();c.file.Value=cfg.geometry.stlFile;c.scale.Value=cfg.geometry.scale;
press(c.load);press(fig.UserData.controls.pca);press(fig.UserData.controls.confirm);
press(fig.UserData.controls.generate);press(fig.UserData.workflowButtons(3));press(fig.UserData.controls.confirm);
stubDir=fullfile(out,'fixture');mkdir(stubDir);
fixture=fullfile(stubDir,'runAeroHeatSingleCase.m');
copyfile(fullfile(fileparts(mfilename('fullpath')),'partial_result_fixture.txt'),fixture);
original=pwd;restore=onCleanup(@()cd(original)); %#ok<NASGU>
cd(stubDir);clear runAeroHeatSingleCase;
assert(strcmp(which('runAeroHeatSingleCase'),fixture));
press(fig.UserData.controls.run);s=fig.UserData.workflowState;
assert(strcmp(s.error.identifier,'RocketAeroHeat:WorkflowIncompleteSurface'),fig.UserData.status.Text);
assert(s.boundaryConfirmed && ~s.busy && s.currentStep==4 && ~s.steps(5).available);
assert(isempty(s.result) && isempty(fig.UserData.caseReport) && isempty(fig.UserData.webUI));
assert(strcmp(fig.UserData.workflowButtons(5).Enable,'off') && ~isfield(fig.UserData.controls,'export'));
assert(strcmp(fig.UserData.resultHost.Visible,'off') && contains(s.error.message,'1 个缺少'));
proof=struct('stateRegressionPassed',true,'archiveFaces',r.jointComputedFaceCount, ...
    'archiveExportEnabled',true,'partialTestFaces',size(s.mesh.F,1),'partialMissingFaces',1, ...
    'visibleStep',s.currentStep,'resultPageEnabled',s.steps(5).available,'partialResultStored',~isempty(s.result));
writelines(string(jsonencode(proof,'PrettyPrint',true)),fullfile(out,'ui.json'),'Encoding','UTF-8');
end
function press(control)
control.ButtonPushedFcn(control,[]);drawnow;
end
function choose(control,value)
control.Value=value;control.ValueChangedFcn(control,[]);drawnow;
end
