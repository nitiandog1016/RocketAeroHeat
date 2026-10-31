function testKinneyWorkbench(out)
%TESTKINNEYWORKBENCH 验证统一五步骤中的方法选择、当前计算与失效流程。
% 模块一：真实STL、朝向、重网格在同一宿主完成，不创建新结果窗口。
% 模块二：切换侧壁公式只修改边界配置；显式计算才产生新场量。
% 模块三：切回默认方法必须清除旧结果，重算后与原积分对照逐值一致。
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
assert(state.currentStep==5 && strcmp(state.status,'solved'),fig.UserData.status.Text);
assert(isfield(state.result.heating,'kinneyAcreage'));
assert(ancestor(fig.UserData.webUI.html,'figure')==fig);
selected=state.result.heating;
data=buildAerothermalViewerData(state.result);
assert(any(strcmp({data.fields.key},'heatingReferenceLength')) && any(strcmp({data.fields.key},'baselineHeatFlux')));
wall=data.fields(strcmp({data.fields.key},'wallSolved'));assert(all(wall.status==0));
report=fig.UserData.caseReport;
press(fig.UserData.workflowButtons(3));control=fig.UserData.controls.sidewall;
control.Value='boundaryLayerIntegral';control.ValueChangedFcn(control,[]);drawnow;
state=fig.UserData.workflowState;
assert(isempty(state.result) && ~state.boundaryConfirmed && ~state.steps(5).available);
press(fig.UserData.controls.confirm);press(fig.UserData.controls.run);
state=fig.UserData.workflowState;
assert(strcmp(state.status,'solved') && ~isfield(state.result.heating,'kinneyAcreage'));
assert(isequaln(state.result.heating.heatFlux,selected.kinneyAcreage.baselineHeatFlux));
summary=struct('status','passed','steps',5,'faceCount',size(state.mesh.F,1), ...
    'selectedSidewallFaces',nnz(selected.kinneyAcreage.acceptedFace), ...
    'heatAreaFraction',report.heatComputedAreaFraction,'methodsRecomputed',2, ...
    'singleHost',true,'nativeWebGLRenderClaimed',false, ...
    'scope','隐藏宿主验证真实控件回调/计算/数据，实际三维显示另用浏览器核验');
fid=fopen(fullfile(out,'summary.json'),'w','n','UTF-8');fprintf(fid,'%s\n',jsonencode(summary));fclose(fid);
nativeHeat=selected.heatFlux;defaultHeat=state.result.heating.heatFlux;
save(fullfile(out,'workbench_evidence.mat'),'summary','report','nativeHeat','defaultHeat','-v7');
fprintf('Kinney workbench passed: 5 steps, %d faces, 2 explicit method solves, stale-result invalidation and exact default replay in one host.\n',summary.faceCount);
end
function press(control)
control.ButtonPushedFcn(control,[]);drawnow;
end
