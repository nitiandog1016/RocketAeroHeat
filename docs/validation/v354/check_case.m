function check_case(out,action)
% 模块一：在普通五步侧栏选择同一模型、尾迹扰动初态与固定物理时刻。
evalin('base','main');app=evalin('base','appFigure');guard=onCleanup(@()delete(app)); %#ok<NASGU>
app.Visible='off';press(app.UserData.controls.load);press(app.UserData.controls.confirm);
press(app.UserData.controls.generate);press(app.UserData.workflowButtons(3));
selector=app.UserData.controls.launcherState;assert(strcmp(selector.Value,'steady'));
selector.Value='wake';selector.ValueChangedFcn(selector,[]);press(app.UserData.controls.confirm);
timer=tic;
if strcmp(action,'archive'),press(app.UserData.controls.readLauncher);else,press(app.UserData.controls.run);end
s=app.UserData.workflowState;assert(s.currentStep==5 && isfield(s.result,'launcherSummary'),app.UserData.status.Text);
ui=app.UserData.webUI;renderTimer=tic;
while isempty(ui.html.UserData.rendered) && toc(renderTimer)<40,drawnow;pause(.1);end
gl=ui.html.UserData.rendered;assert(gl.webglError==0 && gl.visiblePixelSamples>0 && gl.faces==9216);
waitSeconds=toc(timer);
% 模块二：核对时刻、源帧、给定壁温、面元拾取和两个实际截面的导出。
r=s.result;theme=r.launcherSummary;d=buildAerothermalViewerData(r);
assert(isequal(r.hybridCondition,s.freestream) && strcmp(r.launcherAction,action));
assert(abs(d.resultCondition.elapsedPhysicalTimeS-.6)<1e-7 && contains(d.caseNote,'0.600 s'));
assert(~contains(d.caseNote,'稳态轴对称') && ~isfield(r,'launcherAngularHeat'));
assert(strcmp(d.evidenceIdentity.sourceFrame.thermalSHA256,theme.thermalFrameSHA256));
assert(d.wallCoverage.fullSurfaceComputed && abs(d.wallCoverage.wholeSurfaceHeatPowerW+theme.outwardW)<1e-12);
v=ui.evidenceView;ui.evidenceTab.Parent.SelectedTab=ui.evidenceTab;
for key={'pressureGauge','heatFlux'}
    label=key{1};v.selectField(label);
    for id=[1 1537 2017 6913 9121 9169]
        v.inspectFace(id,[]);snapshot=v.getState();assert(snapshot.selectedFaceID==id);
    end
    v.setSectionX(.7);t=tic;
    while toc(t)<35
        drawnow;cut=v.getSection();
        if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,label) && abs(cut.x-.7)<1e-10,break;end
        pause(.1);
    end
    assert(strcmp(cut.fieldKey,label) && abs(cut.x-.7)<1e-10);
    file=fullfile(out,[label '_section.csv']);v.exportSectionTo(file);
    tab=readtable(file,'Delimiter',',','NumHeaderLines',0,'ReadVariableNames',true,'Encoding','UTF-8');
    f=d.fields(strcmp({d.fields.key},label));
    assert(height(tab)>0 && max(abs(tab.value-f.values(tab.face_id)))<1e-11);
end
v.exportWallPairTo(fullfile(out,'wall_pair.csv'));exportAerothermalWeb(r,fullfile(out,'workbench'));
% 模块三：切回稳态必须清除上一时刻结果。
press(app.UserData.workflowButtons(3));selector=app.UserData.controls.launcherState;selector.Value='steady';
selector.ValueChangedFcn(selector,[]);assert(isempty(app.UserData.workflowState.result));
proof=struct('mode',action,'resultWaitSeconds',waitSeconds,'faces',gl.faces,'visiblePixelSamples',gl.visiblePixelSamples, ...
    'sourceFlow',theme.sourceFrameSHA256,'sourceHeat',theme.thermalFrameSHA256,'wallSHA256',theme.wallSHA256, ...
    'timeS',d.resultCondition.elapsedPhysicalTimeS,'outwardW',theme.outwardW,'switchInvalidated',true);
writelines(string(jsonencode(proof,'PrettyPrint',true)),fullfile(out,'ui.json'),'Encoding','UTF-8');
end

function press(control)
% 控件回调与手动点击走相同的工作台状态转换。
control.ButtonPushedFcn(control,[]);drawnow;
end
