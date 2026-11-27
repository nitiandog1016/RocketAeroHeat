function check_workbench(out,action)
% 模块一：从普通五步侧栏选择启动瞬态，区分读取整件归档与点击计算。
evalin('base','main');app=evalin('base','appFigure');
guard=onCleanup(@()delete(app)); %#ok<NASGU>
app.Visible='off';press(app.UserData.controls.load);press(app.UserData.controls.confirm);
press(app.UserData.controls.generate);press(app.UserData.workflowButtons(3));
sel=app.UserData.controls.launcherState;assert(strcmp(sel.Value,'steady'));
sel.Value='wake';sel.ValueChangedFcn(sel,[]);press(app.UserData.controls.confirm);
timer=tic;
if strcmp(action,'archive'),press(app.UserData.controls.readLauncher);else,press(app.UserData.controls.run);end
state=app.UserData.workflowState;
assert(state.currentStep==5 && isfield(state.result,'launcherWakeHeatGrid'),app.UserData.status.Text);
ui=app.UserData.webUI;
while isempty(ui.html.UserData.rendered) && toc(timer)<230,drawnow;pause(.1);end
gl=ui.html.UserData.rendered;assert(gl.webglError==0 && gl.visiblePixelSamples>0 && gl.faces==9216);
waitSeconds=toc(timer);result=state.result;assert(strcmp(result.launcherAction,action));
% 模块二：逐面核对数值、单位和色标。点选六区，在同一实际截面导出三种热流读数。
ui.evidenceTab.Parent.SelectedTab=ui.evidenceTab;view=ui.evidenceView;
d=view.getData();r=result.launcherWakeHeatGrid;lo=result.hybridWall.heatFluxWm2;hi=r.heatFluxWm2;
assert(abs(d.resultCondition.elapsedPhysicalTimeS-.6)<1e-7);
assert(strcmp(d.wakeGridSources.flowSHA256,result.evidenceIdentity.sourceFrame.flowSHA256));
both=[lo;hi];bounds=[min(both) max(both)];c=d.comparisonColorRanges;
assert(isequal(c.heatFlux,bounds) && isequal(c.wakeFineHeat,bounds));
assert(isequal(c.wakeHeatDifference,max(abs(hi-lo))*[-1 1]));
assert(d.wallCoverage.fullSurfaceComputed && abs(d.wallCoverage.jointComputedAreaM2-r.areaM2)<1e-10);
keys={'heatFlux','wakeFineHeat','wakeHeatDifference'};vals={lo,hi,hi-lo};rows=zeros(1,3);
for i=1:numel(keys)
    key=keys{i};view.selectField(key);field=d.fields(strcmp({d.fields.key},key));
    assert(isequal(field.values(:),vals{i}) && all(field.status==1) && strcmp(field.unit,'W/m²'));
    for id=[1 1537 2017 6913 9102 9169]
        view.inspectFace(id,[]);current=view.getState();
        assert(current.selectedFaceID==id && strcmp(current.fieldKey,key));
    end
    view.setSectionX(.7);poll=tic;
    while toc(poll)<35
        drawnow;cut=view.getSection();
        if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,key) && abs(cut.x-.7)<1e-10,break;end
        pause(.1);
    end
    assert(strcmp(cut.fieldKey,key) && abs(cut.x-.7)<1e-10);
    path=fullfile(out,[key '_section.csv']);view.exportSectionTo(path);
    tab=readtable(path,'Encoding','UTF-8');rows(i)=height(tab);
    assert(rows(i)>0 && max(abs(tab.value-field.values(tab.face_id)))<1e-13);
end
view.exportWallPairTo(fullfile(out,'wall_pair.csv'));exportAerothermalWeb(result,fullfile(out,'workbench'));
job=app.UserData.caseReport.runDirectory;
% 模块三：切回稳态清除瞬态复核，并实际读取既有稳态场。
press(app.UserData.workflowButtons(3));sel=app.UserData.controls.launcherState;
sel.Value='steady';sel.ValueChangedFcn(sel,[]);assert(isempty(app.UserData.workflowState.result));
press(app.UserData.controls.confirm);press(app.UserData.controls.readLauncher);
steady=app.UserData.workflowState;assert(steady.currentStep==5 && isfield(steady.result,'launcherAngularHeat'));
assert(~isfield(steady.result,'launcherWakeHeatGrid') && ~isfield(steady.freestream,'elapsedPhysicalTimeS'));
receipt=struct('action',action,'resultWaitSeconds',waitSeconds,'faces',gl.faces,'visiblePixelSamples',gl.visiblePixelSamples, ...
    'regionsProbedPerField',6,'fields',{keys},'sectionX',.7,'sectionRows',rows,'restoredSteady',true, ...
    'coarseWallSHA256',r.coarseWallSHA256,'fineWallSHA256',r.fineWallSHA256,'flowSHA256',r.flowSHA256, ...
    'fineHeatOrigin','archived comparison','runDirectory',job);
writelines(string(jsonencode(receipt,'PrettyPrint',true)),fullfile(out,'ui.json'),'Encoding','UTF-8');
end

function press(control)
% 此处执行用户控件的实际回调，并处理图形事件。
control.ButtonPushedFcn(control,[]);drawnow;
end
