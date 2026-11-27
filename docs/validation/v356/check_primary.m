function check_primary(out,action,expectCoarse)
% 模块一：通过普通侧栏载入简化箭体并执行指定的计算或归档操作。
evalin('base','main');app=evalin('base','appFigure');guard=onCleanup(@()delete(app)); %#ok<NASGU>
app.Visible='off';press(app.UserData.controls.load);press(app.UserData.controls.confirm);
press(app.UserData.controls.generate);press(app.UserData.workflowButtons(3));
selector=app.UserData.controls.launcherState;selector.Value='wake';selector.ValueChangedFcn(selector,[]);
press(app.UserData.controls.confirm);timer=tic;
if strcmp(action,'solve'),press(app.UserData.controls.run);else,press(app.UserData.controls.readLauncher);end
state=app.UserData.workflowState;
assert(state.currentStep==5 && isfield(state.result,'launcherSummary'),app.UserData.status.Text);
ui=app.UserData.webUI;
while isempty(ui.html.UserData.rendered) && toc(timer)<290,drawnow;pause(.1);end
gl=ui.html.UserData.rendered;assert(gl.webglError==0 && gl.visiblePixelSamples>0 && gl.faces==9216);
seconds=toc(timer);r=state.result;s=r.launcherSummary;
assert(strcmp(r.launcherAction,action) && strcmp(s.primaryHeatGrid,'fine'));
assert(isfield(r,'launcherWakeCoarse')==expectCoarse);
job=app.UserData.caseReport.runDirectory;
% 模块二：工作台主热流、首层温度及来源身份必须逐面匹配本次终帧。
if strcmp(action,'solve')
    source=fullfile(job,'heat','disturbed_half_2','wall.csv');
    assert(isfile(fullfile(job,'completed.json')) && ...
        strcmp(cm011FileSHA256(fullfile(job,'heat','disturbed_half_2','thermal.npz')),s.thermalFrameSHA256));
else
    source=fullfile(job,'wall','fine_wall.csv');
end
assert(strcmp(cm011FileSHA256(source),s.wallSHA256));w=readtable(source);
ui.evidenceTab.Parent.SelectedTab=ui.evidenceTab;view=ui.evidenceView;d=view.getData();
assert(strcmp(d.defaultField,'heatFlux') && strcmp(d.evidenceIdentity.sourceFrame.wallSHA256,s.wallSHA256));
assert(strcmp(d.evidenceIdentity.sourceFrame.thermalSHA256,s.thermalFrameSHA256));
assert(d.wallCoverage.fullSurfaceComputed && abs(d.wallCoverage.wholeSurfaceHeatPowerW+s.outwardW)<1e-12);
keys={'heatFlux','pressure','nearWallCellTemperature','wallInput'};
values={w.heatFluxWm2,w.pressurePa/1000,w.nearWallK,repmat(300.001,height(w),1)};
if expectCoarse
    coarse=r.launcherWakeCoarse.heatFluxWm2;
    keys=[keys {'wakeCoarseHeat','wakeHeatDifference'}];values=[values {coarse,w.heatFluxWm2-coarse}];
    c=d.comparisonColorRanges;both=[w.heatFluxWm2;coarse];
    assert(isequal(c.heatFlux,[min(both) max(both)]) && isequal(c.wakeCoarseHeat,c.heatFlux));
    assert(isequal(c.wakeHeatDifference,[-1 1]*max(abs(w.heatFluxWm2-coarse))));
end
probes=[1 1537 2017 6913 9102 9169];rows=zeros(size(keys));
for k=1:numel(keys)
    key=keys{k};view.selectField(key);field=d.fields(strcmp({d.fields.key},key));
    assert(max(abs(field.values(:)-values{k}))<1e-12);
    if strcmp(key,'heatFlux'),assert(contains(field.label,'细网') && contains(field.methodNames{1},'细热网'));end
    for face=probes
        view.inspectFace(face,[]);v=view.getState();assert(v.selectedFaceID==face && strcmp(v.fieldKey,key));
    end
    view.setSectionX(.7);poll=tic;
    while toc(poll)<35
        drawnow;cut=view.getSection();
        if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,key) && abs(cut.x-.7)<1e-10,break;end
        pause(.1);
    end
    assert(strcmp(cut.fieldKey,key) && abs(cut.x-.7)<1e-10);
    path=fullfile(out,[key '_section.csv']);view.exportSectionTo(path);
    tab=readtable(path,'Delimiter',',','Encoding','UTF-8');rows(k)=height(tab);
    assert(rows(k)>0 && max(abs(tab.value-values{k}(tab.face_id)))<1e-12);
end
% 模块三：导出实际配对面表和可直接打开的离线页，保存时间、网格和来源回执。
view.exportWallPairTo(fullfile(out,'wall_pair.csv'));exportAerothermalWeb(r,fullfile(out,'workbench'));
receipt=struct('action',action,'resultWaitSeconds',seconds,'visiblePixelSamples',gl.visiblePixelSamples, ...
    'faces',gl.faces,'probedFaceIDs',probes,'fields',{keys},'sectionX',.7,'sectionRows',rows, ...
    'wallSHA256',s.wallSHA256,'thermalFrameSHA256',s.thermalFrameSHA256,'sourceFrameSHA256',s.sourceFrameSHA256, ...
    'heatMeshSHA256',s.heatMeshSHA256,'outwardW',s.outwardW,'regionW',s.regionW,'runDirectory',job);
writelines(string(jsonencode(receipt,'PrettyPrint',true)),fullfile(out,'ui.json'),'Encoding','UTF-8');
end

function press(button)
% 调用按钮本身的回调并处理界面事件。
button.ButtonPushedFcn(button,[]);drawnow;
end
