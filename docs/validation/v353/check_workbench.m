function check_workbench(out,action)
% 模块一：通过正常五步侧栏加载米制模型并确认低速风洞边界。
evalin('base','main');app=evalin('base','appFigure');
cleanup=onCleanup(@()delete(app)); %#ok<NASGU>
app.Visible='off';press(app.UserData.controls.load);press(app.UserData.controls.confirm);
press(app.UserData.controls.generate);press(app.UserData.workflowButtons(3));press(app.UserData.controls.confirm);
timer=tic;
if strcmp(action,'archive'),press(app.UserData.controls.readLauncher);else,press(app.UserData.controls.run);end
state=app.UserData.workflowState;
assert(state.currentStep==5 && isfield(state.result,'launcherAngularHeat'),app.UserData.status.Text);
ui=app.UserData.webUI;
while isempty(ui.html.UserData.rendered) && toc(timer)<900,drawnow;pause(.1);end
gl=ui.html.UserData.rendered;
assert(isstruct(gl) && gl.webglError==0 && gl.visiblePixelSamples>0 && gl.faces==9216);
seconds=toc(timer);

% 模块二：选择新增的全周复核场，拾取六区，在实际交线上逐行回查导出值。
ui.evidenceTab.Parent.SelectedTab=ui.evidenceTab;ev=ui.evidenceView;
key='launcherAngularHeat';ev.selectField(key);
data=ev.getData();field=data.fields(strcmp({data.fields.key},key));r=state.result.launcherAngularHeat;
assert(isequal(field.values(:),r.heatFluxWm2(:)) && strcmp(field.unit,'W/m²') && all(field.status==1));
assert(isequal(data.comparisonColorRanges.heatFlux,data.comparisonColorRanges.(key)));
definition=jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')),'..','..','..','cases','simple_launcher','geometry.json')));
for region=definition.regions'
    id=region.firstFace;ev.inspectFace(id,[]);s=ev.getState();
    assert(s.selectedFaceID==id && strcmp(s.fieldKey,key));
end
ev.setSectionX(.7);timer=tic;
while toc(timer)<40
    drawnow;cut=ev.getSection();
    if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,key) && abs(cut.x-.7)<1e-10,break;end
    pause(.1);
end
assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,key) && abs(cut.x-.7)<1e-10);
path=fullfile(out,'section.csv');ev.exportSectionTo(path);
tab=readtable(path,'Encoding','UTF-8');assert(height(tab)>0 && max(abs(tab.value-field.values(tab.face_id)))<1e-13);
ev.exportWallPairTo(fullfile(out,'wall_pair.csv'));
exportAerothermalWeb(state.result,fullfile(out,'workbench'));
receipt=struct('mode',action,'faces',gl.faces,'visiblePixelSamples',gl.visiblePixelSamples,'resultWaitSeconds',seconds, ...
    'regionsProbed',numel(definition.regions),'sectionRows',height(tab),'sectionX',cut.x,'angularCSV',r.csvSHA256, ...
    'flowSHA256',r.flowSHA256,'angularFieldSHA256',r.fieldSHA256,'runDirectory',app.UserData.caseReport.runDirectory);
writelines(string(jsonencode(receipt,'PrettyPrint',true)),fullfile(out,'ui.json'),'Encoding','UTF-8');
end

function press(control)
% 模块三：执行同一控件的实际回调并提交当前界面事件。
control.ButtonPushedFcn(control,[]);drawnow;
end
